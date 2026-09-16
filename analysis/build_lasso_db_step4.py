"""Route D database, step 4: pull the RNAP-lasso complexes the full-text search missed.

6N60/6N61/6N62 were found only through the target-driven query ("RNA polymerase microcin"),
not through any "lasso peptide" search - the titles never use the word lasso. They are
precisely the structures route D was missing: MccJ25 and capistruin bound to a bacterial
RNA polymerase, i.e. the bacterial counterpart of human POLR2A.

Also pulls 5UI6/5UI7 (other RNAP-microcin hits) for completeness.
"""
import json
import os
import time
import urllib.request

import pandas as pd

ROOT = r'D:\deepseek_harness\prp49'
DB = os.path.join(ROOT, 'lasso_target_db')

IDS = ['6N60', '6N61', '6N62', '5UI6', '5UI7', '4CU4', '9KDF', '8IBO', '8IBP', '3HBL', '2B5U', '1JCH']
GQL = """{
  entries(entry_ids: [%s]) {
    rcsb_id
    struct { title }
    rcsb_entry_info { resolution_combined experimental_method }
    polymer_entities {
      entity_poly { pdbx_seq_one_letter_code_can type }
      rcsb_polymer_entity { pdbx_description }
      entity_src_gen { pdbx_gene_src_scientific_name }
    }
  }
}"""


def gql(ids):
    q = GQL % ','.join(f'"{i}"' for i in ids)
    req = urllib.request.Request('https://data.rcsb.org/graphql',
                                 data=json.dumps({'query': q}).encode(),
                                 headers={'Content-Type': 'application/json'})
    with urllib.request.urlopen(req, timeout=90) as r:
        return json.loads(r.read()).get('data', {}).get('entries', []) or []


ents = gql(IDS)
print(f'fetched {len(ents)} entries\n')

rows = []
for e in ents:
    if not e:
        continue
    info = e.get('rcsb_entry_info') or {}
    res = (info.get('resolution_combined') or [None])[0]
    print(f"=== {e['rcsb_id']} ({info.get('experimental_method')}, "
          f"{res if res is None else round(res,2)} A) ===")
    print(f"    {(e.get('struct') or {}).get('title','')[:110]}")
    polys = []
    for pe in e.get('polymer_entities') or []:
        poly = (pe or {}).get('entity_poly') or {}
        s = (poly.get('pdbx_seq_one_letter_code_can') or '').replace('\n', '')
        if not s:
            continue
        desc = ((pe or {}).get('rcsb_polymer_entity') or {}).get('pdbx_description') or ''
        raw = (pe or {}).get('entity_src_gen')
        if isinstance(raw, list):
            org = (raw[0] or {}).get('pdbx_gene_src_scientific_name') if raw else None
        elif isinstance(raw, dict):
            org = raw.get('pdbx_gene_src_scientific_name')
        else:
            org = None
        polys.append((len(s), s, desc, org))
    polys.sort(key=lambda x: x[0])
    for L, s, desc, org in polys:
        tag = 'PEPTIDE' if L <= 80 else 'protein'
        print(f'      [{tag:7s}] {L:5d} aa  {str(desc)[:58]:58s} {str(org)[:28]}')
        rows.append(dict(pdb_id=e['rcsb_id'], title=(e.get('struct') or {}).get('title', ''),
                         method=info.get('experimental_method'), resolution=res,
                         entity_len=L, entity_kind=tag, seq=s,
                         description=desc, organism=org))
    print()

df = pd.DataFrame(rows)
df.to_csv(os.path.join(DB, 'rnap_lasso_entities.csv'), index=False, lineterminator='\n')
print(f'wrote rnap_lasso_entities.csv: {len(df)} entities')

# --- the pairs that matter for route D ---
peps = df[df.entity_kind == 'PEPTIDE']
prots = df[df.entity_kind == 'protein']
pairs = []
for _, p in peps.iterrows():
    for _, t in prots.iterrows():
        if p.pdb_id != t.pdb_id:
            continue
        pairs.append(dict(pdb_id=p.pdb_id, peptide_len=p.entity_len, peptide_seq=p.seq,
                          peptide_desc=p.description, target_len=t.entity_len,
                          target_seq=t.seq, target_desc=t.description,
                          target_organism=t.organism, resolution=p.resolution))
pr = pd.DataFrame(pairs).drop_duplicates(subset=['pdb_id', 'peptide_seq', 'target_seq'])
pr.to_csv(os.path.join(DB, 'rnap_lasso_pairs.csv'), index=False, lineterminator='\n')
print(f'\nwrote rnap_lasso_pairs.csv: {len(pr)} (peptide, target) pairs')
if len(pr):
    print(pr[['pdb_id', 'peptide_len', 'target_len', 'target_desc', 'target_organism']]
          .to_string(index=False))
