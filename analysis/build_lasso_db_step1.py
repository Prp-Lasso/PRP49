"""Route D database, step 1: harvest lasso-peptide structures from the PDB.

Propedia has no RNA polymerase receptor, which is why route D could not cover the
MccJ25/Capistruin -> POLR2A correspondence. The PDB does contain lasso peptides in
complex with their (mostly bacterial) targets, so it is the right source to build the
missing auxiliary pairs.

Queries run:
  * full-text "lasso peptide"
  * per-name searches for the well-characterised lasso peptides
Output: lasso_target_db/pdb_lasso_entries.csv
"""
import json
import os
import time
import urllib.request

ROOT = r'D:\deepseek_harness\prp49'
DB = os.path.join(ROOT, 'lasso_target_db')
os.makedirs(DB, exist_ok=True)

KNOWN = ['microcin J25', 'MccJ25', 'capistruin', 'lassomycin', 'lariatin', 'klebsazolicin',
         'siamycin', 'anantin', 'propeptin', 'RES-701', 'sphingopyxin', 'ubonodin',
         'chaxapeptin', 'caulosegnin', 'astexin', 'pseudomycoidin', 'lassopeptide',
         'lasso peptide', 'Lasso peptide']


def rcsb_search(query_text, rows=200):
    q = {"query": {"type": "terminal", "service": "full_text",
                   "parameters": {"value": query_text}},
         "return_type": "entry",
         "request_options": {"paginate": {"start": 0, "rows": rows}}}
    req = urllib.request.Request('https://search.rcsb.org/rcsbsearch/v2/query',
                                 data=json.dumps(q).encode(),
                                 headers={'Content-Type': 'application/json'})
    with urllib.request.urlopen(req, timeout=60) as r:
        d = json.loads(r.read())
    return d.get('total_count', 0), [x['identifier'] for x in d.get('result_set', [])]


found = {}
cache = os.path.join(DB, 'search_hits.json')
if os.path.exists(cache):
    found = {k: set(v) for k, v in json.load(open(cache)).items()}
    print(f'loaded {len(found)} cached search hits')
else:
    for term in KNOWN:
        try:
            total, ids = rcsb_search(term)
            for i in ids:
                found.setdefault(i, set()).add(term)
            print(f'{term:20s}: {total:4d} hits ({len(ids)} fetched)')
        except Exception as exc:
            print(f'{term:20s}: FAILED {type(exc).__name__} (rate limit?) - retrying once')
            time.sleep(3)
            try:
                total, ids = rcsb_search(term)
                for i in ids:
                    found.setdefault(i, set()).add(term)
                print(f'{term:20s}: {total:4d} hits on retry')
            except Exception as exc2:
                print(f'{term:20s}: gave up ({type(exc2).__name__})')
        time.sleep(0.6)
    json.dump({k: sorted(v) for k, v in found.items()}, open(cache, 'w'), indent=2)

print(f'\nunique PDB entries: {len(found)}')
json.dump({k: sorted(v) for k, v in found.items()},
          open(os.path.join(DB, 'search_hits.json'), 'w'), indent=2)

# --- pull metadata via GraphQL (batched) ---
ids = sorted(found)
GQL = """{
  entries(entry_ids: [%s]) {
    rcsb_id
    struct { title }
    rcsb_entry_info { resolution_combined experimental_method deposited_polymer_entity_instance_count }
    polymer_entities {
      entity_poly { pdbx_seq_one_letter_code_can type }
      rcsb_polymer_entity { pdbx_description }
      entity_src_gen { pdbx_gene_src_scientific_name }
    }
  }
}"""


def gql(batch):
    q = GQL % ','.join(f'"{i}"' for i in batch)
    req = urllib.request.Request('https://data.rcsb.org/graphql',
                                 data=json.dumps({'query': q}).encode(),
                                 headers={'Content-Type': 'application/json'})
    with urllib.request.urlopen(req, timeout=90) as r:
        return json.loads(r.read()).get('data', {}).get('entries', []) or []

rows = []
for i in range(0, len(ids), 40):
    batch = ids[i:i + 40]
    try:
        ents = gql(batch)
    except Exception as exc:
        print(f'graphql batch {i} failed: {exc}')
        continue
    for e in ents:
        if not e:
            continue
        ents_poly = e.get('polymer_entities') or []
        seqs = []
        for pe in ents_poly:
            poly = (pe or {}).get('entity_poly') or {}
            s = (poly.get('pdbx_seq_one_letter_code_can') or '').replace('\n', '')
            desc = ((pe or {}).get('rcsb_polymer_entity') or {}).get('pdbx_description')
            src_raw = (pe or {}).get('entity_src_gen')
            if isinstance(src_raw, list):          # GraphQL returns a list here
                src = (src_raw[0] or {}).get('pdbx_gene_src_scientific_name') if src_raw else None
            elif isinstance(src_raw, dict):
                src = src_raw.get('pdbx_gene_src_scientific_name')
            else:
                src = None
            if s:
                seqs.append(dict(len=len(s), seq=s, desc=desc, organism=src))
        info = e.get('rcsb_entry_info') or {}
        rows.append(dict(pdb_id=e.get('rcsb_id'),
                         title=((e.get('struct') or {}).get('title') or '')[:150],
                         method=info.get('experimental_method'),
                         resolution=(info.get('resolution_combined') or [None])[0],
                         n_chains=info.get('deposited_polymer_entity_instance_count'),
                         n_entities=len(seqs),
                         chains=json.dumps(seqs)[:4000]))
    print(f'  metadata {min(i+40, len(ids))}/{len(ids)}', flush=True)
    time.sleep(0.3)

import pandas as pd
df = pd.DataFrame(rows)
out = os.path.join(DB, 'pdb_lasso_entries.csv')
df.to_csv(out, index=False, lineterminator='\n')
print(f'\nwrote {out}: {len(df)} entries')
if len(df):
    print('methods:', df.method.value_counts().to_dict())
    print('\nsample titles:')
    for _, r in df.head(8).iterrows():
        print(f"  {r.pdb_id}: {r.title[:90]}")
