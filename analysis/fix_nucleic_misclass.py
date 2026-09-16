"""Fix mistake #42: nucleic acid chains were classified as peptides.

The entity filter used only length (<=80 aa -> "peptide"), but 6N60/6N61/6N62 also
contain short DNA/RNA chains, which sailed through and entered the auxiliary set as
"lasso peptides". The fix uses the deposited entity TYPE ('polypeptide(L)') rather
than length, then rebuilds the auxiliary set and the verification inputs.

Detection aid: a sequence made only of A/C/G/T/U with no other residue letters is
nucleic acid, not protein.
"""
import json
import os
import re
import urllib.request

import pandas as pd

ROOT = r'D:\deepseek_harness\prp49'
DB = os.path.join(ROOT, 'lasso_target_db')
IDS = ['6N60', '6N61', '6N62', '8IBO', '8IBP', '4CU4', '9KDF']
GQL = """{
  entries(entry_ids: [%s]) {
    rcsb_id
    polymer_entities {
      entity_poly { pdbx_seq_one_letter_code_can type }
      rcsb_polymer_entity { pdbx_description }
      entity_src_gen { pdbx_gene_src_scientific_name }
    }
  }
}"""
q = GQL % ','.join(f'"{i}"' for i in IDS)
req = urllib.request.Request('https://data.rcsb.org/graphql',
                             data=json.dumps({'query': q}).encode(),
                             headers={'Content-Type': 'application/json'})
with urllib.request.urlopen(req, timeout=90) as r:
    ents = json.loads(r.read())['data']['entries']

NUC = re.compile(r'^[ACGTU]+$')
rows = []
for e in ents:
    for pe in e.get('polymer_entities') or []:
        poly = (pe or {}).get('entity_poly') or {}
        s = (poly.get('pdbx_seq_one_letter_code_can') or '').replace('\n', '')
        etype = poly.get('type') or ''
        if not s:
            continue
        raw = (pe or {}).get('entity_src_gen')
        if isinstance(raw, list):
            org = (raw[0] or {}).get('pdbx_gene_src_scientific_name') if raw else None
        elif isinstance(raw, dict):
            org = raw.get('pdbx_gene_src_scientific_name')
        else:
            org = None
        rows.append(dict(pdb_id=e['rcsb_id'], entity_type=etype, seq=s, length=len(s),
                         description=((pe or {}).get('rcsb_polymer_entity') or {})
                         .get('pdbx_description'),
                         organism=org, looks_nucleic=bool(NUC.match(s))))

df = pd.DataFrame(rows)
df.to_csv(os.path.join(DB, 'rnap_lasso_entities.csv'), index=False, lineterminator='\n')
print('entity types found:')
print(df.groupby(['entity_type', 'looks_nucleic']).size().to_string())

prot = df[df.entity_type.str.startswith('polypeptide')].copy()
print(f'\npolypeptide entities: {len(prot)} (was {len(df)} raw)')
print(f'dropped as nucleic acid: {int(df.looks_nucleic.sum())}')

peps = prot[prot.length <= 80]
print(f'\nreal peptide chains (<=80 aa, polypeptide): {len(peps)}')
for _, r in peps.iterrows():
    print(f"  {r.pdb_id}  {r.length:3d} aa  {r.seq[:44]:44s} {str(r.description)[:40]}")
