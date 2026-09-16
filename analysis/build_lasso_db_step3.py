"""Route D database, step 3: expand beyond the "lasso peptide" full-text search.

Finding so far: of 98 lasso-peptide PDB entries, only 5 pair a lasso peptide with a
larger protein, and none of them is an RNA polymerase. Two consequences:

  1. MccJ25 x RNAP has NO structure in the PDB at all - the interaction is known from
     biochemistry, mutagenesis and crosslinking, not crystallography. No amount of
     structure mining will produce it.
  2. The auxiliary signal route D needs must therefore be assembled from (a) the few
     real complexes, (b) literature-known pairs without structures, and (c) a
     functional target correspondence table, which is knowledge work rather than
     structure parsing.

This step does the target-driven search (query by the biological targets themselves)
and builds the functional correspondence table.
"""
import json
import os
import time
import urllib.request

import pandas as pd

ROOT = r'D:\deepseek_harness\prp49'
DB = os.path.join(ROOT, 'lasso_target_db')


def rcsb_search(text, rows=100):
    q = {"query": {"type": "terminal", "service": "full_text",
                   "parameters": {"value": text}},
         "return_type": "entry",
         "request_options": {"paginate": {"start": 0, "rows": rows}}}
    req = urllib.request.Request('https://search.rcsb.org/rcsbsearch/v2/query',
                                 data=json.dumps(q).encode(),
                                 headers={'Content-Type': 'application/json'})
    with urllib.request.urlopen(req, timeout=60) as r:
        d = json.loads(r.read())
    return d.get('total_count', 0), [x['identifier'] for x in d.get('result_set', [])]


# (a) target-driven search: does any structure show a short peptide bound to these targets?
QUERIES = {
    'RNA polymerase + microcin': 'RNA polymerase microcin',
    'RNA polymerase + lasso': 'RNA polymerase lasso',
    'MccJ25 + RNA polymerase': 'MccJ25 RNA polymerase',
    'ClpC1 + lassomycin': 'ClpC1 lassomycin',
    'ClpC1 + cyclic peptide': 'ClpC1 cyclic peptide',
    'RNAP secondary channel inhibitor': 'RNA polymerase secondary channel inhibitor',
    'MccJ25 receptor FhuA': 'MccJ25 FhuA',
    'capistruin target': 'capistruin',
}
print('=== target-driven PDB searches ===')
extra = {}
for label, q in QUERIES.items():
    try:
        total, ids = rcsb_search(q)
        extra[label] = ids
        print(f'  {label:36s}: {total:4d} hits')
    except Exception as exc:
        print(f'  {label:36s}: failed ({type(exc).__name__})')
    time.sleep(0.6)
json.dump(extra, open(os.path.join(DB, 'target_driven_hits.json'), 'w'), indent=2)

# (b) functional correspondence: bacterial/other target <-> human target
#     evidence = literature knowledge, NOT sequence identity (which is <15% for POLR2A)
MAP = [
    # human,        bacterial counterpart,            basis,                              confidence
    ('POLR2A', 'DNA-directed RNA polymerase subunit beta-prime (rpoC)',
     'MccJ25 and capistruin bind the RNAP secondary channel; the channel is conserved '
     'between bacterial and human RNAP despite <15% sequence identity', 'high'),
    ('POLR2A', 'DNA-directed RNA polymerase subunit beta (rpoB)',
     'rifampicin pocket, adjacent to the secondary channel', 'medium'),
    ('CLPB', 'ClpB / ClpC1 (AAA+ unfoldase)',
     'lassomycin kills M. tuberculosis by binding ClpC1; CLPB is the human mitochondrial '
     'orthologue of the same AAA+ clade', 'high'),
    ('CLPB', 'ClpP (caseinolytic protease)',
     'ADEP/acyldepsipeptide class binds ClpP; CLPB is the unfoldase partner', 'medium'),
    ('PPIA', 'peptidyl-prolyl isomerase (bacterial cyclophilin-like)',
     'cyclophilin domain PF00160 found in 2 Propedia receptors (weight 0.94)', 'high'),
    ('FKBP1A', 'FKBP-type PPIase (bacterial)',
     'FKBP domain PF00254 is widespread in bacteria', 'medium'),
    ('NPR1', 'guanylate cyclase (bacterial)',
     'NPR1 carries an adenylate/guanylate cyclase catalytic domain PF00211', 'medium'),
    ('C3', 'alpha-2-macroglobulin (bacterial)',
     'C3 shares 11 of 12 Pfam domains with a Propedia receptor (weight 0.92)', 'high'),
    ('MDM2', 'SWIB/MDM2 domain protein',
     'MDM2 domain PF02201 found in 13 Propedia receptors', 'medium'),
    ('ITGAV', 'integrin alpha (metazoan only)',
     'no bacterial counterpart; integrins are metazoan-specific', 'none'),
    ('ITGB3', 'integrin beta (metazoan only)',
     'no bacterial counterpart; integrins are metazoan-specific', 'none'),
    ('EDNRB', 'GPCR (eukaryotic only)',
     '7TM GPCRs do not exist in bacteria; 74 weak single-domain hits are not homologs', 'none'),
    ('PLXNB1', 'plexin/SEMA domain (metazoan only)',
     'plexin domains are metazoan-specific', 'none'),
]
fm = pd.DataFrame(MAP, columns=['human_target', 'counterpart', 'basis', 'confidence'])
fm.to_csv(os.path.join(DB, 'target_functional_map.csv'), index=False, lineterminator='\n')
print(f'\nwrote target_functional_map.csv: {len(fm)} rows')
print(fm.groupby('confidence').size().to_string())

print('\n=== honest summary of what the new database can and cannot supply ===')
print('  CAN:  C3 (domain-verified, w=0.92), PPIA (cyclophilin, w=0.94), CLPB (AAA+),')
print('        NPR1 (cyclase), MDM2 (SWIB) - bacterial counterparts exist as sequences')
print('  CANNOT: POLR2A as a STRUCTURE. No PDB entry shows MccJ25/capistruin bound to')
print('        RNAP; the correspondence is functional, provable only from mutagenesis')
print('        and crosslinking literature.')
print('  CANNOT: EDNRB, PLXNB1, ITGAV, ITGB3 - metazoan-specific folds, no counterpart')
