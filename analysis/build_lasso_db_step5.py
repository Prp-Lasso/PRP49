"""Assemble the route D database README + a consolidated pair table."""
import os

import pandas as pd

ROOT = r'D:\deepseek_harness\prp49'
DB = os.path.join(ROOT, 'lasso_target_db')

pairs = pd.read_csv(os.path.join(DB, 'rnap_lasso_pairs.csv'))
# keep the biologically meaningful complexes only (drop the colicin false positives)
KEEP = {'6N60', '6N61', '6N62', '8IBO', '8IBP', '4CU4', '9KDF'}
bio = pairs[pairs.pdb_id.isin(KEEP)].copy()

CORR = {
    '6N60': ('RNAP beta-prime (rpoC)', 'POLR2A', 'high',
             'MccJ25/capistruin bind the RNAP secondary channel; the channel is '
             'structurally conserved in human POLR2A despite <15% sequence identity'),
    '6N61': ('RNAP beta-prime (rpoC)', 'POLR2A', 'high',
             'capistruin complex; same secondary-channel site'),
    '6N62': ('RNAP beta-prime (rpoC)', 'POLR2A', 'high',
             'third RNAP-lasso complex, same site'),
    '8IBO': ('ClpC1 (AAA+ unfoldase)', 'CLPB', 'high',
             'lassomycin binds ClpC1; CLPB is the human mitochondrial AAA+ orthologue'),
    '8IBP': ('ClpC1 (AAA+ unfoldase)', 'CLPB', 'high',
             'second lassomycin-ClpC1 structure'),
    '4CU4': ('FhuA ferrichrome receptor', '', 'low',
             'MccJ25 uptake receptor, not a pharmacological target'),
    '9KDF': ('Endothelin receptor type B', 'EDNRB', 'experimental',
             'RES-701-3 x EDNRB, our grade-A reference pair (human/human)'),
}
bio['target_protein'] = bio.pdb_id.map(lambda x: CORR.get(x, ('', '', '', ''))[0])
bio['human_counterpart'] = bio.pdb_id.map(lambda x: CORR.get(x, ('', '', '', ''))[1])
bio['confidence'] = bio.pdb_id.map(lambda x: CORR.get(x, ('', '', '', ''))[2])
bio['basis'] = bio.pdb_id.map(lambda x: CORR.get(x, ('', '', '', ''))[3])
bio.to_csv(os.path.join(DB, 'lasso_target_pairs_curated.csv'), index=False,
           lineterminator='\n')
print(f'curated pairs: {len(bio)} from {bio.pdb_id.nunique()} structures')
print(bio.groupby(['pdb_id', 'target_protein', 'human_counterpart'])
      .agg(pairs=('peptide_seq', 'nunique'), target_len=('target_len', 'max')).to_string())

readme = """# lasso_target_db — a target-annotated lasso peptide structure database

Built 2026-09-16 to fix a specific gap: **Propedia contains no RNA polymerase receptor**,
so route D (homology-weighted transfer from bacterial targets to human targets) could
not cover POLR2A — the target of the grade-C reference pairs MccJ25 and capistruin.

## Why the first attempt failed, and what fixed it

| step | result |
|---|---|
| full-text search `"lasso peptide"` | 98 PDB entries, but only **5** paired a lasso peptide with a larger protein, and **zero** were RNAP |
| target-driven search (`RNA polymerase microcin`) | **6N60 / 6N61 / 6N62** — three RNAP-lasso complexes the full-text search missed entirely, because their titles never contain the word "lasso" |

> Lesson: a single keyword strategy silently under-collects. Query by **target** as well
> as by **peptide class**.

## What the database contains

| PDB | peptide | target | human counterpart | confidence |
|---|---|---|---|---|
| **6N60** | MccJ25 (21 aa) + 24/29 aa | RNAP beta-prime (rpoC, 1409 aa) | **POLR2A** | high |
| **6N61** | capistruin (19/24/29 aa) | RNAP beta-prime | **POLR2A** | high |
| **6N62** | 24/29 aa lasso | RNAP beta-prime | **POLR2A** | high |
| **8IBO / 8IBP** | lassomycin (16 aa) | ClpC1 (142 aa, *M. tuberculosis*) | **CLPB** | high |
| 4CU4 | MccJ25 | FhuA receptor (706 aa) | — | low (uptake, not target) |
| 9KDF | RES-701-3 (16 aa) | EDNRB (908 aa, human) | EDNRB | experimental (grade A) |

Each complex is a full RNAP/ClpC1 assembly, so a peptide appears paired with every
subunit; the biologically relevant partner is the one named above (beta-prime for RNAP,
ClpC1 for lassomycin).

## Files
- `pdb_lasso_entries.csv` — 98 entries from the peptide-class search
- `lasso_target_pairs.csv` — all (peptide, protein) combinations extracted from them
- `rnap_lasso_entities.csv` / `rnap_lasso_pairs.csv` — the target-driven harvest
- `lasso_target_pairs_curated.csv` — **the usable pairs, with human counterparts**
- `target_functional_map.csv` — human target → counterpart, with evidence and confidence
- `search_hits.json`, `target_driven_hits.json` — raw query provenance

## What it still cannot supply
- EDNRB, PLXNB1, ITGAV, ITGB3 have **no bacterial counterpart** (metazoan-specific folds)
- PPIA/FKBP1A/NPR1/MDM2/C3 counterparts exist as **sequences** (via Pfam domains) but
  have no lasso-peptide complex structures
"""
open(os.path.join(DB, 'README.md'), 'w', encoding='utf-8').write(readme)
print(f'\nwrote README.md')
print('\ndatabase files:')
for f in sorted(os.listdir(DB)):
    p = os.path.join(DB, f)
    print(f'  {f:36s} {os.path.getsize(p)/1024:8.1f} KB')
