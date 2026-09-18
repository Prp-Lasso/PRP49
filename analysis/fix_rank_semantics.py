"""Fix the rank semantics: `rank_target` was a GLOBAL rank over 650 rows computed from
within-target-standardised scores - not a rank inside the target. Reporting it as a
within-target rank produced nonsense like "110/13".

Correct definitions now:
  rank_global        : rank over all 650 rows using globally standardised scores
  rank_targetglobal  : rank over all 650 rows using within-target-standardised scores
  rank_in_target     : rank of this peptide AMONG the peptides tested on ITS target
                       (this is the number that answers "which peptide for this target")
"""
import os

import numpy as np
import pandas as pd

ROOT = r'D:\deepseek_harness\prp49'
p = os.path.join(ROOT, 'results', 'candidates_screen50_within_target.csv')
d = pd.read_csv(p)
pd.set_option('display.width', 220)

d['rank_in_target'] = (d.groupby('target_id').composite_target
                        .rank(ascending=False, method='min'))
d['n_peptides_for_target'] = d.groupby('target_id').peptide_id.transform('count')

print('=== sanity: is rank_in_target actually bounded by the per-target peptide count? ===')
bad = d[d.rank_in_target > d.n_peptides_for_target]
print(f'  rows violating the bound: {len(bad)}')
print(f'  e.g. POLR2A has {int(d[d.target_id=="POLR2A"].n_peptides_for_target.iloc[0])} peptides, '
      f'ranks run 1..{int(d[d.target_id=="POLR2A"].rank_in_target.max())}')
assert len(bad) == 0, 'rank_in_target exceeds the per-target count - still wrong'

print('\n=== positive controls, now with correct within-target ranks ===')
for pep, tgt, lvl in [('RES-701-3', 'EDNRB', 'A: 9KDF cryo-EM'), ('RES-701-1', 'EDNRB', 'B: IC50 10 nM'),
                      ('Anantin', 'NPR1', 'B: Kd 0.6 uM'), ('MccJ25', 'POLR2A', 'C: RpoC homology'),
                      ('Capistruin', 'POLR2A', 'C: RpoC homology')]:
    s = d[(d.peptide_id == pep) & (d.target_id == tgt)]
    if len(s):
        r = s.iloc[0]
        print(f'  {pep:11s} x {tgt:8s} [{lvl:20s}] rank {int(r.rank_in_target)}'
              f'/{int(r.n_peptides_for_target)}  composite {r.composite_target:+.3f}')

print('\n=== per-peptide: how often is it the BEST peptide for a target? ===')
prof = (d[d.rank_in_target == 1].groupby('peptide_id').size()
        .rename('n_best_for_a_target').sort_values(ascending=False))
print(prof.to_string() if len(prof) else '  (none)')

print('\n=== the deliverable: best peptide per target ===')
best = d[d.rank_in_target == 1].sort_values('composite_target', ascending=False)
KNOWN = {'POLR2A', 'CLPB', 'EDNRB', 'NPR1', 'C3', 'MDM2', 'ITGAV', 'ITGB3', 'PPIA',
         'FKBP1A', 'PLXNB1'}
best = best.assign(is_known_target=best.target_id.isin(KNOWN))
print(best[['target_id', 'peptide_id', 'composite_target', 'binding_prob', 'dock_score',
            'is_known_target']].round(3).to_string(index=False))

d.sort_values(['target_id', 'rank_in_target']).to_csv(p, index=False, lineterminator='\n')
print(f'\nupdated {p} with rank_in_target')
