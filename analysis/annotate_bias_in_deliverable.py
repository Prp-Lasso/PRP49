"""Write the screening corrections and peptide-bias annotation into the deliverable.

Two things go into the shortlist itself (not just into prose):
  1. the within-target ranking (already computed) - the ranking that answers
     "which peptide for THIS target"
  2. an explicit peptide-bias flag: how many distinct targets a peptide ranks first for.
     A peptide topping many unrelated targets is more likely a peptide-level constant
     than a genuine multi-target binder, and the deliverable must say so per row.
"""
import os
import sys
sys.path.insert(0, r'D:\deepseek_harness\prp49\hpc_deploy\ops')
from opslib import connect, run, put_file

import numpy as np
import pandas as pd

ROOT = r'D:\deepseek_harness\prp49'
p = os.path.join(ROOT, 'results', 'candidates_screen50_within_target.csv')
d = pd.read_csv(p)

KNOWN = {'POLR2A', 'CLPB', 'EDNRB', 'NPR1', 'C3', 'MDM2', 'ITGAV', 'ITGB3', 'PPIA',
         'FKBP1A', 'PLXNB1'}
d['is_known_target'] = d.target_id.isin(KNOWN)

# ---- peptide bias: in how many targets does this peptide come first? ----
n_best = d[d.rank_in_target <= 1].groupby('peptide_id').size()
d['peptide_n_targets_ranked_1st'] = d.peptide_id.map(n_best).fillna(0).astype(int)

def bias_level(n):
    """Thresholds chosen from the observed distribution, and stated so they can be argued with."""
    if n >= 5:
        return 'high'        # tops 5+ unrelated targets -> likely an identity constant
    if n >= 2:
        return 'medium'
    return 'low'

d['peptide_bias'] = d.peptide_n_targets_ranked_1st.map(bias_level)

# also record whether the docking term was usable for this row
d['dock_used'] = d.get('dock_ok', pd.Series(True, index=d.index)).fillna(False)
d['dock_clash'] = d.dock_score.notna() & (d.dock_score > 0)

print('=== peptide bias summary ===')
prof = (d.drop_duplicates('peptide_id')
          [['peptide_id', 'peptide_n_targets_ranked_1st', 'peptide_bias']]
          .sort_values('peptide_n_targets_ranked_1st', ascending=False))
print(prof.to_string(index=False))

print('\n=== how many rows are affected ===')
print(d.peptide_bias.value_counts().to_string())

print('\n=== the 39 drug-panel targets: best peptide each, WITH bias flags ===')
panel = d[(~d.is_known_target) & (d.rank_in_target <= 2)].sort_values(
    ['target_id', 'rank_in_target'])
print(panel[['target_id', 'rank_in_target', 'peptide_id', 'composite_target',
             'binding_prob', 'dock_score', 'peptide_bias']].round(3).to_string(index=False))

out = os.path.join(ROOT, 'results', 'candidates_screen50_within_target.csv')
d.sort_values(['target_id', 'rank_in_target']).to_csv(out, index=False, lineterminator='\n')
print(f'\nupdated {out} with peptide_n_targets_ranked_1st / peptide_bias / dock_used / dock_clash')

# ---- ship it: local releases + cluster ----
import shutil
for rel in ['v2_2026-09-18_full', 'FINAL_2026-09-18_v2_full']:
    shutil.copy2(out, os.path.join(ROOT, 'releases', rel, 'candidates_screen50.csv'))
print('copied into both releases')

cli = connect()
st, h, _ = run(cli, 'echo $HOME')
put_file(cli, out, f'{h.strip()}/LassoPep/results/candidates_screen50_annotated.csv', progress=False)
cli.close()

print('\n=== P0-3 fix: why batch 8 + max_len 512 should fit now ===')
print('  failed run : batch 16 x max_len 1024, TWO forward passes (pos+neg) -> >40 GB')
print('  fixed run  : batch  8 x max_len  512, same two passes')
print('  activation memory scales ~linearly in batch and ~quadratically in length for attention')
print('  -> roughly 8/16 x (512/1024)^2 = 1/8 of the attention footprint: ~5 GB, safe on a 40 GB A100')
