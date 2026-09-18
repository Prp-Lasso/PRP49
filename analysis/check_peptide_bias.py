"""Peptide-bias check on the corrected shortlist.

The released model's known weakness is that scores are dominated by peptide identity
(62.8% of variance on seen targets, 99.7% on unseen). So a peptide topping many targets
is a warning sign, not a result: check each peptide's profile ACROSS targets.
"""
import os

import numpy as np
import pandas as pd

ROOT = r'D:\deepseek_harness\prp49'
d = pd.read_csv(os.path.join(ROOT, 'results', 'candidates_screen50_within_target.csv'))
pd.set_option('display.width', 220)

print('=== per-peptide profile across the 50 targets ===')
prof = d.groupby('peptide_id').agg(
    mean_composite=('composite_target', 'mean'),
    std_composite=('composite_target', 'std'),
    n_targets=('target_id', 'nunique'),
    n_top1=('rank_target', lambda s: int((s <= 1).sum())),
    n_top5=('rank_target', lambda s: int((s <= 5).sum())),
    mean_binding=('binding_logit', 'mean'),
    mean_dock=('dock_score', 'mean'),
).round(3).sort_values('mean_composite', ascending=False)
print(prof.to_string())

print('\n=== interpretation ===')
hi_mean = prof[prof.mean_composite > 0.8]
if len(hi_mean):
    print('  peptides with a high MEAN composite across all targets (bias suspects):')
    for p, r in hi_mean.iterrows():
        print(f'    {p:16s} mean {r.mean_composite:.2f} | top-1 in {int(r.n_top1):2d} targets | '
              f'top-5 in {int(r.n_top5):2d} | mean binding logit {r.mean_binding:.2f}')
    print('  -> a peptide that ranks first for many UNRELATED targets is more likely a')
    print('     peptide-level constant than a genuine multi-target binder.')

print('\n=== how much of the within-target ranking is explained by binding vs docking? ===')
sub = d.dropna(subset=['dock_z_target'])
if len(sub) > 10:
    r_b = np.corrcoef(sub.binding_z_target, sub.composite_target)[0, 1]
    r_d = np.corrcoef(sub.dock_z_target, sub.composite_target)[0, 1]
    r_bd = np.corrcoef(sub.binding_z_target, sub.dock_z_target)[0, 1]
    print(f'  corr(composite, binding_z) {r_b:+.3f}')
    print(f'  corr(composite, dock_z)    {r_d:+.3f}   (expected ~0.60 by construction)')
    print(f'  corr(binding_z, dock_z)    {r_bd:+.3f}   (near 0 => the two signals are complementary)')

print('\n=== positive controls under the corrected ranking ===')
for pep, tgt, lvl in [('RES-701-3', 'EDNRB', 'A: 9KDF cryo-EM'), ('RES-701-1', 'EDNRB', 'B: IC50 10 nM'),
                      ('Anantin', 'NPR1', 'B: Kd 0.6 uM'), ('MccJ25', 'POLR2A', 'C: RpoC homology'),
                      ('Capistruin', 'POLR2A', 'C: RpoC homology')]:
    s = d[(d.peptide_id == pep) & (d.target_id == tgt)]
    if len(s):
        r = s.iloc[0]
        n_in_tgt = int((d.target_id == tgt).sum())
        print(f'  {pep:11s} x {tgt:8s} [{lvl:20s}] within-target rank '
              f'{int(r.rank_target):2d}/{n_in_tgt}  composite {r.composite_target:+.3f}')
