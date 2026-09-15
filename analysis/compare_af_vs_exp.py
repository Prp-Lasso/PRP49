"""Compare AF-source vs experimental-source docking: magnitude, outliers, ranking impact."""
import os

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = r'D:\deepseek_harness\prp49'
af = pd.read_csv(os.path.join(ROOT, 'docking_af', 'af_scores.csv'))
exp = pd.read_csv(os.path.join(ROOT, 'docking', 'dock_scores_long.csv'))
print(f'AF pairs {len(af)} | experimental pairs {len(exp)}')

# ---- outliers in the AF set ----
pos = af[af.af_score > 0].sort_values('af_score', ascending=False)
print(f'\n=== AF positives (score > 0) : {len(pos)} ===')
if len(pos):
    print(pos.to_string(index=False))
    print('\nper-target outlier counts:')
    print(pos.groupby('target').size().to_string())

# ---- which pair is missing ----
tasks = pd.read_csv(os.path.join(ROOT, 'docking_af', 'tasks_af.csv'))
want = set(zip(tasks.pep, tasks.rec.str.replace('af_', '')))
have = set(zip(af.peptide, af.target))
print(f'\nmissing pair(s): {sorted(want - have)}')

# ---- join with the experimental panel ----
m = af.merge(exp.rename(columns={'dock_score': 'exp_score'}),
             on=['peptide', 'target'], how='inner')
# drop the clash outliers from the comparison (they are box artefacts, not binding)
clean = m[m.af_score <= 0].copy()
clean['delta'] = clean.af_score - clean.exp_score

print(f'\n=== paired comparison: {len(m)} pairs, {clean.shape[0]} after dropping AF clashes ===')
print(f'delta (AF - exp): mean {clean.delta.mean():+.2f}  median {clean.delta.median():+.2f} '
      f'sd {clean.delta.std():.2f}  range [{clean.delta.min():+.2f}, {clean.delta.max():+.2f}]')
print(f'correlation of the two score sets: pearson {clean.af_score.corr(clean.exp_score):+.3f}')

print('\nper target:')
rows = []
for t, g in clean.groupby('target'):
    rho = spearmanr(g.af_score, g.exp_score)[0] if len(g) > 3 else np.nan
    rows.append(dict(target=t, n=len(g), mean_af=g.af_score.mean(), mean_exp=g.exp_score.mean(),
                     mean_delta=g.delta.mean(), spearman=rho))
pt = pd.DataFrame(rows).sort_values('mean_delta')
print(pt.round(3).to_string(index=False))

# ---- does the ranking of peptides within a target change? ----
print('\n=== within-target peptide ranking agreement ===')
same_top = 0
for t, g in clean.groupby('target'):
    if len(g) < 3:
        continue
    best_af = g.loc[g.af_score.idxmin(), 'peptide']
    best_exp = g.loc[g.exp_score.idxmin(), 'peptide']
    agree = best_af == best_exp
    same_top += agree
    print(f'  {t:9s} best by AF: {best_af:14s} | best by exp: {best_exp:14s} '
          f'{"AGREE" if agree else "differ"}')
print(f'top-1 agreement: {same_top} targets')

# ---- literature pairs: where do they land in the AF set? ----
print('\n=== literature-supported pairs under AF structures ===')
KNOWN = [('RES-701-3', 'EDNRB', 'A: 9KDF'), ('RES-701-1', 'EDNRB', 'B: IC50 10 nM'),
         ('MccJ25', 'POLR2A', 'C: RpoC homology'), ('Capistruin', 'POLR2A', 'C: RpoC homology'),
         ('Lassomycin', 'CLPB', 'C: ClpC1 homology'), ('Anantin', 'NPR1', 'B: Kd 0.6 uM')]
for pep, tgt, note in KNOWN:
    g = af[af.target == tgt]
    if g.empty:
        print(f'  {pep} x {tgt}: target not in AF set'); continue
    gx = exp[exp.target == tgt]
    if g[g.peptide == pep].empty:
        af_r = None
    else:
        af_r = int(g.af_score.rank().loc[g[g.peptide == pep].index[0]])
    if gx[gx.peptide == pep].empty:
        exp_r = None
    else:
        exp_r = int(gx.dock_score.rank().loc[gx[gx.peptide == pep].index[0]])
    print(f'  {pep:12s} x {tgt:8s} AF rank {af_r}/{len(g)} | exp rank {exp_r}/{len(gx)}  ({note})')

out = os.path.join(ROOT, 'docking_af', 'af_vs_exp_comparison.csv')
clean.to_csv(out, index=False, lineterminator='\n')
pt.to_csv(os.path.join(ROOT, 'docking_af', 'af_vs_exp_per_target.csv'), index=False, lineterminator='\n')
print(f'\nwrote {out} and af_vs_exp_per_target.csv')
