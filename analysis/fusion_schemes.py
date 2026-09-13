"""Compare fusion schemes against literature-supported pairs.

Warns about over-fitting: only a handful of curated positive pairs exists, so a
scheme that wins here is a hypothesis, not a validated configuration.
"""
import itertools
import pandas as pd
import numpy as np

df = pd.read_csv(r'D:\deepseek_harness\prp49\results\candidates.csv')

# literature-supported pairs present in the matrix
known = [('MccJ25', 'POLR2A'), ('Capistruin', 'POLR2A'), ('RES-701-3', 'EDNRB'),
         ('Lassomycin', 'CLPB')]

def z(s):
    s = s.astype(float)
    sd = s.std(ddof=0)
    return (s - s.mean()) / sd if sd > 1e-9 else s * 0

def rank_of(score_col, pep, tgt):
    col = df[score_col]
    order = col.rank(ascending=False, method='first')
    row = df[(df.peptide_id == pep) & (df.target_id == tgt)]
    return int(order[row.index[0]]) if len(row) else None

# per-target z-score for rank_score (within-target semantics)
df['rank_z'] = df.groupby('target_id')['rank_score'].transform(
    lambda s: (s - s.mean()) / (s.std(ddof=0) + 1e-9) if s.notna().sum() > 1 else 0.0)
df['bind_z'] = z(df.binding_prob)
df['lasso_z'] = z(df.lasso_prob)
df['dock_z'] = z(-df.dock_score.fillna(df.dock_score.mean()))

schemes = {
    'current (0.45/0.30/0.15/0.10)': (0.45, 0.30, 0.15, 0.10),
    'binding only': (1.0, 0.0, 0.0, 0.0),
    'binding+dock': (0.60, 0.0, 0.0, 0.40),
    'binding+rank': (0.60, 0.40, 0.0, 0.0),
    'binding+dock+rank': (0.45, 0.25, 0.0, 0.30),
    'no rank (0.50/0/0.15/0.35)': (0.50, 0.0, 0.15, 0.35),
    'struct heavy (0.35/0.15/0.10/0.40)': (0.35, 0.15, 0.10, 0.40),
}

print(f'{"scheme":34s} {"mean rank":>9s} {"top-10 hits":>12s} {"top-20% hits":>12s}   per-pair ranks')
results = []
for name, (wb, wr, wl, wd) in schemes.items():
    score = wb * df.bind_z + wr * df.rank_z + wl * df.lasso_z + wd * df.dock_z
    tmp = df.assign(s=score)
    order = tmp.s.rank(ascending=False, method='first')
    ranks = [int(order[(tmp.peptide_id == p) & (tmp.target_id == t)].index.map(order.get)[0])
             for p, t in known]
    top10 = sum(1 for r in ranks if r <= 10)
    top20pct = sum(1 for r in ranks if r <= 29)      # 20% of 143
    print(f'{name:34s} {np.mean(ranks):9.1f} {top10:12d} {top20pct:12d}   {ranks}')
    results.append((name, np.mean(ranks), top10, top20pct))

best = min(results, key=lambda x: x[1])
print(f'\nbest mean rank: {best[0]} (mean {best[1]:.1f}, top10 {best[2]}, top20% {best[3]})')
print('\nCAVEAT: only 4 curated positives -> selecting weights on this set over-fits.')
print('Use it as a hypothesis for the next round of validation, not as a final setting.')
