"""Extended fusion analysis: 11 peptides x 10 targets, 8 known positive pairs."""
import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = r'D:\deepseek_harness\prp49'
dock = pd.read_csv(os.path.join(ROOT, 'docking', 'dock_matrix_full.csv'), index_col=0)
scan = pd.read_csv(os.path.join(ROOT, 'results', 'scan_matrix_ext.csv'), index_col=0)

target_map = {
    'EDNRB': ['EDNRB'], 'PLXNB1': ['PLXNB1'], 'POLR2A': ['POLR2A'], 'PPIA': ['PPIA'],
    'FKBP1A': ['FKBP1A'], 'NPR1': ['NPR1'], 'C3': ['C3'], 'MDM2': ['MDM2'],
    'ITGAVB3': ['ITGAV', 'ITGB3'], 'CLPB': ['CLPB'],
}
scan_cols = {c.split('|')[0]: c for c in scan.columns}
model = pd.DataFrame(index=scan.index)
for tgt, keys in target_map.items():
    cols = [scan_cols[k] for k in keys if k in scan_cols]
    if cols:
        model[tgt] = scan[cols].mean(axis=1)

peps = [p for p in dock.index if p in model.index]
tgts = [c for c in dock.columns if c in model.columns and c != 'CTRL_9KDF']
print(f'peptides: {len(peps)}, targets: {len(tgts)}, pairs: {len(peps)*len(tgts)}')
M = model.loc[peps, tgts]
D = -dock.loc[peps, tgts]

a, b = M.values.flatten(), D.values.flatten()
ok = ~np.isnan(a) & ~np.isnan(b)
rho, p = spearmanr(a[ok], b[ok])
print(f'[Spearman] model vs -dock: rho={rho:.3f} p={p:.3f} (n={ok.sum()})')

known = [
    ('RES-701-3', 'EDNRB'), ('RES-701-1', 'EDNRB'),
    ('PB1m7', 'PLXNB1'),
    ('MccJ25', 'POLR2A'), ('Capistruin', 'POLR2A'), ('Capi-var1', 'POLR2A'),
    ('Lassomycin', 'CLPB'), ('Ubonodin', 'POLR2A'),
]

def z(df):
    return (df - df.mean()) / (df.std(ddof=0) + 1e-9)

def rank_of(sdf, pep, tgt):
    col = sdf[tgt].dropna()
    if pep not in col.index:
        return None
    return list(col.sort_values(ascending=False).index).index(pep) + 1

def top20(sdf, pairs):
    hits = 0
    for pep, tgt in pairs:
        r = rank_of(sdf, pep, tgt)
        col = sdf[tgt].dropna()
        if r and r <= max(1, int(round(0.2 * len(col)))):
            hits += 1
    return hits

print('\n[per-target ranking of 8 known pairs]')
for label, sdf, sign in [('model', z(M), 1), ('dock', z(D), 1)]:
    ranks = [rank_of(sdf, p, t) for p, t in known]
    print(f'  {label:6s} ranks={ranks} mean={np.mean([r for r in ranks if r]):.2f}')

print('\n[fusion sweep] z(model) + w*z(dock)')
best = None
for w in [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0]:
    F = z(M) + w * z(D)
    ranks = [rank_of(F, p, t) for p, t in known]
    rr = [r for r in ranks if r]
    h = top20(F, known)
    print(f'  w={w:<4} mean rank={np.mean(rr):.2f} top20%={h}/{len(known)} ranks={ranks}')
    if best is None or np.mean(rr) < best[1]:
        best = (w, np.mean(rr), h)
print(f'\nbest w={best[0]} (mean rank {best[1]:.2f}, top20% {best[2]}/{len(known)})')
