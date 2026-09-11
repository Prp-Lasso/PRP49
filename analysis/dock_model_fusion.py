"""Fuse model scores (sequence pLM) with docking scores (structure) and evaluate.

Inputs:
  docking/dock_matrix.csv        vina scores, 8 peptides x 10 targets (+ control)
  results/scan_matrix.csv        model logits, 10 peptides x 11 chains (from job_scan)
Outputs:
  - Spearman correlation(model, -dock) overall and per target
  - ranking of known positive pairs under model-only / dock-only / fused scores
  - fusion weight sweep
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = r'D:\deepseek_harness\prp49'
dock = pd.read_csv(os.path.join(ROOT, 'docking', 'dock_matrix.csv'), index_col=0)
scan_path = os.path.join(ROOT, 'results', 'scan_matrix.csv')
if not os.path.exists(scan_path):
    raise SystemExit(f'missing {scan_path} (run job_scan first)')
scan = pd.read_csv(scan_path, index_col=0)

# align target columns: ITGAVB3 = mean(ITGAV, ITGB3)
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

common_peps = [p for p in dock.index if p in model.index]
common_tgts = [c for c in dock.columns if c in model.columns and c != 'CTRL_9KDF']
print(f'common peptides: {len(common_peps)}, targets: {len(common_tgts)}')
M = model.loc[common_peps, common_tgts]
D = -dock.loc[common_peps, common_tgts]     # higher = better binding

# ---- Spearman ----
flat_m = M.values.flatten()
flat_d = D.values.flatten()
ok = ~np.isnan(flat_m) & ~np.isnan(flat_d)
rho, p = spearmanr(flat_m[ok], flat_d[ok])
print(f'\n[Spearman] model vs -dock (all pairs): rho={rho:.3f} p={p:.3f} (n={ok.sum()})')
print('\nper-target Spearman:')
for t in common_tgts:
    a, b = M[t].values, D[t].values
    ok2 = ~np.isnan(a) & ~np.isnan(b)
    if ok2.sum() >= 4:
        r2, p2 = spearmanr(a[ok2], b[ok2])
        print(f'  {t:9s} rho={r2:+.3f} p={p2:.3f}')

# ---- ranking of known positives ----
known = [('RES-701-3', 'EDNRB'), ('MccJ25', 'POLR2A'), ('Capistruin', 'POLR2A'),
         ('Lassomycin', 'CLPB')]
def rank_of(score_df, pep, tgt):
    col = score_df[tgt].dropna()
    if pep not in col.index:
        return None, None
    order = col.sort_values(ascending=False)
    return list(order.index).index(pep) + 1, len(order)

def z(df):
    return (df - df.mean()) / (df.std(ddof=0) + 1e-9)

print('\n[known positive pairs] rank within peptide column (10 targets):')
for label, sdf in [('model', z(M)), ('dock', z(D))]:
    ranks = [rank_of(sdf, p, t) for p, t in known]
    print(f'  {label:6s}: ' + ', '.join(f'{p}@{t}={r[0]}/{r[1]}' for (p, t), r in zip(known, ranks)))

print('\n[fusion sweep] score = z(model) + w * z(dock), rank of known pairs in target column:')
for w in [0.0, 0.25, 0.5, 0.75, 1.0, 1.5]:
    F = z(M) + w * z(D)
    ranks = [rank_of(F, p, t)[0] for p, t in known]
    ranks = [r for r in ranks if r]
    print(f'  w={w:<4} mean rank={np.mean(ranks):.2f}  ranks={ranks}')
