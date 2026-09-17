"""Full analysis of the real-complex verification.

Key diagnostic: is the score driven by PEPTIDE identity or by TARGET compatibility?
A model that learned transferable binding chemistry should vary strongly across
targets for a fixed peptide. A model that merely recognises "this is a known lasso
peptide" gives one peptide high scores everywhere.
"""
import os

import numpy as np
import pandas as pd

ROOT = r'D:\deepseek_harness\prp49'
m = pd.read_csv(os.path.join(ROOT, 'results', 'verify_real_complexes.csv')).set_index('id')
truth = pd.read_csv(os.path.join(ROOT, 'mvp_cpu', 'verify_truth.csv'))

pd.set_option('display.width', 200)
print('=== full matrix (logit) ===')
print(m.round(3).to_string())

print('\n=== variance decomposition ===')
pep_mean = m.mean(axis=1)
tgt_mean = m.mean(axis=0)
grand = m.values.mean()
ss_pep = len(m.columns) * ((pep_mean - grand) ** 2).sum()
ss_tgt = len(m.index) * ((tgt_mean - grand) ** 2).sum()
ss_tot = ((m.values - grand) ** 2).sum()
print(f'  between-peptide SS : {ss_pep:9.1f}  ({100*ss_pep/ss_tot:5.1f}% of total)')
print(f'  between-target  SS : {ss_tgt:9.1f}  ({100*ss_tgt/ss_tot:5.1f}% of total)')
print(f'  residual (interaction) : {ss_tot-ss_pep-ss_tgt:9.1f}  ({100*(ss_tot-ss_pep-ss_tgt)/ss_tot:5.1f}%)')
print(f'  per-peptide spread across targets (std):')
for k, v in m.std(axis=1).items():
    print(f'    {k:20s} {v:6.2f}   (range {m.loc[k].min():7.2f} .. {m.loc[k].max():7.2f})')

print('\n=== rank of each TRUE complex within its target column ===')
print('(rank 1 = the model considers this the best peptide for that target)')
rows = []
for _, t in truth.iterrows():
    col = t.target_id
    if col not in m.columns:
        cand = [c for c in m.columns if c.startswith(col[:18])]
        if not cand:
            continue
        col = cand[0]
    s = m[col].sort_values(ascending=False)
    if t.peptide_id not in s.index:
        continue
    r = int(list(s.index).index(t.peptide_id)) + 1
    rows.append(dict(peptide=t.peptide_id, target=col, pdb=t.pdb_id,
                     human=t.human_counterpart, score=round(float(s[t.peptide_id]), 2),
                     rank_in_column=f'{r}/{len(s)}', column_max=round(float(s.max()), 2),
                     column_min=round(float(s.min()), 2)))
res = pd.DataFrame(rows)
print(res.to_string(index=False))

print('\n=== the decisive question: does the peptide beat its own decoys? ===')
for _, t in truth.iterrows():
    pep = t.peptide_id
    col = t.target_id if t.target_id in m.columns else None
    if col is None or pep not in m.index:
        continue
    row = m.loc[pep]
    true_score = row[col]
    others = row.drop(col)
    # how many targets score higher than the TRUE partner?
    n_higher = int((others > true_score).sum())
    print(f'  {pep:20s} true={col[:26]:26s} {true_score:7.2f} | '
          f'targets scoring higher: {n_higher}/{len(others)} | '
          f'row max {row.max():7.2f} on {row.idxmax()[:24]}')

print('\n=== interpretation aids ===')
print('  MccJ25 / capistruin are known lasso peptides; lassomycin is too, but in the')
print('  CLPB direction. If those first two score highly on nearly EVERY target, the')
print('  score reflects peptide familiarity, not target-specific binding chemistry.')
known = [i for i in m.index if 'Microcin' in i or 'Capistruin' in i]
for k in known:
    row = m.loc[k]
    print(f'  {k:20s} min {row.min():7.2f} median {row.median():7.2f} max {row.max():7.2f} '
          f'-> {int((row > 5).sum())}/{len(row)} targets above logit 5')
