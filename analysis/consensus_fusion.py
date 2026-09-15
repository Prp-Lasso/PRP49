"""Dual-structure consensus fusion + re-evaluation.

Why: the same (peptide, target) pair scores very differently depending on the
receptor structure used - Pearson +0.195 between the two versions, and the top-1
peptide agrees in only 1/9 targets. The shortlist's composite is 60% docking, so
the previous ranking was an artefact of which structures happened to be used.

Consensus rules (from docs/PRP49_当前任务与交接.md appendix A):
  * both versions available -> take the WORSE of the two (max, since more negative
    is better) => a conservative score that cannot be inflated by a lucky structure
  * only one version        -> use it, but flag the row as single-source
  * either version > 0      -> flag unreliable (steric clash / bad box)
  * per-target rank agreement between the two versions -> consistent / conflict
"""
import json
import os

import numpy as np
import pandas as pd

ROOT = r'D:\deepseek_harness\prp49'
HARD_POS = [('RES-701-3', 'EDNRB'), ('RES-701-1', 'EDNRB'), ('Anantin', 'NPR1')]
HARD_NEG = [('MccJ25', 'ITGAV'), ('MccJ25', 'ITGB3')]
MISSING_PENALTY = -0.5
W_BIND, W_DOCK = 0.40, 0.60

af = pd.read_csv(os.path.join(ROOT, 'docking_af', 'af_scores.csv'))
ex = pd.read_csv(os.path.join(ROOT, 'docking', 'dock_scores_long.csv'))
cand = pd.read_csv(os.path.join(ROOT, 'results', 'candidates_grp.csv'))
print(f'af {len(af)} | exp {len(ex)} | candidates {len(cand)}')

# --- replicate the pipeline's heterodimer mapping: ITGAVB3 -> ITGAV and ITGB3 ---
dimer = ex[ex.target == 'ITGAVB3'].copy()
ex = ex[ex.target != 'ITGAVB3'].copy()
ex.loc[ex.index, 'target'] = ex['target']
dimer_a, dimer_b = dimer.copy(), dimer.copy()
dimer_a['target'] = 'ITGAV'
dimer_b['target'] = 'ITGB3'
ex = pd.concat([ex, dimer_a, dimer_b], ignore_index=True)
ex = ex[ex.target != 'CTRL_9KDF']
print(f'exp after heterodimer split: {len(ex)} rows, targets {sorted(ex.target.unique())}')

# ---------------------------------------------------------------- consensus
rows = []
keys = set(zip(cand.peptide_id, cand.target_id))
for pep, tgt in sorted(keys):
    a = af[(af.peptide == pep) & (af.target == tgt)].af_score
    e = ex[(ex.peptide == pep) & (ex.target == tgt)].dock_score
    av = float(a.iloc[0]) if len(a) else None
    ev = float(e.iloc[0]) if len(e) else None
    if av is not None and ev is not None:
        cons, src = max(av, ev), 'both'
    elif ev is not None:
        cons, src = ev, 'exp_only'
    elif av is not None:
        cons, src = av, 'af_only'
    else:
        cons, src = None, 'none'
    rows.append(dict(peptide_id=pep, target_id=tgt, af_score=av, exp_score=ev,
                     dock_consensus=cons, dock_source=src,
                     delta=(av - ev) if (av is not None and ev is not None) else None))
dock = pd.DataFrame(rows)

# reliability flag: any positive score means clash / bad box
dock['unreliable'] = dock.apply(
    lambda r: bool((r.af_score is not None and r.af_score > 0) or
                   (r.exp_score is not None and r.exp_score > 0)), axis=1)

# per-target rank agreement between the two versions
def target_ranks(df, col):
    out = {}
    for tgt, g in df.dropna(subset=[col]).groupby('target_id'):
        order = g.sort_values(col, ascending=True)          # more negative = better
        for i, (_, r) in enumerate(order.iterrows(), start=1):
            out[(r.peptide_id, r.target_id)] = i
    return out

ra, re_ = target_ranks(dock, 'af_score'), target_ranks(dock, 'exp_score')
def agreement(r):
    ka, ke = (r.peptide_id, r.target_id), (r.peptide_id, r.target_id)
    if ka in ra and ke in re_:
        d = abs(ra[ka] - re_[ke])
        return 'consistent' if d <= 2 else 'conflict'
    return 'single_source'
dock['struct_agreement'] = dock.apply(agreement, axis=1)

m = cand.merge(dock, on=['peptide_id', 'target_id'], how='left')
print(f'\nconsensus coverage: {m.dock_consensus.notna().mean()*100:.0f}%  '
      f'(both {int((m.dock_source=="both").sum())}, exp-only {int((m.dock_source=="exp_only").sum())}, '
      f'none {int((m.dock_source=="none").sum())})')
print(f'structure agreement: {m.struct_agreement.value_counts().to_dict()}')
print(f'unreliable (positive score in either version): {int(m.unreliable.sum())}')

# ---------------------------------------------------------------- refusion
def z(x):
    x = np.asarray(x, dtype=float)
    sd = np.nanstd(x)
    return (x - np.nanmean(x)) / sd if sd > 1e-9 else np.zeros_like(x)

m['composite_old'] = m['composite']
m['rank_old'] = m['rank']
bind_z = z(m.binding_prob.values)
has = m.dock_consensus.notna().values
dock_z = z(-m.dock_consensus.fillna(m.dock_consensus.mean()).values)
m['composite_consensus'] = W_BIND * bind_z + W_DOCK * np.where(has, dock_z, MISSING_PENALTY)
m['rank_consensus'] = m['composite_consensus'].rank(ascending=False).astype(int)
n = len(m)

# ---------------------------------------------------------------- evaluation
report = {'n_rows': n, 'protocol': 'dual-structure consensus (worse-of-two dock)',
          'fusion': {'w_binding': W_BIND, 'w_dock': W_DOCK, 'missing_penalty': MISSING_PENALTY},
          'coverage': {'both': int((m.dock_source == 'both').sum()),
                       'exp_only': int((m.dock_source == 'exp_only').sum()),
                       'none': int((m.dock_source == 'none').sum())},
          'positives': {}, 'negatives': {}}

print('\n' + '=' * 78)
print(f'{"pair":26s} {"old rank":>9s} {"new rank":>9s} {"af":>7s} {"exp":>7s} {"cons":>7s}  agreement')
print('=' * 78)
hp_old = hp_new = hn_old = hn_new = 0
for pep, tgt in HARD_POS + HARD_NEG:
    s = m[(m.peptide_id == pep) & (m.target_id == tgt)]
    if s.empty:
        continue
    r = s.iloc[0]
    kind = 'positives' if (pep, tgt) in HARD_POS else 'negatives'
    report[kind][f'{pep}x{tgt}'] = dict(
        rank_old=int(r.rank_old), rank_new=int(r.rank_consensus),
        pct_new=round(100 * r.rank_consensus / n, 1),
        af=None if pd.isna(r.af_score) else round(float(r.af_score), 2),
        exp=None if pd.isna(r.exp_score) else round(float(r.exp_score), 2),
        consensus=None if pd.isna(r.dock_consensus) else round(float(r.dock_consensus), 2),
        agreement=r.struct_agreement, unreliable=bool(r.unreliable))
    f = lambda v: '  n/a ' if pd.isna(v) else f'{v:6.2f}'
    print(f'{pep+" x "+tgt:26s} {int(r.rank_old):9d} {int(r.rank_consensus):9d} '
          f'{f(r.af_score)} {f(r.exp_score)} {f(r.dock_consensus)}  {r.struct_agreement}')
    if kind == 'positives':
        hp_old += r.rank_old <= 0.2 * n
        hp_new += r.rank_consensus <= 0.2 * n
    else:
        hn_old += r.rank_old > 0.5 * n
        hn_new += r.rank_consensus > 0.5 * n

report['scorecard'] = {'positives_top20_old': f'{hp_old}/3', 'positives_top20_new': f'{hp_new}/3',
                       'negatives_controlled_old': f'{hn_old}/2', 'negatives_controlled_new': f'{hn_new}/2'}
print('=' * 78)
print(f"positives in top-20%:  old {hp_old}/3   ->  new {hp_new}/3")
print(f"negatives controlled:  old {hn_old}/2   ->  new {hn_new}/2")

# ---------------------------------------------------------------- outputs
keep = ['peptide_id', 'target_id', 'composite_old', 'rank_old', 'composite_consensus', 'rank_consensus',
        'binding_prob', 'lasso_prob', 'af_score', 'exp_score', 'dock_consensus', 'delta',
        'dock_source', 'struct_agreement', 'unreliable']
out = m[keep].sort_values('rank_consensus')
out_csv = os.path.join(ROOT, 'results', 'candidates_consensus.csv')
out.round(4).to_csv(out_csv, index=False, lineterminator='\n')
json.dump(report, open(os.path.join(ROOT, 'results', 'consensus_evaluation.json'), 'w'),
          indent=2, ensure_ascii=False)
print(f'\nwrote {out_csv}')
print(f'wrote results/consensus_evaluation.json')

print('\ntop 12 by consensus:')
print(out.head(12)[['rank_consensus', 'peptide_id', 'target_id', 'composite_consensus',
                    'binding_prob', 'af_score', 'exp_score', 'dock_consensus',
                    'struct_agreement']].round(3).to_string(index=False))
