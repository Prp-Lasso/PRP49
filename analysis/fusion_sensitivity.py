"""Fusion-weight sensitivity analysis - the honest version of "optimising the weights".

Why not just fit the weights: the usable human-target literature pairs number about six
(RES-701-1/3 x EDNRB, Anantin x NPR1, BI-32169 x GCGR, plus the two functional RNAP
correspondences). Six points cannot identify two free parameters without overfitting a
number we would then quote as if it were validated.

What CAN be done, and is more useful for a deliverable: show how the selection responds
to the weighting. A pair that survives every weighting is robust to a choice we cannot
justify from data; a pair that appears only at one weighting is an artefact of that choice.
"""
import os

import numpy as np
import pandas as pd

ROOT = r'D:\deepseek_harness\prp49'
d = pd.read_csv(os.path.join(ROOT, 'results', 'candidates_screen50_audited.csv'))
pd.set_option('display.width', 240)

KNOWN = {'POLR2A', 'CLPB', 'EDNRB', 'NPR1', 'C3', 'MDM2', 'ITGAV', 'ITGB3', 'PPIA',
         'FKBP1A', 'PLXNB1'}
d['is_known_target'] = d.target_id.isin(KNOWN)


def z(s):
    sd = s.std()
    return (s - s.mean()) / sd if sd and sd > 1e-9 else s * 0.0


d['bz'] = d.groupby('target_id').binding_logit.transform(z)
d['dz'] = d[d.dock_score.notna() & (d.dock_score <= 0)].groupby('target_id').dock_score \
    .transform(lambda s: z(-s)).reindex(d.index)

WEIGHTS = [0.4, 0.5, 0.6, 0.7, 0.8]
results = {}
for wd in WEIGHTS:
    wb = 1.0 - wd
    comp = np.where(d.dz.isna(), d.bz, wb * d.bz + wd * d.dz)
    d['c'] = comp
    d['r_in'] = d.groupby('target_id').c.rank(ascending=False, method='min')
    el = d[(~d.is_known_target) & (d.peptide_bias != 'high') &
           d.dock_score.notna() & (d.dock_score <= 0) & (d.r_in <= 3)].copy()
    el['agree'] = ((el.bz > 0) & (el.dz > 0)).astype(int)
    el['sel'] = 0.40 * z(el.c) + 0.30 * z(-el.r_in) + 0.30 * el.agree
    el = el.sort_values('sel', ascending=False)
    picked, cnt = [], {}
    for _, r in el.iterrows():
        if cnt.get(r.target_id, 0) >= 2:
            continue
        cnt[r.target_id] = cnt.get(r.target_id, 0) + 1
        picked.append(f'{r.peptide_id} x {r.target_id}')
        if len(picked) >= 10:
            break
    results[wd] = picked
    print(f'w_dock={wd:.1f}: {len(picked)} selected')

print('\n=== stability across weightings ===')
from collections import Counter
c = Counter(p for v in results.values() for p in v)
n_weights = len(WEIGHTS)
print(f'{"pair":32s} {"appears in":>11s}  verdict')
stable, marginal = [], []
for pair, k in c.most_common():
    tag = 'STABLE (all weightings)' if k == n_weights else (
        'mostly' if k >= n_weights - 1 else 'weighting-dependent')
    (stable if k == n_weights else marginal).append(pair)
    print(f'{pair:32s} {k:>5d}/{n_weights}   {tag}')

print(f'\nstable pairs: {len(stable)}')
for p in stable:
    print('  ', p)
print(f'\nweighting-dependent pairs: {len(marginal)}')
for p in marginal:
    print('  ', p)

print('\n=== the default weighting (0.40/0.60) list, with stability annotation ===')
for p in results[0.6]:
    k = c[p]
    print(f'  {p:32s} in {k}/{n_weights} weightings'
          + ('   <- robust' if k == n_weights else '   <- sensitive'))

out = os.path.join(ROOT, 'results', 'fusion_weight_sensitivity.csv')
pd.DataFrame({'pair': list(c.keys()),
              'n_weightings_selected': [c[k] for k in c.keys()],
              'n_weightings_tested': n_weights}).to_csv(out, index=False, lineterminator='\n')
print(f'\nwrote {out}')
