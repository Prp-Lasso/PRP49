"""Select the final ~10 candidate peptide-target pairs for delivery.

The selection must be defensible, so the rules are explicit and applied in order:

  FILTER (hard, each excludes rather than down-weights)
    1. druggable panel target only - the 11 Lasso targets are positive controls, not
       screening output (they are reported separately)
    2. peptide_bias != high - a peptide that ranks first for 5+ unrelated targets is more
       likely a peptide-level constant than a genuine binder (measured: RES-701-1 tops 15)
    3. docking must be usable - a numeric score that is not a clash (>0). Rows without a
       docking score rest on the model alone, whose absolute scores we showed are not
       interpretable, so they are excluded from the shortlist (but kept in the full table)
    4. rank_in_target <= 3 - must be among the top three peptides for that target

  RANK (soft, among survivors)
    z-scored combination of: within-target rank, composite, and signal agreement
    (binding and docking agreeing is stronger than either alone - they correlate only +0.103)

  DIVERSITY
    at most 2 peptides per target, so the list is not one target repeated

Everything excluded is reported with the reason, so the shortlist can be audited.
"""
import os

import numpy as np
import pandas as pd

ROOT = r'D:\deepseek_harness\prp49'
d = pd.read_csv(os.path.join(ROOT, 'results', 'candidates_screen50_within_target.csv'))

KNOWN = {'POLR2A', 'CLPB', 'EDNRB', 'NPR1', 'C3', 'MDM2', 'ITGAV', 'ITGB3', 'PPIA',
         'FKBP1A', 'PLXNB1'}
d['is_known_target'] = d.target_id.isin(KNOWN)
d['dock_usable'] = d.dock_score.notna() & (d.dock_score <= 0)

print(f'total rows: {len(d)}')
rows = []
for _, r in d.iterrows():
    reasons = []
    if r.is_known_target:
        reasons.append('known Lasso target (positive control, reported separately)')
    if r.peptide_bias == 'high':
        reasons.append(f'peptide bias high ({int(r.peptide_n_targets_ranked_1st)} targets)')
    if not r.dock_usable:
        reasons.append('no usable docking score')
    if r.rank_in_target > 3:
        reasons.append(f'within-target rank {int(r.rank_in_target)} > 3')
    rows.append('; '.join(reasons))
d['exclusion'] = rows
d['eligible'] = d.exclusion == ''

print(f'eligible after filters: {int(d.eligible.sum())} / {len(d)}')
print('\nwhy rows were excluded (counts):')
from collections import Counter
c = Counter()
for e in d[~d.eligible].exclusion:
    for part in e.split('; '):
        c[part.split(' (')[0]] += 1
for k, v in c.most_common():
    print(f'  {v:4d}  {k}')

el = d[d.eligible].copy()
if len(el):
    # signal agreement: both signals on the same side of their target's median
    el['agree'] = ((el.binding_z_target > 0) & (el.dock_z_target > 0)).astype(int)
    def z(s):
        sd = s.std()
        return (s - s.mean()) / sd if sd and sd > 1e-9 else s * 0
    el['sel_score'] = (0.40 * z(el.composite_target)
                       + 0.30 * z(-el.rank_in_target.astype(float))
                       + 0.30 * el.agree)
    el = el.sort_values('sel_score', ascending=False)

    # diversity: at most 2 per target
    picked, counts = [], {}
    for _, r in el.iterrows():
        if counts.get(r.target_id, 0) >= 2:
            continue
        counts[r.target_id] = counts.get(r.target_id, 0) + 1
        picked.append(r)
        if len(picked) >= 10:
            break
    top = pd.DataFrame(picked)
else:
    top = el

print(f'\n=== SELECTED TOP {len(top)} ===')
cols = ['target_id', 'peptide_id', 'rank_in_target', 'composite_target', 'binding_prob',
        'dock_score', 'peptide_bias', 'agree']
cols = [c for c in cols if c in top.columns]
pd.set_option('display.width', 220)
print(top[cols].round(3).to_string(index=False))

# reasons for inclusion, per row
def reason(r):
    bits = [f'rank {int(r.rank_in_target)}/{int(r.n_peptides_for_target)} for {r.target_id}',
            f'bias {r.peptide_bias}',
            f'docking {r.dock_score:.2f}']
    if r.get('agree', 0) == 1:
        bits.append('binding and docking agree')
    return '; '.join(bits)

if len(top):
    top = top.assign(selection_reason=top.apply(reason, axis=1))
    out = os.path.join(ROOT, 'results', 'final_candidates_top10.csv')
    top.drop(columns=[c for c in ['exclusion', 'eligible'] if c in top.columns]) \
       .to_csv(out, index=False, lineterminator='\n')
    print(f'\nwrote {out}')
    print('\nselection reasons:')
    for _, r in top.iterrows():
        print(f'  {r.target_id:8s} x {r.peptide_id:14s} {r.selection_reason}')

# also save the audit trail
d.to_csv(os.path.join(ROOT, 'results', 'candidates_screen50_audited.csv'),
         index=False, lineterminator='\n')
print('\nwrote results/candidates_screen50_audited.csv (every row with its exclusion reason)')
