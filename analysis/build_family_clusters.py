"""Route B (corrected): family grouping via sequence clustering.

The LassoPred `Lasso_Peptide_Family` column is unusable - 4,442 of 4,749 rows are
NaN and only 49 families are named. So families are defined here by clustering the
training peptides themselves at an identity threshold, which is what actually
matters for CV leakage.

Protocol ladder this completes:
  random KFold          -> loosest (not run)
  group by peptide      -> 0.8333  "new peptide x known target"   (current headline)
  group by family       -> ???     "new family  x known target"   (this run)
  group by target       -> 0.643   "known peptide x new target"   (A+C affinity task)
"""
import os
from itertools import combinations

import numpy as np
import pandas as pd

ROOT = r'D:\deepseek_harness\prp49'
pairs = pd.read_csv(os.path.join(ROOT, 'mvp_cpu', 'train_pairs_hard.csv'))
peps = sorted(pairs.pep_seq.unique(), key=len, reverse=True)
print(f'peptides: {len(peps)} | length {min(map(len,peps))}-{max(map(len,peps))}')


def identity(a, b):
    """Best identity over the overlap (ungapped sliding window)."""
    if len(a) > len(b):
        a, b = b, a
    if not a:
        return 0.0
    best = 0.0
    for off in range(len(b) - len(a) + 1):
        win = b[off:off + len(a)]
        best = max(best, sum(1 for x, y in zip(a, win) if x == y) / len(a))
        if best == 1.0:
            break
    return best


for THRESH in (0.7, 0.5, 0.3):
    clusters = []
    for s in peps:
        placed = False
        for c in clusters:
            if identity(s, c[0]) >= THRESH:
                c.append(s)
                placed = True
                break
        if not placed:
            clusters.append([s])
    sizes = sorted((len(c) for c in clusters), reverse=True)
    multi = [c for c in clusters if len(c) > 1]
    covered = sum(len(c) for c in multi)
    print(f'\nidentity >= {THRESH}: {len(clusters)} families '
          f'({len(multi)} with >1 member, covering {covered} peptides)')
    print(f'  size distribution: {sizes[:12]}')
    for c in sorted(multi, key=len, reverse=True)[:4]:
        print(f'    family of {len(c)}: {[s[:14] for s in c[:5]]}')
    if THRESH == 0.5:
        best = clusters

fam_of = {}
for i, c in enumerate(best):
    for s in c:
        fam_of[s] = f'fam{i:03d}' if len(c) > 1 else f'solo_{s[:8]}'

out = pairs.copy()
out['family_id'] = out.pep_seq.map(fam_of)
fp = os.path.join(ROOT, 'mvp_cpu', 'train_pairs_family.csv')
out.to_csv(fp, index=False, lineterminator='\n')

pos = out[out.label == 1]
sizes = pos.groupby('family_id').size().sort_values(ascending=False)
print(f'\nwrote {fp}: {len(out)} pairs | {out.family_id.nunique()} families (threshold 0.5)')
print(f'positives: {len(pos)} in {pos.family_id.nunique()} families')
print('largest positive families:')
print(sizes.head(8).to_string())
print(f'\npositives in multi-member families: {int(pos.family_id.str.startswith("fam").sum())} / {len(pos)}')
print('  -> these are the pairs that COULD have leaked across folds in the old peptide-grouped CV')
