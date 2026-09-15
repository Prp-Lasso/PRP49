"""Rebuild route-A propagation on top of the sequence-keyed pairs (v2).

Differences vs the first attempt:
  * targets keyed by SEQUENCE (238 keys) instead of ambiguous chain descriptors
  * identity recomputed with a proper global-ish alignment
  * soft labels graded by identity, and rows whose source target is a chain
    descriptor are simply impossible now
"""
import os

import numpy as np
import pandas as pd

ROOT = r'D:\deepseek_harness\prp49'
v2 = pd.read_csv(os.path.join(ROOT, 'mvp_cpu', 'train_pairs_v2.csv'))
print(f'v2 pairs {len(v2)} | targets {v2.target_key.nunique()} | families {v2.family_id.nunique()}')

members = v2.drop_duplicates('pep_seq').groupby('family_id')['pep_seq'].apply(list).to_dict()
known = set(zip(v2.pep_seq, v2.target_key))
tgt_seq = v2.drop_duplicates('target_key').set_index('target_key')['prot_seq'].to_dict()

pos = v2[v2.label == 1]
print(f'positives {len(pos)} | families with >1 member: {sum(1 for m in members.values() if len(m)>1)}')


def identity(a, b):
    """Ungapped best-overlap identity.

    NOTE: the comparison MUST carry `if x == y` - without it the generator counts
    zipped elements (always len(a)) and every sequence pair scores 1.0, which
    silently pushed every propagated row into the top soft-label tier.
    """
    if len(a) > len(b):
        a, b = b, a
    if not a:
        return 0.0
    best = 0.0
    for o in range(len(b) - len(a) + 1):
        win = b[o:o + len(a)]
        best = max(best, sum(1 for x, y in zip(a, win) if x == y) / len(a))
        if best == 1.0:
            break
    return best


rows = []
for _, r in pos.iterrows():
    fam_id, tkey, src = r['family_id'], r['target_key'], r['pep_seq']
    if len(members.get(fam_id, [])) < 2:
        continue
    for p2 in members[fam_id]:
        if p2 == src or (p2, tkey) in known:
            continue
        ident = identity(src, p2)
        if ident < 0.5:
            continue
        soft = 0.90 if ident >= 0.8 else (0.70 if ident >= 0.6 else 0.55)
        rows.append(dict(pep_id=f'prop-{len(rows)}', pep_seq=p2,
                         prot_seq=tgt_seq[tkey], target_key=tkey, label=soft,
                         family_id=fam_id, source_peptide=src,
                         identity=round(float(ident), 3), propagated=True,
                         neg_kind='propagated_positive'))

prop = pd.DataFrame(rows)
orig = v2.assign(propagated=False, source_peptide='', identity=1.0)
out = pd.concat([orig, prop], ignore_index=True).sample(frac=1.0, random_state=42)
out = out[['pep_id', 'pep_seq', 'prot_seq', 'target_key', 'label', 'family_id',
           'propagated', 'source_peptide', 'identity', 'neg_kind']]
dst = os.path.join(ROOT, 'mvp_cpu', 'train_pairs_prop2.csv')
out.to_csv(dst, index=False, lineterminator='\n')

n_orig_pos = int(v2.label.sum())
n_prop_pos = int((prop.label >= 0.5).sum()) if len(prop) else 0
print(f'\nwrote {dst}: {len(out)} pairs')
print(f'  original positives : {n_orig_pos}')
print(f'  propagated positives: {n_prop_pos}  (families {prop.family_id.nunique() if len(prop) else 0}, '
      f'targets {prop.target_key.nunique() if len(prop) else 0})')
if len(prop):
    print(f'  soft tiers: {prop.label.value_counts().to_dict()}')
    print(f'  identity: mean {prop.identity.mean():.3f} median {prop.identity.median():.3f} '
          f'min {prop.identity.min():.3f}')
    truly_new = set(zip(prop.pep_seq, prop.target_key)) - known
    print(f'  truly new (peptide,target) combos: {len(truly_new)}')
    print(prop.head(5)[['pep_seq', 'target_key', 'label', 'source_peptide', 'identity']].to_string(index=False))
print(f'\n  => total positives for training: {n_orig_pos + n_prop_pos}')
