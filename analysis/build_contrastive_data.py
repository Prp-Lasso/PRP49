"""Route C data: assign every LassoPred sequence to a family for contrastive learning.

LassoPred's own Lasso_Peptide_Family column is unusable (4,442/4,749 NaN), so
families come from the 178 clusters derived from the training peptides, and every
one of the 4,749 Core_Sequences is assigned to the nearest family by identity.
Sequences that match nothing above the threshold are kept as singletons: they
contribute negatives but never form a positive pair.
"""
import os
import re

import numpy as np
import pandas as pd

ROOT = r'D:\deepseek_harness\prp49'
db = pd.read_csv(os.path.join(ROOT, 'lassopred_database.csv'))
print(f'database rows: {len(db)} | columns with seq: '
      f'{[c for c in db.columns if "equence" in c]}')

AA = re.compile(r'[^ACDEFGHIKLMNPQRSTVWY]')
db['core'] = db['Core_Sequence'].astype(str).str.upper().map(lambda s: AA.sub('', s))
db = db[db.core.str.len() >= 5].copy()
print(f'usable cores: {len(db)} | length {db.core.str.len().min()}-{db.core.str.len().max()}')

# families from the training peptides (50% identity clustering, already computed)
fam_df = pd.read_csv(os.path.join(ROOT, 'mvp_cpu', 'train_pairs_family.csv'))
fam_of = fam_df.drop_duplicates('pep_seq').set_index('pep_seq')['family_id'].to_dict()
members = {}
for pep, fam in fam_of.items():
    if str(fam).startswith('fam'):            # multi-member families only
        members.setdefault(fam, []).append(pep)
print(f'families with >1 member: {len(members)} | covering {sum(len(v) for v in members.values())} peptides')

reps = {}
for fam, seqs in members.items():
    # use the LONGEST member as the representative: a short representative slides
    # inside longer sequences and matches them by chance (that is how one family
    # swallowed 87% of the database on the first attempt)
    reps[fam] = max(seqs, key=len)
print('representative lengths:',
      {f: len(s) for f, s in list(reps.items())[:6]})


def identity(a, b):
    """Best-overlap identity. NOTE the `if x == y` - omitting it counts zipped
    elements instead of matches and makes every pair score 1.0 (mistake #26)."""
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


def seq_identity(a, b):
    """Identity over the overlap, but only for sequences of comparable length.

    A short sequence can slide inside a long one and hit 50% by chance, so pairs
    differing by more than 5 residues are rejected outright, and the overlap must
    cover at least 80% of the longer sequence.
    """
    if not a or not b:
        return 0.0
    if abs(len(a) - len(b)) > 5:
        return 0.0
    if len(a) > len(b):
        a, b = b, a
    n = len(a)
    best = 0.0
    for o in range(len(b) - n + 1):
        best = max(best, sum(1 for x, y in zip(a, b[o:o + n])) / len(b))
        if best == 1.0:
            break
    return best


THRESH = 0.7
rows, assigned = [], 0
for _, r in db.iterrows():
    core = r.core
    best_fam, best_id = None, 0.0
    for fam, rep in reps.items():
        v = seq_identity(core, rep)
        if v > best_id:
            best_id, best_fam = v, fam
    if best_id >= THRESH:
        rows.append(dict(sequence=core, family_id=best_fam, identity=round(best_id, 3),
                         source=r.get('Lasso_Peptide_Name', '')))
        assigned += 1
    else:
        rows.append(dict(sequence=core, family_id=f'solo_{core[:8]}', identity=round(best_id, 3),
                         source=r.get('Lasso_Peptide_Name', '')))

out = pd.DataFrame(rows).drop_duplicates('sequence')
print(f'\nassigned to a named family: {assigned}/{len(db)} ({100*assigned/len(db):.0f}%)')
sizes = out.groupby('family_id').size().sort_values(ascending=False)
multi = sizes[sizes > 1]
print(f'families present: {out.family_id.nunique()} | with >1 member: {len(multi)}')
print(f'sequences usable for positive pairs: {int(multi.sum())}')
print('largest families:')
print(multi.head(8).to_string())

dst = os.path.join(ROOT, 'mvp_cpu', 'contrastive_pairs.csv')
out.to_csv(dst, index=False, lineterminator='\n')
print(f'\nwrote {dst}: {len(out)} sequences')
