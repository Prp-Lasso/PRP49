"""Build grouped CV splits for the teammate's dataset.

Why this instead of downloading their 2.64 GB checkpoint: the mirror delivers 0.04 MB/s
(ETA 16.7 h). They already have the weights and a working environment; what they lack is a
leakage-controlled split. Producing that split is minutes of work here and lets them obtain
the comparison number themselves.

Grouping key: the file names are "<protein>vs<peptide>.pdbqt". Their 243 pairs cover only
32 proteins (~7.6 pairs per protein), so a random split puts the same protein in both train
and test. Grouping by protein removes that leak - the same fix we applied to our own
pipeline (stratified 0.8145 -> grouped 0.8333).
"""
import os
import re
from collections import Counter

import numpy as np
import pandas as pd

ROOT = r'D:\deepseek_harness\prp49'
TEAM = os.path.join(ROOT, 'external', 'LioXue_PRP49', 'data')

frames = []
for f, split in [('out_train.csv', 'their_train'), ('out_test.csv', 'their_test')]:
    p = os.path.join(TEAM, f)
    d = pd.read_csv(p)
    d['source_file'] = split
    frames.append(d)
d = pd.concat(frames, ignore_index=True)
print(f'combined: {len(d)} rows (their train {sum(d.source_file=="their_train")} '
      f'+ test {sum(d.source_file=="their_test")})')

# protein / peptide ids from "1BI7vs2mw3.pdbqt"
pat = re.compile(r'^([0-9A-Za-z]{4})vs(.+)\.pdbqt$', re.I)


def parse(fn):
    m = pat.match(str(fn))
    if m:
        return m.group(1).upper(), m.group(2).upper()
    parts = str(fn).split('vs')
    return (parts[0].upper() if parts else str(fn)), (parts[1].replace('.pdbqt', '').upper() if len(parts) > 1 else '')


parsed = d.filename.map(parse)
d['protein_id'] = [p[0] for p in parsed]
d['peptide_id'] = [p[1] for p in parsed]

print(f'unique proteins: {d.protein_id.nunique()} | unique peptides: {d.peptide_id.nunique()}')
per_prot = d.groupby('protein_id').size()
print(f'pairs per protein: min {per_prot.min()} median {per_prot.median():.0f} max {per_prot.max()}')
per_pep = d.groupby('peptide_id').size()
print(f'pairs per peptide: min {per_pep.min()} median {per_pep.median():.0f} max {per_pep.max()}')

# ---- grouped 5-fold by PROTEIN (the main leak source) ----
rng = np.random.RandomState(42)
prots = np.array(sorted(d.protein_id.unique()))
rng.shuffle(prots)
fold_of_prot = {}
for i, p in enumerate(prots):
    fold_of_prot[p] = i % 5
d['fold_by_protein'] = d.protein_id.map(fold_of_prot)

# ---- grouped 5-fold by PEPTIDE (secondary check) ----
peps = np.array(sorted(d.peptide_id.unique()))
rng.shuffle(peps)
fold_of_pep = {p: i % 5 for i, p in enumerate(peps)}
d['fold_by_peptide'] = d.peptide_id.map(fold_of_pep)

print('\n=== fold sizes (by protein) ===')
for f, g in d.groupby('fold_by_protein'):
    print(f'  fold {f}: {len(g):4d} pairs | {g.protein_id.nunique():3d} proteins | '
          f'{g.peptide_id.nunique():3d} peptides')

# ---- verify no protein spans folds ----
leak = d.groupby('protein_id').fold_by_protein.nunique()
print(f'\nproteins appearing in more than one fold: {int((leak > 1).sum())} (must be 0)')
assert (leak > 1).sum() == 0, 'protein leakage in the protein-grouped split'
leak2 = d.groupby('peptide_id').fold_by_peptide.nunique()
print(f'peptides appearing in more than one fold (peptide split): {int((leak2 > 1).sum())} (must be 0)')
assert (leak2 > 1).sum() == 0

# ---- what the OLD random split would have leaked (for the paper) ----
from sklearn.model_selection import train_test_split
idx_tr, idx_te = train_test_split(np.arange(len(d)), test_size=0.2, random_state=42)
leaked_prot = set(d.protein_id.iloc[idx_tr]) & set(d.protein_id.iloc[idx_te])
leaked_pep = set(d.peptide_id.iloc[idx_tr]) & set(d.peptide_id.iloc[idx_te])
print(f'\n=== what a random 80/20 split shares between train and test ===')
print(f'  proteins in both: {len(leaked_prot)} of {d.protein_id.nunique()} '
      f'({100*len(leaked_prot)/d.protein_id.nunique():.0f}%)')
print(f'  peptides in both: {len(leaked_pep)} of {d.peptide_id.nunique()} '
      f'({100*len(leaked_pep)/d.peptide_id.nunique():.0f}%)')

out = os.path.join(ROOT, 'external', 'teammate_grouped_cv_splits.csv')
d[['filename', 'protein seq', 'lasso seq', 'energy', 'protein_id', 'peptide_id',
   'fold_by_protein', 'fold_by_peptide', 'source_file']].to_csv(out, index=False, lineterminator='\n')
print(f'\nwrote {out}')
print('columns: fold_by_protein / fold_by_peptide are 0-4 for 5-fold grouped CV')
