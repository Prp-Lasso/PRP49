"""Precise target-ID audit.

Two distinct questions:
  1) which prot_id values are genuinely malformed (PDB chain descriptors) vs
     legitimate UniProt accessions, possibly with a domain range?
  2) do several different prot_id values denote the SAME biological target
     (e.g. 'EDNRB', 'P24530', '9KDF:R')? If so, group-by-prot_id CV leaks the
     target across folds and the 0.643 "new target" number would be optimistic.
"""
import os
import re

import pandas as pd

ROOT = r'D:\deepseek_harness\prp49'
df = pd.read_csv(os.path.join(ROOT, 'mvp_cpu', 'train_pairs_hard.csv'))

UNIPROT = re.compile(r'^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$')
CHAINISH = re.compile(r'^(chain|chains)\b', re.I)


def classify(v):
    v = str(v)
    base = v.split('[')[0].split(':')[0]
    if CHAINISH.match(v):
        return 'pdb_chain_descriptor'
    if UNIPROT.match(base):
        return 'uniprot' + ('+range' if '[' in v else '')
    if re.match(r'^[0-9][A-Z0-9]{3}$', base):
        return 'pdb_id' + ('+chain' if ':' in v else '')
    return 'other'


df['id_class'] = df.prot_id.map(classify)
print('=== prot_id classes (train_pairs_hard.csv) ===')
print(df.id_class.value_counts().to_string())
print()
for c in df.id_class.unique():
    ex = df[df.id_class == c].prot_id.unique()[:6]
    print(f'  {c:24s} n={len(df[df.id_class==c]):5d}  examples {list(ex)}')

print('\n=== do distinct prot_ids share the same prot_seq (i.e. same target)? ===')
seq2ids = df.groupby('prot_seq')['prot_id'].nunique().sort_values(ascending=False)
multi = seq2ids[seq2ids > 1]
print(f'  prot_seq values: {df.prot_seq.nunique()} | prot_id values: {df.prot_id.nunique()}')
print(f'  protein sequences carrying MORE THAN ONE prot_id: {len(multi)}')
for seq, n in multi.head(8).items():
    ids = sorted(df[df.prot_seq == seq].prot_id.unique())
    labs = df[df.prot_seq == seq].label.value_counts().to_dict()
    print(f'    {n} ids, labels {labs}: {ids[:6]}')

print('\n=== same check for the A+C affinity data (used with group-by-prot_id CV) ===')
aff = pd.read_csv(os.path.join(ROOT, 'mvp_cpu', 'affinity_pairs.csv'))
aff['id_class'] = aff.prot_id.map(classify)
print(aff.id_class.value_counts().to_string())
s2i = aff.groupby('prot_seq')['prot_id'].nunique()
print(f'  prot_seq {aff.prot_seq.nunique()} vs prot_id {aff.prot_id.nunique()} '
      f'| sequences with >1 id: {int((s2i>1).sum())}')
print('\n  impact:')
print('  train_pairs_hard is used with group-by-peptide CV -> target-ID mixing does NOT leak.')
print('  affinity_pairs is used with group-by-prot_id CV -> any target carrying 2+ ids')
print('  splits that target across folds, so the reported 0.643 is an upper bound.')
