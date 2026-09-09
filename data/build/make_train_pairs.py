"""Merge expanded_pairs with augmented positives -> train_pairs.csv.

Augmented variants inherit the parent row's prot_seq/pdb/label so that
grouped CV (by pdb) keeps variants in the same fold as the parent.
Usage: python scratch/make_train_pairs.py
"""
import os
import pandas as pd

root = r'D:\deepseek_harness\prp49'
exp = pd.read_csv(os.path.join(root, 'mvp_cpu', 'expanded_pairs.csv'))
aug = pd.read_csv(os.path.join(root, 'mvp_cpu', 'augmented_positives.csv'))

rows = [exp]
n_added = 0
for _, a in aug.iterrows():
    parents = exp[exp['pep_id'] == a['orig']]
    for _, p in parents.iterrows():
        rows.append(pd.DataFrame([{
            'pep_id': a['pep_id'], 'pdb': p['pdb'],
            'pep_seq': a['seq'], 'prot_seq': p['prot_seq'],
            'prot_id': p['prot_id'], 'label': 1}]))
        n_added += 1

out = pd.concat(rows, ignore_index=True)
out.to_csv(os.path.join(root, 'mvp_cpu', 'train_pairs.csv'), index=False)
print(f'train_pairs.csv: {len(out)} rows (expanded {len(exp)} + variants {n_added}); '
      f'positives {int(out["label"].sum())}')
