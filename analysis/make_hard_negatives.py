"""Build hard negatives (I5): same target x peptides known to bind OTHER targets.

Rationale: current negatives are random re-pairings + scrambled peptides, which
let the model exploit composition shortcuts. Pairs that keep the target fixed
while swapping in a real lasso peptide force the model to use peptide-target
compatibility rather than peptide composition alone.
"""
import os
import pandas as pd

ROOT = r'D:\deepseek_harness\prp49'
exp = pd.read_csv(os.path.join(ROOT, 'mvp_cpu', 'expanded_pairs.csv'))
pos = exp[exp['label'] == 1].copy()
print(f'positives: {len(pos)}, targets: {pos["prot_id"].nunique()}')

# target -> set of peptide sequences known to bind it
known = pos.groupby('prot_id')['pep_seq'].apply(set).to_dict()
# peptide -> its target(s)
pep_targets = pos.groupby('pep_seq')['prot_id'].apply(set).to_dict()

rows = []
MAX_PER_TARGET = 40
for tgt, tgt_peps in known.items():
    prot_seq = pos[pos['prot_id'] == tgt].iloc[0]['prot_seq']
    cand = [p for p, targets in pep_targets.items()
            if tgt not in targets and p not in tgt_peps][:MAX_PER_TARGET]
    for p in cand:
        rows.append(dict(pep_id=f'hard-{tgt}-{len(rows)}', pdb='hard-neg',
                         pep_seq=p, prot_seq=prot_seq, prot_id=tgt, label=0))
hard = pd.DataFrame(rows)
hard.to_csv(os.path.join(ROOT, 'mvp_cpu', 'hard_negatives.csv'), index=False)
print(f'hard negatives: {len(hard)} pairs across {hard["prot_id"].nunique()} targets')

# merge into a training set: train_pairs (expanded + augmented) + hard negatives
train = pd.read_csv(os.path.join(ROOT, 'mvp_cpu', 'train_pairs.csv'))
merged = pd.concat([train, hard], ignore_index=True)
merged.to_csv(os.path.join(ROOT, 'mvp_cpu', 'train_pairs_hard.csv'), index=False)
print(f'train_pairs_hard.csv: {len(merged)} rows '
      f'(pos {int(merged["label"].sum())}, neg {int((merged["label"] == 0).sum())})')
