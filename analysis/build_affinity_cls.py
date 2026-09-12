"""Build A+C dataset: target-wise z-scored affinity, binned into a clean binary task.

A (classification): drop the ambiguous middle band, keep strong vs weak binders.
C (normalisation): z-score pAffinity *within each target* so Ki/Kd/IC50 scales and
   cross-target affinity differences no longer inflate the loss.
"""
import os
import numpy as np
import pandas as pd

ROOT = r'D:\deepseek_harness\prp49'
src = os.path.join(ROOT, 'mvp_cpu', 'affinity_pairs.csv')
df = pd.read_csv(src)
print(f'input: {len(df)} pairs, {df.prot_id.nunique()} targets')

# --- C: keep targets with enough pairs for a meaningful within-target z-score
sizes = df.groupby('prot_id').size()
keep = sizes[sizes >= 5].index
d = df[df.prot_id.isin(keep)].copy()
print(f'after target filter (>=5 pairs): {len(d)} pairs, {d.prot_id.nunique()} targets')

g = d.groupby('prot_id')['energy']
d['z'] = (d['energy'] - g.transform('mean')) / g.transform('std').replace(0, np.nan)
d = d.dropna(subset=['z'])

# --- A: bin into strong / weak, drop the ambiguous band
strong = d[d.z >= 0.5].copy()
weak = d[d.z <= -0.5].copy()
strong['label'] = 1
weak['label'] = 0
out = pd.concat([strong, weak]).sample(frac=1.0, random_state=42).reset_index(drop=True)
out['pep_id'] = [f'affc-{i}' for i in range(len(out))]

cols = ['pep_id', 'pep_seq', 'prot_seq', 'prot_id', 'energy', 'z', 'label',
        'affinity_type', 'source']
out = out[cols]
dst = os.path.join(ROOT, 'mvp_cpu', 'affinity_cls_pairs.csv')
out.to_csv(dst, index=False, lineterminator='\n')

pos, neg = int((out.label == 1).sum()), int((out.label == 0).sum())
print(f'output: {len(out)} pairs -> {dst}')
print(f'  positive (z>=+0.5): {pos}   negative (z<=-0.5): {neg}   ratio {pos/max(neg,1):.2f}')
print(f'  targets: {out.prot_id.nunique()}   peptides: {out.pep_seq.nunique()}')
print(f'  affinity types: {out.affinity_type.value_counts().to_dict()}')
print(f'  z range: {out.z.min():.2f} .. {out.z.max():.2f}')
tsz = out.groupby('prot_id').size()
print(f'  targets with both classes: {out.groupby("prot_id").label.nunique().eq(2).sum()}')
print(f'  pairs per target: min {tsz.min()}, median {int(tsz.median())}, max {tsz.max()}')
