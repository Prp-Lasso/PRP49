"""Build within-target pairwise ranking data for drug-candidate prioritisation.

Task: given two peptides assayed against the SAME target, which one binds stronger?
Rationale for a druggability-screening goal: absolute Kd prediction is not what a
screening funnel needs — it needs a correct *ranking* so the top-k candidates go
to the bench. Pairwise comparison is also the most natural use of affinity data
and cancels cross-target scale differences (Ki/Kd/IC50) by construction.

Sampling:
  - all within-target peptide pairs, |delta pAffinity| >= MIN_DELTA (drop near-ties)
  - cap per target to avoid a few 100+ pair targets dominating
Output: prot_id, prot_seq, pep_a_seq, pep_b_seq, energy_a, energy_b, label (1 = a stronger)
"""
import os
import numpy as np
import pandas as pd

ROOT = r'D:\deepseek_harness\prp49'
MIN_DELTA = 0.3          # pAffinity units (~2x affinity) — below this is assay noise
MAX_PER_TARGET = 200
GLOBAL_CAP = 60000

df = pd.read_csv(os.path.join(ROOT, 'mvp_cpu', 'affinity_pairs.csv'))
df = df.dropna(subset=['energy', 'prot_seq', 'pep_seq'])
print(f'input pairs: {len(df)} | targets: {df.prot_id.nunique()}')

rng = np.random.default_rng(42)
rows = []
for pid, grp in df.groupby('prot_id'):
    if len(grp) < 2:
        continue
    recs = grp[['pep_seq', 'energy']].drop_duplicates('pep_seq').to_dict('records')
    n = len(recs)
    pairs = []
    for i in range(n):
        for j in range(i + 1, n):
            ea, eb = recs[i]['energy'], recs[j]['energy']
            if abs(ea - eb) < MIN_DELTA:
                continue
            if ea > eb:
                pairs.append((recs[i]['pep_seq'], recs[j]['pep_seq'], ea, eb, 1))
            else:
                pairs.append((recs[j]['pep_seq'], recs[i]['pep_seq'], eb, ea, 1))
    if not pairs:
        continue
    if len(pairs) > MAX_PER_TARGET:
        idx = rng.choice(len(pairs), size=MAX_PER_TARGET, replace=False)
        pairs = [pairs[k] for k in idx]
    prot_seq = grp.prot_seq.iloc[0]
    for a, b, ea, eb, lab in pairs:
        rows.append(dict(prot_id=pid, prot_seq=prot_seq, pep_a_seq=a, pep_b_seq=b,
                         energy_a=ea, energy_b=eb, label=lab))

out = pd.DataFrame(rows).sample(frac=1.0, random_state=42).reset_index(drop=True)
if len(out) > GLOBAL_CAP:
    out = out.iloc[:GLOBAL_CAP]
print(f'pairwise samples: {len(out)} | targets: {out.prot_id.nunique()}')
print(f'  pairs per target: min {out.groupby("prot_id").size().min()}, '
      f'median {int(out.groupby("prot_id").size().median())}, max {out.groupby("prot_id").size().max()}')
print(f'  delta pAffinity: mean {(out.energy_a-out.energy_b).mean():.2f}, '
      f'min {(out.energy_a-out.energy_b).min():.2f}, max {(out.energy_a-out.energy_b).max():.2f}')

dst = os.path.join(ROOT, 'mvp_cpu', 'affinity_pairs_rank.csv')
out.to_csv(dst, index=False, lineterminator='\n')
print(f'wrote {dst}')
