"""Route D v2: merge the real complexes from the new database into the auxiliary set.

D' (running) used only Propedia domain homologs - 199 receptors, no POLR2A.
D v2 keeps that signal but adds the six genuine lasso/protein complexes found in the
PDB (MccJ25+capistruin x RNAP beta-prime, lassomycin x ClpC1), oversampled so that a
handful of high-confidence pairs is not drowned by 4,000 Propedia pairs.

Design kept identical to D' (same config, same lambda_aux) so the only change is the
auxiliary data, which is what makes the comparison meaningful.
"""
import os

import pandas as pd

ROOT = r'D:\deepseek_harness\prp49'

homolog = pd.read_csv(os.path.join(ROOT, 'mvp_cpu', 'homolog_aux_resampled.csv'))
db_aux = pd.read_csv(os.path.join(ROOT, 'mvp_cpu', 'db_aux_pairs.csv'))
print(f"domain-homolog aux pairs : {len(homolog)}")
print(f"real-complex aux pairs   : {len(db_aux)} "
      f"({int((db_aux.label == 1).sum())} pos / {int((db_aux.label == 0).sum())} neg)")

# Oversample the real complexes: 26 pairs against 4,000 would contribute ~0.6% of the
# auxiliary gradient. x20 gives them ~12% weight while keeping the broader Propedia
# signal that D' showed to be useful for cross-domain generalisation.
OVERSAMPLE = 20
boosted = pd.concat([db_aux] * OVERSAMPLE, ignore_index=True)
print(f"oversampled real complexes x{OVERSAMPLE}: {len(boosted)}")

# align schemas on the columns train_mtl.py actually reads
KEEP = ['pep_seq', 'prot_seq', 'label', 'prot_id']
for df in (homolog, boosted):
    for c in KEEP:
        assert c in df.columns, f'missing {c}'

merged = pd.concat([homolog[KEEP], boosted[KEEP]], ignore_index=True)
# NO drop_duplicates: both sets encode weight through repetition (the domain-homolog
# set was resampled by homology weight, the real complexes are oversampled here).
# Deduplicating silently cancelled the oversampling - 520 boosted rows collapsed to 26.
merged = merged[merged.pep_seq.str.len().between(8, 80)]
out = os.path.join(ROOT, 'mvp_cpu', 'homolog_plus_db_aux.csv')
merged.to_csv(out, index=False, lineterminator='\n')
print(f'\nwrote {out}: {len(merged)} rows')
print(f'  label balance: {merged.label.value_counts().to_dict()}')
print(f'  unique peptides {merged.pep_seq.nunique()} | receptors {merged.prot_seq.nunique()}')

# how much of the set now comes from the real complexes?
real_keys = set(zip(db_aux.pep_seq, db_aux.prot_seq))
n_real = sum(1 for p, t in zip(merged.pep_seq, merged.prot_seq) if (p, t) in real_keys)
print(f'  rows traceable to a real PDB complex: {n_real} ({100*n_real/len(merged):.1f}%)')

# confirm the POLR2A-related pairs survive
rnap = merged[merged.prot_seq.str.len() > 1000]
print(f'  rows using a >1000 aa target (RNAP beta-prime): {len(rnap)}')
print(f'  of which labelled positive: {int(rnap.label.sum())}')
