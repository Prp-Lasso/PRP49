"""Rebuild the training pairs with an unambiguous target key.

Problem found by the audit: 81% of prot_id values in train_pairs_hard.csv are PDB
chain descriptors ("Chains A, B"), and 238 distinct protein sequences collapse
onto only 49 ids. The hard negatives were built as "same target, different
peptide family", so an ambiguous target key means some negatives are not really
same-target pairs.

Fix: key every target by its SEQUENCE (the biologically correct identity).
  * positives keep their peptide/target combination
  * negatives are rebuilt as "same target sequence, different peptide family"
  * an explicit target_id -> sequence map is written for traceability
"""
import os

import numpy as np
import pandas as pd

ROOT = r'D:\deepseek_harness\prp49'
src = pd.read_csv(os.path.join(ROOT, 'mvp_cpu', 'train_pairs_hard.csv'))
fam = pd.read_csv(os.path.join(ROOT, 'mvp_cpu', 'train_pairs_family.csv'))
print(f'input: {len(src)} pairs | prot_id {src.prot_id.nunique()} | prot_seq {src.prot_seq.nunique()}')

fam_map = fam.drop_duplicates('pep_seq').set_index('pep_seq')['family_id'].to_dict()
src['family_id'] = src.pep_seq.map(fam_map).fillna('unknown')

# stable target id from the sequence itself
tgt_ids = {s: f'T{i:04d}' for i, s in enumerate(sorted(src.prot_seq.unique()))}
src['target_key'] = src.prot_seq.map(tgt_ids)
print(f'target keys (by sequence): {src.target_key.nunique()}')

pos = src[src.label == 1].copy()
print(f'positives: {len(pos)} | distinct (pep, target_key): {pos.groupby(["pep_seq","target_key"]).ngroups}')

# ---- rebuild negatives: same target sequence, DIFFERENT peptide family ----
keep_pos = pos[['pep_id', 'pep_seq', 'prot_seq', 'target_key', 'label', 'family_id']].copy()
keep_pos['neg_kind'] = 'positive'

rng = np.random.default_rng(42)
neg_rows = []
for tkey, grp in src.groupby('target_key'):
    fams = grp.drop_duplicates('family_id')[['family_id', 'pep_seq']]
    peps_by_fam = grp.groupby('family_id')['pep_seq'].apply(list).to_dict()
    fam_list = list(peps_by_fam.keys())
    if len(fam_list) < 2:
        continue
    for fam_a, peps_a in peps_by_fam.items():
        for p in peps_a:
            other_fams = [f for f in fam_list if f != fam_a]
            if not other_fams:
                continue
            f_b = other_fams[rng.integers(len(other_fams))]
            q = peps_by_fam[f_b][rng.integers(len(peps_by_fam[f_b]))]
            neg_rows.append(dict(pep_id=f'hn-{len(neg_rows)}', pep_seq=q,
                                 prot_seq=grp.prot_seq.iloc[0], target_key=tkey,
                                 label=0, family_id=f_b, neg_kind='same_target_other_family'))
neg = pd.DataFrame(neg_rows)
print(f'rebuilt negatives: {len(neg)} (was {int((src.label==0).sum())})')

out = pd.concat([keep_pos, neg], ignore_index=True).sample(frac=1.0, random_state=42)
out = out[['pep_id', 'pep_seq', 'prot_seq', 'target_key', 'label', 'family_id', 'neg_kind']]
dst = os.path.join(ROOT, 'mvp_cpu', 'train_pairs_v2.csv')
out.to_csv(dst, index=False, lineterminator='\n')

# map for traceability
pd.DataFrame([{'target_key': k, 'prot_seq': s} for s, k in tgt_ids.items()]).to_csv(
    os.path.join(ROOT, 'mvp_cpu', 'target_key_map.csv'), index=False, lineterminator='\n')

print(f'\nwrote {dst}: {len(out)} pairs')
print(f'  positives {int(out.label.sum())} | negatives {int((out.label==0).sum())}')
print(f'  targets {out.target_key.nunique()} | families {out.family_id.nunique()}')
print(f'  neg kinds: {out.neg_kind.value_counts().to_dict()}')
ok = out.groupby('target_key').family_id.nunique()
print(f'  targets with >=2 families (valid negatives): {int((ok>=2).sum())} / {len(ok)}')
