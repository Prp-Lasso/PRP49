"""Build a specificity-enhanced training set.

Diagnosis: the model assigns 0.77-0.88 binding probability to wild-type MccJ25
against integrins, where it is experimentally inactive (>10 uM). The current
negatives are "same target, different peptide" (hard negatives for target
discrimination) but nothing teaches "THIS peptide binds only ITS OWN target".

Fix: add cross-target negatives - each positive pair (p, t) is paired with the
same peptide against unrelated targets t' != t. Combined with the audit's
confirmed experimental negatives, this directly supervises specificity.

Output: mvp_cpu/train_pairs_xneg.csv
"""
import os
import numpy as np
import pandas as pd

ROOT = r'D:\deepseek_harness\prp49'
src = pd.read_csv(os.path.join(ROOT, 'mvp_cpu', 'train_pairs_hard.csv'))
print(f'input: {len(src)} pairs, labels {src.label.value_counts().to_dict()}')

pos = src[src.label == 1].copy()
neg = src[src.label == 0].copy()
print(f'positives {len(pos)} | existing negatives {len(neg)}')

rng = np.random.default_rng(42)

# ---- cross-target negatives: same peptide, a different (unrelated) target
targets = pos[['prot_id', 'prot_seq']].drop_duplicates('prot_id').to_dict('records')
print(f'distinct targets among positives: {len(targets)}')

rows = []
for _, p in pos.iterrows():
    cand = [t for t in targets if t['prot_id'] != p['prot_id']]
    if not cand:
        continue
    # one cross-target negative per positive (keeps the set balanced)
    t = cand[rng.integers(len(cand))]
    rows.append(dict(pep_id=f'xneg-{len(rows)}', pep_seq=p['pep_seq'], prot_seq=t['prot_seq'],
                     prot_id=t['prot_id'], label=0, neg_kind='cross_target'))
xneg = pd.DataFrame(rows)
print(f'cross-target negatives generated: {len(xneg)}')

# ---- audit-confirmed experimental negatives (wild-type MccJ25 on integrins)
MCCJ25 = 'GGAGHVPEYFVGIGTPISFYG'
exp_neg = pd.DataFrame([
    dict(pep_id='expneg-0', pep_seq=MCCJ25,
         prot_seq=open(os.path.join(ROOT, 'mvp_cpu', 'scan_targets.fasta')).read(),  # placeholder
         prot_id='ITGAV', label=0, neg_kind='experimental'),
])
exp_neg = exp_neg.iloc[0:0]      # built properly below if target sequences are available
fasta = os.path.join(ROOT, 'mvp_cpu', 'scan_targets.fasta')
if os.path.exists(fasta):
    seqs, name, buf = {}, None, []
    for line in open(fasta):
        line = line.strip()
        if line.startswith('>'):
            if name:
                seqs[name] = ''.join(buf)
            name, buf = line[1:].split()[0].split('|')[0], []
        elif line:
            buf.append(line)
    if name:
        seqs[name] = ''.join(buf)
    recs = []
    for tgt in ('ITGAV', 'ITGB3'):
        if tgt in seqs:
            recs.append(dict(pep_id=f'expneg-{tgt}', pep_seq=MCCJ25, prot_seq=seqs[tgt],
                             prot_id=tgt, label=0, neg_kind='experimental'))
    exp_neg = pd.DataFrame(recs)
print(f'audit-confirmed experimental negatives: {len(exp_neg)}')

out = pd.concat([pos.assign(neg_kind='positive'), neg.assign(neg_kind='hard_negative'),
                 xneg, exp_neg], ignore_index=True)
out = out.sample(frac=1.0, random_state=42).reset_index(drop=True)
for c in ('pep_id', 'pep_seq', 'prot_seq', 'prot_id', 'label', 'neg_kind'):
    if c not in out.columns:
        out[c] = ''
out = out[['pep_id', 'pep_seq', 'prot_seq', 'prot_id', 'label', 'neg_kind']]
dst = os.path.join(ROOT, 'mvp_cpu', 'train_pairs_xneg.csv')
out.to_csv(dst, index=False, lineterminator='\n')
print(f'\nwrote {dst}: {len(out)} pairs')
print(out.neg_kind.value_counts().to_dict())
print(out.label.value_counts().to_dict())
