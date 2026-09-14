"""Build a Propedia multi-task TRAINING set, disjoint from the cross-domain TEST set.

Task 3 motivation: the pairwise ranking head trained on ordinary peptides actively
hurts on lasso peptides (known-pair mean rank 23 -> 71) while the lasso classifier
is useless on ordinary peptides (AUC 0.47). Both point at domain-specific
representations. Multi-task training asks whether one shared BAN trunk can learn
something transferable when supervised on both domains at once.

Design (keeps the input format identical for both heads so one forward pass serves
both): each domain is a peptide-target binary classification problem.
  head A (lasso)    : train_pairs_hard.csv        - the existing in-domain set
  head B (ordinary) : Propedia complexes vs shuffled peptides, receptors split by
                      group so no receptor appears in both train and test

Split rule: receptors are partitioned 80/20; the 20% held-out receptors become the
new cross-domain TEST set (the previous crossdomain_test.csv used all receptors, so
it cannot be reused for training).
"""
import os
import numpy as np
import pandas as pd

ROOT = r'D:\deepseek_harness\prp49'
src = pd.read_csv(os.path.join(ROOT, 'docking_propedia', 'all_scores.csv'))
src = src[(src.n_pep >= 5) & (src.n_pep <= 40)]
src = src[~src.pep_seq.str.contains('X') & ~src.rec_seq.str.contains('X')]
src = src[src.vina_score <= -12.0]
src = src.drop_duplicates(['rec_seq', 'pep_seq'])
print(f'Propedia positives (experimental complexes): {len(src)}')

rng = np.random.default_rng(7)
recs = src.rec_seq.unique()
rng.shuffle(recs)
n_test = max(1, int(round(0.2 * len(recs))))
test_recs, train_recs = set(recs[:n_test]), set(recs[n_test:])
train_pos = src[src.rec_seq.isin(train_recs)].copy()
test_pos = src[src.rec_seq.isin(test_recs)].copy()
print(f'receptors: {len(recs)} -> train {len(train_recs)} / test {len(test_recs)}')
print(f'positives: train {len(train_pos)} / test {len(test_pos)}')


def shuffled(seq):
    a = list(seq)
    rng.shuffle(a)
    return ''.join(a)


train_neg = train_pos[['rec_seq', 'pep_seq']].copy()
train_neg['pep_seq'] = train_neg.pep_seq.map(shuffled)
test_neg = test_pos[['rec_seq', 'pep_seq']].copy()
test_neg['pep_seq'] = test_neg.pep_seq.map(shuffled)

# a second, harder negative family for the test set: another real peptide on the
# same receptor (only for receptors with >= 2 complexes)
hard = []
by_rec = test_pos.groupby('rec_seq')['pep_seq'].apply(list).to_dict()
for rec, peps in by_rec.items():
    if len(peps) < 2:
        continue
    for i, p in enumerate(peps):
        hard.append(dict(rec_seq=rec, pep_seq=rng.choice([q for j, q in enumerate(peps) if j != i])))
hard = pd.DataFrame(hard)


def frame(pos, neg, labels):
    a = pos[['rec_seq', 'pep_seq']].copy(); a['label'] = 1
    b = neg[['rec_seq', 'pep_seq']].copy(); b['label'] = 0
    out = pd.concat([a, b], ignore_index=True).rename(columns={'rec_seq': 'prot_seq'})
    # full-sequence hash as the group id: a 12-char prefix collides between
    # unrelated receptors and produced 91 phantom train/test overlaps
    out['prot_id'] = out.prot_seq.map(lambda s: 'r' + str(abs(hash(s)) % 10**10))
    return out.sample(frac=1.0, random_state=42).reset_index(drop=True)


train = frame(train_pos, train_neg, None)
test = frame(test_pos, test_neg, None)
train.to_csv(os.path.join(ROOT, 'mvp_cpu', 'propedia_train_mtl.csv'), index=False, lineterminator='\n')
test.to_csv(os.path.join(ROOT, 'mvp_cpu', 'propedia_test_mtl.csv'), index=False, lineterminator='\n')
if len(hard):
    h = hard.rename(columns={'rec_seq': 'prot_seq'}).copy()
    h['label'] = 0
    h['prot_id'] = h.prot_seq.map(lambda x: 'r' + str(abs(hash(x)) % 10**10))
    h.to_csv(os.path.join(ROOT, 'mvp_cpu', 'propedia_test_hard.csv'), index=False, lineterminator='\n')
    print(f'hard-negative test pairs: {len(h)}')

# verify no receptor overlap
tr_recs = set(train_pos.rec_seq) | set(train_neg.rec_seq)
te_recs = set(test_pos.rec_seq) | set(test_neg.rec_seq)
overlap = tr_recs & te_recs
print(f'\ntrain {len(train)} pairs / {train.prot_id.nunique()} receptors')
print(f'test  {len(test)} pairs / {test.prot_id.nunique()} receptors')
print(f'receptor overlap train/test (full seqs): {len(overlap)} (must be 0)')
print('wrote mvp_cpu/propedia_train_mtl.csv, propedia_test_mtl.csv')
