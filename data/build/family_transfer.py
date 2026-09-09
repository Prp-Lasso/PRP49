"""Family-internal transfer test: train on one variant, test on same-family variant.

Pairs: (train_variant, test_variant, family)
  RES701-1 -> RES701-3 (same target ETB)
  MCoTI-II x bovine trypsin -> human trypsin (homolog target)
  Anantin -> Anantin-v2
  Cilengitide-b3 (ITGB3) -> Cilengitide-av (ITGAV)
  Capistruin -> Capi-var1 -> Capi-var2 (RNAP)
Each test: positives of test variant + clean negatives; model trained on rest.
Report rank of each positive among negatives (percentile) and hit@k.
"""
import sys
import csv
import numpy as np
import torch
import torch.nn as nn

sys.stdout.reconfigure(encoding='utf-8')

SEEDS = [0, 1, 2, 3, 4]


def load_data():
    rows = list(csv.DictReader(open('mvp_cpu/pairs.csv', encoding='utf-8')))
    z = np.load('mvp_cpu/emb.npz', allow_pickle=True)
    lookup = {k: v for k, v in zip(z['keys'], z['emb'])}
    for r in rows:
        r['p'] = lookup[f'P|{r["pep_seq"]}'].astype(np.float32)
        r['q'] = lookup[f'Q|{r["prot_seq"]}'].astype(np.float32)
        r['label'] = int(r['label'])
    return rows


class MLP(nn.Module):
    def __init__(self, d=320):
        super().__init__()
        self.net = nn.Sequential(nn.Linear(d * 3, 256), nn.ReLU(), nn.Dropout(0.3),
                                 nn.Linear(256, 64), nn.ReLU(), nn.Dropout(0.3),
                                 nn.Linear(64, 1))

    def forward(self, x):
        return self.net(x).squeeze(-1)


def train(X, y, seed):
    torch.manual_seed(seed)
    np.random.seed(seed)
    m = MLP()
    opt = torch.optim.Adam(m.parameters(), lr=2e-3, weight_decay=1e-4)
    pos = (y == 1).sum()
    lossf = nn.BCEWithLogitsLoss(pos_weight=torch.tensor((len(y) - pos) / max(pos, 1)))
    Xt = torch.tensor(X, dtype=torch.float32)
    yt = torch.tensor(y, dtype=torch.float32)
    best, bad = None, 0
    for ep in range(300):
        m.train()
        perm = np.random.permutation(len(y))
        for i in range(0, len(y), 32):
            idx = perm[i:i + 32]
            loss = lossf(m(Xt[idx]), yt[idx])
            opt.zero_grad(); loss.backward(); opt.step()
        with torch.no_grad():
            val = lossf(m(Xt), yt).item()
        if best is None or val < best - 1e-4:
            best, bad = val, 0
        else:
            bad += 1
            if bad >= 40:
                break
    m.eval()
    return m


def feats(r):
    return np.concatenate([r['p'], r['q'], r['p'] * r['q']])


def main():
    rows = load_data()
    pos = {r['pep_id']: r for r in rows if r['label'] == 1}
    negs = [r for r in rows if r['label'] == 0]

    tests = [
        ('RES701-1 -> RES701-3 (same target ETB)', 'RES701-1', 'RES701-3'),
        ('RES701-3 -> RES701-1', 'RES701-3', 'RES701-1'),
        ('MCoTI-II x bovineTRY -> humanTRY', 'MCoTI-II', 'MCoTI-II-h'),
        ('MCoTI-II-h -> bovineTRY', 'MCoTI-II-h', 'MCoTI-II'),
        ('Anantin -> Anantin-v2', 'Anantin', 'Anantin-v2'),
        ('Anantin-v2 -> Anantin', 'Anantin-v2', 'Anantin'),
        ('Cilengitide-ITGB3 -> ITGAV', 'Cilengitide-b3', 'Cilengitide-av'),
        ('Cilengitide-ITGAV -> ITGB3', 'Cilengitide-av', 'Cilengitide-b3'),
        ('Capistruin -> Capi-var1 (RNAP)', 'Capistruin', 'Capi-var1'),
        ('Capi-var1 -> Capi-var2', 'Capi-var1', 'Capi-var2'),
    ]
    print(f'{"transfer":44s} {"test_pos_rank_among_40negs":>30s}  hit@5')
    for name, tr_id, te_id in tests:
        tr = pos[tr_id]
        te = pos[te_id]
        train_rows = [r for r in rows if r['pep_id'] != te_id and r['label'] == 1]
        # train: all positives except test variant + negatives not sharing test peptide
        te_neg = [r for r in negs if r['pep_seq'] != te['pep_seq']][:40]
        train_rows += [r for r in negs if r not in te_neg]
        X = np.stack([feats(r) for r in train_rows])
        y = np.array([r['label'] for r in train_rows])
        Xt = np.stack([feats(te)] + [feats(r) for r in te_neg])
        scores = np.mean([train(X, y, s)(torch.tensor(Xt, dtype=torch.float32)).detach().numpy() for s in SEEDS], axis=0)
        pos_score = scores[0]
        neg_scores = scores[1:]
        # percentile of pos among negs
        pct = np.mean(neg_scores < pos_score)
        hit5 = float(pos_score >= np.sort(neg_scores)[-5])
        print(f'{name:44s} {pct * 100:26.1f}%  {hit5:>5.0f}')

    # control: random positive (should be ~similar distribution)
    print('\ncontrol (unrelated positive as test):')
    ctrl = [('p53pep', 'Lassomycin'), ('p53pep', 'PB1m7'), ('p53pep', 'MccJ25')]
    for a, b in ctrl:
        te = pos[b]
        train_rows = [r for r in rows if r['pep_id'] != b and r['label'] == 1]
        te_neg = [r for r in negs if r['pep_seq'] != te['pep_seq']][:40]
        train_rows += [r for r in negs if r not in te_neg]
        X = np.stack([feats(r) for r in train_rows])
        y = np.array([r['label'] for r in train_rows])
        Xt = np.stack([feats(te)] + [feats(r) for r in te_neg])
        scores = np.mean([train(X, y, s)(torch.tensor(Xt, dtype=torch.float32)).detach().numpy() for s in SEEDS], axis=0)
        pct = np.mean(scores[1:] < scores[0])
        print(f'  {a:20s} -> {b:14s} {pct * 100:5.1f}%')


if __name__ == '__main__':
    main()
