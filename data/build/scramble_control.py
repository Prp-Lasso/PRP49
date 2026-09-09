"""Proper scramble control: re-embed genuinely scrambled peptide sequences.

If the learned signal is sequence-specific, replacing positive peptides with
scrambled sequences (same composition, same protein) should drop AUC towards
chance. Negatives unchanged.
"""
import sys
import csv
import json
import random
import numpy as np
import torch
import torch.nn as nn

sys.stdout.reconfigure(encoding='utf-8')
torch.set_num_threads(4)


class MLP(nn.Module):
    def __init__(self, d=320):
        super().__init__()
        self.net = nn.Sequential(nn.Linear(d * 3, 256), nn.ReLU(), nn.Dropout(0.3),
                                 nn.Linear(256, 64), nn.ReLU(), nn.Dropout(0.3),
                                 nn.Linear(64, 1))

    def forward(self, x):
        return self.net(x).squeeze(-1)


def train(X, y, epochs=100, seed=0):
    torch.manual_seed(seed)
    np.random.seed(seed)
    m = MLP()
    opt = torch.optim.Adam(m.parameters(), lr=2e-3, weight_decay=1e-4)
    pos = (y == 1).sum()
    lossf = nn.BCEWithLogitsLoss(pos_weight=torch.tensor((len(y) - pos) / max(pos, 1)))
    Xt = torch.tensor(X, dtype=torch.float32)
    yt = torch.tensor(y, dtype=torch.float32)
    best, bad = None, 0
    for ep in range(epochs):
        m.train()
        perm = np.random.permutation(len(y))
        for i in range(0, len(y), 64):
            idx = perm[i:i + 64]
            loss = lossf(m(Xt[idx]), yt[idx])
            opt.zero_grad(); loss.backward(); opt.step()
        with torch.no_grad():
            val = lossf(m(Xt), yt).item()
        if best is None or val < best - 1e-4:
            best, bad = val, 0
        else:
            bad += 1
            if bad >= 20:
                break
    m.eval()
    return m


def metrics(scores, labels):
    pos = np.where(labels == 1)[0]
    neg = np.where(labels == 0)[0]
    total = 0.0
    for p in pos:
        for n in neg:
            total += (1.0 if scores[p] > scores[n] else 0.0) + (0.5 if scores[p] == scores[n] else 0.0)
    auc = total / (len(pos) * len(neg))
    order = np.argsort(-scores)
    ranked = labels[order]
    ap = float(np.sum((np.cumsum(ranked) / np.arange(1, len(ranked) + 1)) * ranked) / len(pos))
    return auc, ap


def main():
    rows = list(csv.DictReader(open('mvp_cpu/expanded_pairs.csv', encoding='utf-8')))
    z = np.load('mvp_cpu/emb_expanded.npz', allow_pickle=True)
    lookup = {k: v for k, v in zip(z['keys'], z['emb'])}

    # scrambled sequences for every unique peptide (pos AND neg pools)
    rng = random.Random(0)
    scram_map = {}
    all_peps = sorted(set(r['pep_seq'] for r in rows))
    for seq in all_peps:
        s = list(seq)
        rng.shuffle(s)
        s = ''.join(s)
        if s == seq:
            s = s[::-1]
        scram_map[seq] = s

    from transformers import AutoTokenizer, EsmModel
    tok = AutoTokenizer.from_pretrained('facebook/esm2_t6_8M_UR50D')
    model = EsmModel.from_pretrained('facebook/esm2_t6_8M_UR50D')
    model.eval()
    need = sorted(set(scram_map.values()))
    emb_scram = {}
    with torch.no_grad():
        for i in range(0, len(need), 16):
            batch = need[i:i + 16]
            enc = tok(batch, return_tensors='pt', padding=True, truncation=True, max_length=1022)
            h = model(**enc).last_hidden_state
            mask = enc['attention_mask'].unsqueeze(-1).float()
            mean = (h * mask).sum(1) / mask.sum(1)
            for j, s in enumerate(batch):
                emb_scram[s] = mean[j].numpy()
    print(f'embedded {len(need)} scrambled peptides', flush=True)

    for r in rows:
        r['label'] = int(r['label'])
        r['p'] = emb_scram[scram_map[r['pep_seq']]].astype(np.float32)
        r['q'] = lookup[f'Q|{r["prot_seq"]}'].astype(np.float32)
    X = np.stack([np.concatenate([r['p'], r['q'], r['p'] * r['q']]) for r in rows])
    y = np.array([r['label'] for r in rows])

    from sklearn.model_selection import StratifiedKFold
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    aucs, aps = [], []
    for tr, te in skf.split(X, y):
        m = train(X[tr], y[tr])
        with torch.no_grad():
            sc = m(torch.tensor(X[te], dtype=torch.float32)).numpy()
        a, p = metrics(sc, y[te])
        aucs.append(a); aps.append(p)
    print(f'scrambled-peptide (proper): AUC {np.mean(aucs):.3f}+-{np.std(aucs):.3f} AP {np.mean(aps):.3f}', flush=True)
    json.dump(dict(scrambled=dict(auc=float(np.mean(aucs)), ap=float(np.mean(aps)))),
              open('mvp_cpu/results_scramble.json', 'w'), indent=2)


if __name__ == '__main__':
    main()
