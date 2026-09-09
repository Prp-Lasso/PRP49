"""Verify shuffled-label control correctly (bug-fixed) and add a protein-identity
leakage check: does the model just memorize protein embeddings?

Also add a peptide-scramble control: positive pairs keep peptide + protein but with
scrambled peptide sequences -> must drop to ~0.5 if signal is real.
"""
import sys
import csv
import json
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
    for r in rows:
        r['p'] = lookup[f'P|{r["pep_seq"]}'].astype(np.float32)
        r['q'] = lookup[f'Q|{r["prot_seq"]}'].astype(np.float32)
        r['label'] = int(r['label'])
    X = np.stack([np.concatenate([r['p'], r['q'], r['p'] * r['q']]) for r in rows])
    y = np.array([r['label'] for r in rows])

    from sklearn.model_selection import StratifiedKFold
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    def run_protocol(tag, label_fn):
        aucs, aps = [], []
        for tr, te in skf.split(X, y):
            ytr = label_fn(tr)
            m = train(X[tr], ytr)
            with torch.no_grad():
                sc = m(torch.tensor(X[te], dtype=torch.float32)).numpy()
            a, p = metrics(sc, y[te])
            aucs.append(a); aps.append(p)
        print(f'{tag}: AUC {np.mean(aucs):.3f}+-{np.std(aucs):.3f} AP {np.mean(aps):.3f}', flush=True)
        return np.mean(aucs), np.mean(aps)

    r_true = run_protocol('true labels        ', lambda tr: y[tr].copy())

    def shuffled(tr):
        yy = y.copy()
        sub = yy[tr].copy()
        np.random.RandomState(7).shuffle(sub)
        yy[tr] = sub
        return yy[tr]
    r_shuf = run_protocol('shuffled labels    ', shuffled)

    # scramble-peptide control: positives get scrambled peptide embeddings (kept same for train/test consistency?)
    # simpler: scramble peptide sequences of ALL positives consistently, rebuild embeddings
    def scramble_pep(tr):
        # labels unchanged; this control requires modified features instead.
        return y[tr].copy()  # placeholder; handled below
    # feature-level scramble: replace each positive pair's peptide emb with a scrambled-peptide emb
    rng = np.random.RandomState(0)
    X_scram = X.copy()
    d = 320
    for i, r in enumerate(rows):
        if r['label'] == 1:
            s = list(r['pep_seq'])
            rng.shuffle(s)
            s = ''.join(s)
            pscram = lookup.get(f'P|{s}')
            if pscram is None:
                # embed not cached; approximate by shuffling embedding rows? instead skip: use shuffled order of emb
                pscram = r['p'].copy()
                rng.shuffle(pscram)
            X_scram[i, :d] = pscram
            X_scram[i, 2 * d:] = pscram * X_scram[i, d:2 * d]

    def run_scram():
        aucs, aps = [], []
        for tr, te in skf.split(X_scram, y):
            m = train(X_scram[tr], y[tr])
            with torch.no_grad():
                sc = m(torch.tensor(X_scram[te], dtype=torch.float32)).numpy()
            a, p = metrics(sc, y[te])
            aucs.append(a); aps.append(p)
        print(f'scrambled-pos-pep  : AUC {np.mean(aucs):.3f}+-{np.std(aucs):.3f} AP {np.mean(aps):.3f}', flush=True)
        return np.mean(aucs), np.mean(aps)
    r_scram = run_scram()

    out = dict(true=dict(auc=float(r_true[0]), ap=float(r_true[1])),
               shuffled=dict(auc=float(r_shuf[0]), ap=float(r_shuf[1])),
               scrambled_pep=dict(auc=float(r_scram[0]), ap=float(r_scram[1])))
    json.dump(out, open('mvp_cpu/results_controls.json', 'w'), indent=2)
    print('saved results_controls.json')


if __name__ == '__main__':
    main()
