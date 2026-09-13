"""Probability calibration + decision-threshold selection.

Problem: 29% of screening predictions exceed 0.9 even for pairs known to be
inactive (wild-type MccJ25 on integrins). Raw logits are therefore not usable as
probabilities, and a fixed 0.5 cut-off is meaningless.

Method (standard, cheap):
  1. score a held-out set with peptide-grouped folds (no leakage),
  2. fit a single temperature T minimising NLL (Guo et al. 2017),
  3. report ECE before/after and pick an operating threshold that maximises F1
     (plus the threshold that reaches >=90% precision, for screening use).

Output: JSON with T, ECE before/after, thresholds, and calibrated metrics.
"""
import argparse
import json
import os

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from sklearn.metrics import f1_score, precision_recall_curve, roc_auc_score
from torch.utils.data import DataLoader
from transformers import AutoTokenizer

from .config import Config
from .data import PairDataset, collate_fn
from .model import InteractionPredictor


def ece(scores, labels, bins=15):
    """Expected calibration error."""
    edges = np.linspace(0, 1, bins + 1)
    total = len(scores)
    err = 0.0
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (scores >= lo) & (scores < hi)
        if m.sum() == 0:
            continue
        conf = scores[m].mean()
        acc = labels[m].mean()
        err += (m.sum() / total) * abs(conf - acc)
    return float(err)


@torch.no_grad()
def collect_scores(model, loader, device):
    model.eval()
    logits, labels = [], []
    for batch in loader:
        out = model(batch['seq1_ids'].to(device), batch['seq1_mask'].to(device),
                    batch['seq2_ids'].to(device), batch['seq2_mask'].to(device))
        logits.append(out[0].detach().cpu())
        labels.append(batch['label'])
    return torch.cat(logits), torch.cat(labels)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--config', required=True)
    ap.add_argument('--ckpt_dir', required=True)
    ap.add_argument('--out', default='calibration.json')
    ap.add_argument('--device', default='cuda')
    args = ap.parse_args()

    cfg = Config.from_yaml(args.config)
    cfg._base_dir = os.path.dirname(os.path.abspath(args.config))
    tok = AutoTokenizer.from_pretrained(cfg.model.esm_model)
    df = pd.read_csv(cfg.data.pairs_csv)
    print(f'pairs: {len(df)} | labels {df.label.value_counts().to_dict()}')

    # peptide-grouped folds, mirroring the training protocol
    groups = df[cfg.data.group_col].values
    uniq = pd.unique(groups)
    rng = np.random.default_rng(cfg.train.get('seed', 42))
    rng.shuffle(uniq)
    folds = np.array_split(uniq, cfg.eval.cv_folds)

    ckpts = sorted(f for f in os.listdir(args.ckpt_dir) if f.endswith('.pt') and f.startswith('fold'))
    print(f'checkpoints: {ckpts}')

    all_logits, all_labels = [], []
    for fi, test_groups in enumerate(folds):
        if fi >= len(ckpts):
            break
        te = df[df[cfg.data.group_col].isin(test_groups)]
        if len(te) == 0 or te.label.nunique() < 2:
            continue
        model = InteractionPredictor(
            esm_model_name=cfg.model.esm_model, lassoesm_model_name=cfg.model.lassoesm_model,
            h_dim=cfg.model.h_dim, n_heads=cfg.model.n_heads, k=cfg.model.k,
            dropout=cfg.model.dropout, ban_dropout=cfg.model.ban_dropout,
            freeze_esm=False, lassoesm_unfreeze=0, dual_head=False, use_topo=False).to(args.device)
        sd = torch.load(os.path.join(args.ckpt_dir, ckpts[fi]), map_location=args.device)
        model.load_state_dict(sd.get('state_dict', sd), strict=False)
        loader = DataLoader(PairDataset(te, tok, cfg.data.max_len, 'binary',
                                        energy_col='energy'),
                            batch_size=24, collate_fn=collate_fn)
        lg, lb = collect_scores(model, loader, args.device)
        print(f'  fold {fi}: n={len(lb)} pos={int(lb.sum())}')
        all_logits.append(lg); all_labels.append(lb)
        del model
        torch.cuda.empty_cache()

    logits = torch.cat(all_logits).float()
    labels = torch.cat(all_labels).float().numpy()
    print(f'collected {len(labels)} held-out predictions, positives {int(labels.sum())}')

    # --- temperature scaling on the logit
    T = nn.Parameter(torch.ones(1))
    opt = torch.optim.LBFGS([T], lr=0.1, max_iter=200)
    bce = nn.BCEWithLogitsLoss()

    def closure():
        opt.zero_grad()
        loss = bce(logits / T.clamp(min=1e-3), torch.from_numpy(labels).float())
        loss.backward()
        return loss
    opt.step(closure)

    raw_p = torch.sigmoid(logits).numpy()
    cal_p = torch.sigmoid(logits / T.clamp(min=1e-3)).detach().numpy()
    results = {
        'n': int(len(labels)), 'positives': int(labels.sum()),
        'temperature': float(T.item()),
        'auc_raw': float(roc_auc_score(labels, raw_p)),
        'auc_calibrated': float(roc_auc_score(labels, cal_p)),
        'ece_raw': ece(raw_p, labels),
        'ece_calibrated': ece(cal_p, labels),
        'mean_prob_raw': float(raw_p.mean()), 'mean_prob_calibrated': float(cal_p.mean()),
        'base_rate': float(labels.mean()),
    }
    # thresholds on calibrated probabilities
    prec, rec, thr = precision_recall_curve(labels, cal_p)
    f1 = 2 * prec * rec / np.clip(prec + rec, 1e-9, None)
    best = int(np.nanargmax(f1[:-1])) if len(thr) else 0
    results['threshold_max_f1'] = float(thr[best]) if len(thr) else 0.5
    results['f1_at_max'] = float(f1[best]) if len(f1) else float('nan')
    ok = np.where(prec[:-1] >= 0.90)[0]
    results['threshold_precision_90'] = float(thr[ok[0]]) if len(ok) else None
    results['precision_at_90_thr'] = float(prec[ok[0]]) if len(ok) else None
    results['f1_at_090'] = float(f1[ok[0]]) if len(ok) else None

    print(json.dumps(results, indent=2))
    with open(args.out, 'w') as f:
        json.dump(results, f, indent=2)
    print(f'wrote {args.out}')


if __name__ == '__main__':
    main()
