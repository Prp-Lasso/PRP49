"""Multi-task training: one shared BAN trunk supervised on BOTH domains.

Why: the two domains proved mutually untransferable when trained separately -
the ranking head learned on ordinary peptides actively hurts on lasso peptides
(known-pair mean rank 23 -> 71), and the lasso classifier is at chance on ordinary
peptide-protein complexes (AUC 0.47). Multi-task training is the direct test of
whether a shared trunk can learn anything transferable.

Setup
  task A (in-domain) : lasso peptide x human target binary binding (train_pairs_hard.csv)
  task B (auxiliary) : ordinary peptide x protein binary binding (propedia_train_mtl.csv);
                       evaluated on receptors held out by full-sequence split
Loss: BCE(A) + lambda * BCE(B). Both tasks use the identical input format, so one
forward pass serves both and no architecture change is needed.

Reports per fold: in-domain AUC (peptide-grouped CV) and cross-domain AUC on the
held-out Propedia receptors - the number that must beat the 0.47 single-task baseline.
"""
import argparse
import json
import os
import time

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from sklearn.metrics import roc_auc_score
from torch.utils.data import DataLoader
from transformers import AutoTokenizer

from .config import Config
from .data import PairDataset, collate_fn
from .losses import alignment_loss
from .model import InteractionPredictor


def grouped_folds(df, group_col, k, seed=42):
    groups = pd.unique(df[group_col].values)
    rng = np.random.default_rng(seed)
    rng.shuffle(groups)
    return np.array_split(groups, k)


def loader_for(df, tok, max_len, batch, shuffle):
    return DataLoader(PairDataset(df, tok, max_len, 'binary'),
                      batch_size=batch, shuffle=shuffle, collate_fn=collate_fn,
                      drop_last=shuffle and len(df) > batch)


@torch.no_grad()
def evaluate(model, loader, device):
    model.eval()
    scores, labels = [], []
    for batch in loader:
        logits, _, _ = model(batch['seq1_ids'].to(device), batch['seq1_mask'].to(device),
                             batch['seq2_ids'].to(device), batch['seq2_mask'].to(device))
        scores.append(logits.detach().cpu())
        labels.append(batch['label'])
    s = torch.cat(scores).float().numpy()
    y = torch.cat(labels).float().numpy()
    if len(np.unique(y)) < 2:
        return float('nan')
    return float(roc_auc_score(y, s))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--config', required=True)
    ap.add_argument('--aux_csv', required=True, help='ordinary-peptide training set')
    ap.add_argument('--aux_test_csv', required=True, help='held-out receptors')
    ap.add_argument('--lambda_aux', type=float, default=0.5)
    ap.add_argument('--device', default='cuda')
    ap.add_argument('--out', default='mtl_results.json')
    args = ap.parse_args()

    cfg = Config.from_yaml(args.config)
    cfg._base_dir = os.path.dirname(os.path.abspath(args.config))
    torch.manual_seed(cfg.train.seed)
    np.random.seed(cfg.train.seed)

    tok = AutoTokenizer.from_pretrained(cfg.model.esm_model)
    dom = pd.read_csv(cfg.data.pairs_csv)
    aux = pd.read_csv(args.aux_csv)
    aux_te = pd.read_csv(args.aux_test_csv)
    print(f'in-domain pairs {len(dom)} | aux train {len(aux)} | aux test {len(aux_te)}')
    print(f'lambda_aux = {args.lambda_aux}')

    aux_te_loader = loader_for(aux_te, tok, cfg.data.max_len, 24, False)
    folds = grouped_folds(dom, cfg.data.group_col, cfg.eval.cv_folds, cfg.train.seed)
    results = []

    for fi, test_groups in enumerate(folds):
        tr = dom[~dom[cfg.data.group_col].isin(test_groups)]
        te = dom[dom[cfg.data.group_col].isin(test_groups)]
        if len(te) == 0 or te.label.nunique() < 2:
            continue
        t0 = time.time()
        model = InteractionPredictor(
            esm_model_name=cfg.model.esm_model, lassoesm_model_name=cfg.model.lassoesm_model,
            h_dim=cfg.model.h_dim, n_heads=cfg.model.n_heads, k=cfg.model.k,
            dropout=cfg.model.dropout, ban_dropout=cfg.model.ban_dropout,
            freeze_esm=False, lassoesm_unfreeze=cfg.model.get('lassoesm_unfreeze', 0),
            dual_head=False, use_topo=False).to(args.device)
        opt = torch.optim.Adam([p for p in model.parameters() if p.requires_grad],
                               lr=cfg.train.lr, weight_decay=cfg.train.weight_decay)
        dom_loader = loader_for(tr, tok, cfg.data.max_len, cfg.train.batch_size, True)
        aux_loader = loader_for(aux, tok, cfg.data.max_len, cfg.train.batch_size, True)
        te_loader = loader_for(te, tok, cfg.data.max_len, 24, False)
        bce = nn.BCEWithLogitsLoss()

        best = {'dom': -1.0}
        for ep in range(cfg.train.epochs):
            model.train()
            # Alternate domains rather than running both every step: the in-domain
            # target sequences are up to 1024 tokens, and two batches per step
            # exhausted the A100 (CUDA OOM). Alternating halves peak memory and
            # compute while still updating the shared trunk from both tasks.
            aux_it = iter(aux_loader)
            tot, nb = 0.0, 0
            for step, batch in enumerate(dom_loader):
                opt.zero_grad()
                if step % 2 == 0:
                    lg, att, _ = model(batch['seq1_ids'].to(args.device), batch['seq1_mask'].to(args.device),
                                       batch['seq2_ids'].to(args.device), batch['seq2_mask'].to(args.device))
                    loss = bce(lg, batch['label'].to(args.device))
                else:
                    try:
                        ab = next(aux_it)
                    except StopIteration:
                        aux_it = iter(aux_loader)
                        ab = next(aux_it)
                    alg, alt, _ = model(ab['seq1_ids'].to(args.device), ab['seq1_mask'].to(args.device),
                                        ab['seq2_ids'].to(args.device), ab['seq2_mask'].to(args.device))
                    loss = args.lambda_aux * bce(alg, ab['label'].to(args.device))
                loss.backward()
                torch.nn.utils.clip_grad_norm_([p for p in model.parameters() if p.requires_grad], 5.0)
                opt.step()
                tot += float(loss); nb += 1

            dom_auc = evaluate(model, te_loader, args.device)
            aux_auc = evaluate(model, aux_te_loader, args.device)
            if dom_auc == dom_auc and dom_auc > best['dom']:
                best = {'dom': dom_auc, 'aux': aux_auc, 'epoch': ep, 'loss': tot / max(nb, 1)}
            if ep % 2 == 0:
                print(f'  fold {fi} ep {ep}: loss {tot/max(nb,1):.4f} '
                      f'in-domain AUC {dom_auc:.3f} | cross-domain AUC {aux_auc:.3f}', flush=True)
        print(f'fold {fi}: in-domain {best["dom"]:.3f} | cross-domain {best.get("aux", float("nan")):.3f} '
              f'(ep {best.get("epoch")}, {time.time()-t0:.0f}s)', flush=True)
        results.append(best)
        ck = os.path.join(cfg.paths.checkpoint_dir or 'runs_mtl/checkpoints', f'fold{fi}_best.pt')
        os.makedirs(os.path.dirname(ck), exist_ok=True)
        torch.save({'state_dict': model.state_dict(), 'metrics': best}, ck)
        del model
        torch.cuda.empty_cache()

    dom = [r['dom'] for r in results if r['dom'] == r['dom']]
    auxv = [r.get('aux') for r in results if r.get('aux') == r.get('aux')]
    summary = {'folds': results,
               'in_domain_auc_mean': float(np.mean(dom)) if dom else None,
               'in_domain_auc_std': float(np.std(dom)) if dom else None,
               'cross_domain_auc_mean': float(np.mean(auxv)) if auxv else None,
               'cross_domain_auc_std': float(np.std(auxv)) if auxv else None,
               'lambda_aux': args.lambda_aux}
    print('\n' + json.dumps(summary, indent=2))
    with open(args.out, 'w') as f:
        json.dump(summary, f, indent=2)
    print('MTL_DONE')


if __name__ == '__main__':
    main()
