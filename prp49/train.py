"""Training loop with CV support, early stopping, alignment loss hook."""
import os
import json
import copy
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader
from transformers import AutoTokenizer
from sklearn.model_selection import StratifiedKFold, GroupKFold

from .config import Config
from .data import PairDataset, collate_fn, load_pairs, load_alignment_labels
from .model import InteractionPredictor
from .losses import binary_loss, alignment_loss


def set_seed(seed):
    torch.manual_seed(seed)
    np.random.seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def build_model(cfg):
    return InteractionPredictor(
        esm_model_name=cfg.model.esm_model,
        lassoesm_model_name=cfg.model.lassoesm_model,
        h_dim=cfg.model.h_dim, n_heads=cfg.model.n_heads, k=cfg.model.k,
        dropout=cfg.model.dropout, ban_dropout=cfg.model.ban_dropout,
        freeze_esm=cfg.model.freeze_esm,
        lassoesm_unfreeze=cfg.model.lassoesm_unfreeze,
        dual_head=cfg.model.dual_head)


def train_one(model, loader, optimizer, pos_weight, cfg, device):
    model.train()
    total = 0.0
    align_total = 0.0
    n = 0
    n_align = 0
    for i, batch in enumerate(loader):
        ids1 = batch['seq1_ids'].to(device)
        m1 = batch['seq1_mask'].to(device)
        ids2 = batch['seq2_ids'].to(device)
        m2 = batch['seq2_mask'].to(device)
        lab = batch['label'].to(device)
        logits, att, aff = model(ids1, m1, ids2, m2)
        loss = binary_loss(logits, lab, pos_weight=pos_weight)
        if cfg.train.align_loss_weight:
            al = alignment_loss(att, batch.get('align'))
            if al is not None:
                loss = loss + cfg.train.align_loss_weight * al
                align_total += al.item()
                n_align += 1
        if cfg.model.dual_head and aff is not None:
            loss = loss + nn.functional.mse_loss(aff, lab)
        optimizer.zero_grad()
        loss.backward()
        nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
        optimizer.step()
        total += loss.item() * len(lab)
        n += len(lab)
    return total / max(n, 1), n_align


@torch.no_grad()
def predict(model, loader, device):
    model.eval()
    logits_list, keys_list = [], []
    for batch in loader:
        ids1 = batch['seq1_ids'].to(device)
        m1 = batch['seq1_mask'].to(device)
        ids2 = batch['seq2_ids'].to(device)
        m2 = batch['seq2_mask'].to(device)
        logits, _, _ = model(ids1, m1, ids2, m2)
        logits_list.append(logits.cpu().numpy())
        keys_list.extend(batch['keys'])
    return np.concatenate(logits_list), keys_list


def evaluate_scores(scores, labels):
    from sklearn.metrics import roc_auc_score, average_precision_score
    try:
        auc = roc_auc_score(labels, scores)
    except ValueError:
        auc = float('nan')
    ap = average_precision_score(labels, scores)
    return auc, ap


def run_cv(cfg, device):
    """Cross-validated training per config; returns summary dict."""
    df = load_pairs(cfg.data.pairs_csv)
    tok = AutoTokenizer.from_pretrained(cfg.model.esm_model)
    folds = cfg.eval.cv_folds
    all_scores = np.zeros(len(df))
    all_labels = df['label'].values.astype(float)

    pos = (all_labels == 1).sum()
    pos_weight = max((len(df) - pos) / max(pos, 1), 1.0)

    # alignment labels (L_align): only when weight > 0
    align_map = {}
    if cfg.train.align_loss_weight:
        align_map = load_alignment_labels(cfg.data.align_labels_json)
        _probe = PairDataset(df, tok, cfg.data.max_len, cfg.train.label_type, align_map=align_map)
        n_hit = sum(1 for _, row in df.iterrows() if _probe._match_align(row) is not None)
        print(f'L_align: {len(align_map)} structural matrices, {n_hit}/{len(df)} pairs matched')

    if cfg.eval.cv_mode == 'grouped' and cfg.data.group_col in df.columns:
        splitter = GroupKFold(n_splits=folds).split(np.arange(len(df)), groups=df[cfg.data.group_col].values)
    else:
        splitter = StratifiedKFold(n_splits=folds, shuffle=True, random_state=cfg.train.seed).split(
            np.zeros(len(df)), all_labels)

    fold_metrics = []
    for fold, (tr, te) in enumerate(splitter):
        set_seed(cfg.train.seed)
        model = build_model(cfg).to(device)
        optimizer = torch.optim.Adam(
            [p for p in model.parameters() if p.requires_grad],
            lr=cfg.train.lr, weight_decay=cfg.train.weight_decay)
        tr_loader = DataLoader(PairDataset(df.iloc[tr], tok, cfg.data.max_len, cfg.train.label_type,
                                           align_map=align_map),
                               batch_size=cfg.train.batch_size, shuffle=True, collate_fn=collate_fn)
        te_loader = DataLoader(PairDataset(df.iloc[te], tok, cfg.data.max_len, cfg.train.label_type),
                               batch_size=cfg.train.batch_size, shuffle=False, collate_fn=collate_fn)
        best_auc, bad, best_state = -1.0, 0, None
        for ep in range(cfg.train.epochs):
            loss, n_align = train_one(model, tr_loader, optimizer, pos_weight, cfg, device)
            sc, _ = predict(model, te_loader, device)
            auc, ap = evaluate_scores(sc, all_labels[te])
            if auc > best_auc:
                best_auc, bad = auc, 0
                best_state = copy.deepcopy(model.state_dict())
            else:
                bad += 1
                if bad >= cfg.train.patience:
                    break
            if ep % 10 == 0:
                print(f'  fold {fold} ep {ep}: loss {loss:.4f} val AUC {auc:.3f} AP {ap:.3f}'
                      f' (align batches {n_align})', flush=True)
        model.load_state_dict(best_state)
        ckpt_dir = cfg.paths.checkpoint_dir
        os.makedirs(ckpt_dir, exist_ok=True)
        torch.save(best_state, os.path.join(ckpt_dir, f'fold{fold}_best.pt'))
        sc, _ = predict(model, te_loader, device)
        auc, ap = evaluate_scores(sc, all_labels[te])
        fold_metrics.append(dict(auc=float(auc), ap=float(ap)))
        all_scores[te] = sc
        print(f'fold {fold}: AUC {auc:.3f} AP {ap:.3f}', flush=True)

    overall_auc, overall_ap = evaluate_scores(all_scores, all_labels)
    summary = dict(
        folds=fold_metrics,
        auc_mean=float(np.mean([m['auc'] for m in fold_metrics])),
        auc_std=float(np.std([m['auc'] for m in fold_metrics])),
        ap_mean=float(np.mean([m['ap'] for m in fold_metrics])),
        overall_auc=float(overall_auc), overall_ap=float(overall_ap))
    os.makedirs(cfg.paths.output_dir, exist_ok=True)
    with open(os.path.join(cfg.paths.output_dir, 'cv_summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))
    return summary


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', default='config.yaml')
    parser.add_argument('--device', default='auto')
    args = parser.parse_args()
    cfg = Config.from_yaml(args.config)
    cfg._base_dir = os.path.dirname(os.path.abspath(args.config))
    device = (torch.device('cuda' if torch.cuda.is_available() else 'cpu')
              if args.device == 'auto' else torch.device(args.device))
    print(f'device: {device}')
    run_cv(cfg, device)


if __name__ == '__main__':
    main()
