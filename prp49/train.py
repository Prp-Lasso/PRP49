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
from .data import PairDataset, collate_fn, load_pairs, load_alignment_labels, load_topo_annotations
from .model import InteractionPredictor
from .losses import binary_loss, alignment_loss, ranking_loss


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
        dual_head=cfg.model.dual_head,
        use_topo=bool(cfg.model.get('use_topo')))


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
        topo = batch['topo'].to(device) if batch.get('topo') is not None else None
        logits, att, aff = model(ids1, m1, ids2, m2, seq2_topo=topo)
        if cfg.train.label_type == 'energy':
            # regression mode: the label is a continuous affinity, so BCE is
            # meaningless (it goes negative for labels > 1). Use a robust
            # regression loss on the affinity head (or the main head if single).
            pred = aff if (cfg.model.dual_head and aff is not None) else logits
            loss = nn.functional.huber_loss(pred, lab, delta=1.0)
            rank_w = cfg.train.get('rank_loss_weight')
            if rank_w:
                rl = ranking_loss(pred, (lab > lab.median()).float())
                if rl is not None:
                    loss = loss + rank_w * rl
        else:
            loss = binary_loss(logits, lab, pos_weight=pos_weight)
            rank_w = cfg.train.get('rank_loss_weight')
            if rank_w:
                rl = ranking_loss(logits, lab)
                if rl is not None:
                    loss = loss + rank_w * rl
            if cfg.model.dual_head and aff is not None:
                loss = loss + nn.functional.mse_loss(aff, lab)
        if cfg.train.align_loss_weight:
            al = alignment_loss(att, batch.get('align'))
            if al is not None:
                loss = loss + cfg.train.align_loss_weight * al
                align_total += al.item()
                n_align += 1
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
        topo = batch['topo'].to(device) if batch.get('topo') is not None else None
        logits, _, _ = model(ids1, m1, ids2, m2, seq2_topo=topo)
        logits_list.append(logits.cpu().numpy())
        keys_list.extend(batch['keys'])
    return np.concatenate(logits_list), keys_list


def evaluate_scores(scores, labels, label_type='binary'):
    """Binary: (AUC, AP). Regression: (Spearman rho, -RMSE) - both "higher is better"."""
    if label_type == 'energy':
        from scipy.stats import spearmanr
        from sklearn.metrics import mean_squared_error
        rho = spearmanr(labels, scores).correlation
        if np.isnan(rho):
            rho = 0.0
        rmse = mean_squared_error(labels, scores) ** 0.5
        return float(rho), float(-rmse)
    from sklearn.metrics import roc_auc_score, average_precision_score
    try:
        auc = roc_auc_score(labels, scores)
    except ValueError:
        auc = float('nan')
    ap = average_precision_score(labels, scores)
    return auc, ap


def load_warmup_weights(model, ckpt_path, device):
    """Transfer shape-matching tensors (BAN + heads) from a warmup checkpoint.

    Warmup used dual ESM-2 35M; L2 uses 35M + LassoESM 650M, so encoder
    weights differ in shape and are skipped - only the fusion/head weights
    (the curriculum-learning payload) are transferred.
    """
    sd = torch.load(ckpt_path, map_location=device)
    own = model.state_dict()
    keep = {k: v for k, v in sd.items() if k in own and own[k].shape == v.shape}
    own.update(keep)
    model.load_state_dict(own)
    print(f'init from {ckpt_path}: {len(keep)}/{len(own)} tensors transferred '
          f'(e.g. {sorted(keep)[:3]})', flush=True)
    return model


def run_cv(cfg, device):
    """Cross-validated training per config; returns summary dict."""
    df = load_pairs(cfg.data.pairs_csv)
    tok = AutoTokenizer.from_pretrained(cfg.model.esm_model)
    folds = cfg.eval.cv_folds
    all_scores = np.zeros(len(df))
    if cfg.train.label_type == 'energy':
        ecol = cfg.data.get('energy_col') or 'energy'
        if ecol not in df.columns:
            raise KeyError(f"label_type=energy but no '{ecol}' column in {cfg.data.pairs_csv}")
        all_labels = df[ecol].values.astype(float)
        pos_weight = 1.0
    else:
        all_labels = df['label'].values.astype(float)
        pos = (all_labels == 1).sum()
        pos_weight = max((len(df) - pos) / max(pos, 1), 1.0)

    # alignment labels (L_align): only when weight > 0
    align_map = {}
    if cfg.train.align_loss_weight:
        align_map = load_alignment_labels(cfg.data.align_labels_json)
        _probe = PairDataset(df, tok, cfg.data.max_len, cfg.train.label_type, align_map=align_map,
                             energy_col=cfg.data.get('energy_col') or 'energy')
        n_hit = sum(1 for _, row in df.iterrows() if _probe._match_align(row) is not None)
        print(f'L_align: {len(align_map)} structural matrices, {n_hit}/{len(df)} pairs matched')

    # II2: topology annotations for the peptide side
    topo_map = {}
    if cfg.model.get('use_topo'):
        topo_map = load_topo_annotations(cfg.data.get('topo_csv'))
        n_topo = sum(1 for _, row in df.iterrows() if row['pep_seq'] in topo_map)
        print(f'topo channel: {len(topo_map)} annotated peptides in DB, '
              f'{n_topo}/{len(df)} pairs annotated')

    if cfg.eval.cv_mode == 'grouped' and cfg.data.group_col in df.columns:
        splitter = GroupKFold(n_splits=folds).split(np.arange(len(df)), groups=df[cfg.data.group_col].values)
    elif cfg.train.label_type == 'energy':
        from sklearn.model_selection import KFold
        splitter = KFold(n_splits=folds, shuffle=True,
                         random_state=cfg.train.seed).split(np.arange(len(df)))
    else:
        splitter = StratifiedKFold(n_splits=folds, shuffle=True, random_state=cfg.train.seed).split(
            np.zeros(len(df)), all_labels)

    fold_metrics = []
    for fold, (tr, te) in enumerate(splitter):
        set_seed(cfg.train.seed)
        model = build_model(cfg).to(device)
        if cfg.train.get('init_ckpt'):
            ckpt = cfg.train.init_ckpt
            if '{fold}' in ckpt:
                cand = ckpt.replace('{fold}', str(fold))
            else:
                cand = ckpt
            if os.path.exists(cand):
                model = load_warmup_weights(model, cand, device)
            else:
                print(f'WARN init_ckpt not found: {cand}', flush=True)
        optimizer = torch.optim.Adam(
            [p for p in model.parameters() if p.requires_grad],
            lr=cfg.train.lr, weight_decay=cfg.train.weight_decay)
        tr_loader = DataLoader(PairDataset(df.iloc[tr], tok, cfg.data.max_len, cfg.train.label_type,
                                           align_map=align_map, topo_map=topo_map,
                                           energy_col=cfg.data.get('energy_col') or 'energy'),
                               batch_size=cfg.train.batch_size, shuffle=True, collate_fn=collate_fn,
                               drop_last=len(tr) > cfg.train.batch_size)   # BAN's BatchNorm1d needs batch>1
        # Evaluation must never see propagated (inferred) labels: their supervision
        # came from family propagation, so scoring on them measures the propagation
        # rule, not the model. They stay in the TRAINING split only.
        eval_te = te
        if 'propagated' in df.columns:
            prop_mask = df.iloc[te]['propagated'].fillna(False).astype(bool).values
            if prop_mask.any() and (~prop_mask).sum() >= 10:
                eval_te = te[~prop_mask]
                if fold == 0:
                    print(f'  [eval] excluding {int(prop_mask.sum())} propagated rows '
                          f'from the test fold; {len(eval_te)} original rows kept')
        te_loader = DataLoader(PairDataset(df.iloc[eval_te], tok, cfg.data.max_len, cfg.train.label_type,
                                           topo_map=topo_map,
                                           energy_col=cfg.data.get('energy_col') or 'energy'),
                               batch_size=cfg.train.batch_size, shuffle=False, collate_fn=collate_fn)
        best_auc, bad, best_state = -1.0, 0, None
        for ep in range(cfg.train.epochs):
            loss, n_align = train_one(model, tr_loader, optimizer, pos_weight, cfg, device)
            sc, _ = predict(model, te_loader, device)
            auc, ap = evaluate_scores(sc, all_labels[eval_te], cfg.train.label_type)
            if auc > best_auc:
                best_auc, bad = auc, 0
                best_state = copy.deepcopy(model.state_dict())
            else:
                bad += 1
                if bad >= cfg.train.patience:
                    break
            if ep % 10 == 0:
                m1, m2 = ('Spearman', 'RMSE') if cfg.train.label_type == 'energy' else ('AUC', 'AP')
                v2 = -ap if cfg.train.label_type == 'energy' else ap
                print(f'  fold {fold} ep {ep}: loss {loss:.4f} val {m1} {auc:.3f} '
                      f'{m2} {v2:.3f} (align batches {n_align})', flush=True)
        model.load_state_dict(best_state)
        ckpt_dir = cfg.paths.checkpoint_dir
        os.makedirs(ckpt_dir, exist_ok=True)
        torch.save(best_state, os.path.join(ckpt_dir, f'fold{fold}_best.pt'))
        sc, _ = predict(model, te_loader, device)
        auc, ap = evaluate_scores(sc, all_labels[eval_te], cfg.train.label_type)
        fold_metrics.append(dict(auc=float(auc), ap=float(ap)))
        # sc has len(eval_te) when propagated rows were excluded - must scatter with
        # the SAME index, otherwise numpy raises a shape mismatch (this crashed the
        # first propagation run at the end of fold 0).
        all_scores[eval_te] = sc
        m1, m2 = ('Spearman', 'RMSE') if cfg.train.label_type == 'energy' else ('AUC', 'AP')
        v2 = -ap if cfg.train.label_type == 'energy' else ap
        print(f'fold {fold}: {m1} {auc:.3f} {m2} {v2:.3f}', flush=True)

    # Propagated rows are excluded from evaluation, so they never receive a score
    # and would sit at their initial 0.0 - the overall metric must use the same
    # subset, otherwise it is computed partly on placeholder values.
    if 'propagated' in df.columns:
        eval_mask = ~df['propagated'].fillna(False).astype(bool).values
    else:
        eval_mask = np.ones(len(df), dtype=bool)
    overall_auc, overall_ap = evaluate_scores(all_scores[eval_mask], all_labels[eval_mask],
                                              cfg.train.label_type)
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
    parser.add_argument('--init_ckpt', default=None,
                        help='warmup checkpoint to initialise BAN/heads (curriculum S1->S2); '
                             'may contain {fold}')
    parser.add_argument('--out_suffix', default=None, help='override output_dir/checkpoint_dir suffix')
    args = parser.parse_args()
    cfg = Config.from_yaml(args.config)
    cfg._base_dir = os.path.dirname(os.path.abspath(args.config))
    if args.init_ckpt:
        cfg._d.setdefault('train', {})['init_ckpt'] = args.init_ckpt
    if args.out_suffix:
        cfg._d['paths']['output_dir'] = f'./runs_{args.out_suffix}'
        cfg._d['paths']['checkpoint_dir'] = f'./runs_{args.out_suffix}/checkpoints'
    device = (torch.device('cuda' if torch.cuda.is_available() else 'cpu')
              if args.device == 'auto' else torch.device(args.device))
    print(f'device: {device}')
    run_cv(cfg, device)


if __name__ == '__main__':
    main()
