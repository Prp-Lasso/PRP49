"""Evaluation: held-out metrics, ring-dependence ablation, baselines."""
import json
import os
import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader
from transformers import AutoTokenizer

from .config import Config
from .data import PairDataset, collate_fn, load_pairs
from .model import InteractionPredictor
from .train import predict, evaluate_scores, set_seed


def load_ckpt(model, ckpt_path, device):
    model.load_state_dict(torch.load(ckpt_path, map_location=device))
    return model


def ring_ablation(cfg, model, tok, device, df_base, ablation_csv):
    """Compare scores of (lasso, ring-open, distal-control) triples.

    Requires ablation_set.csv with columns orig/seq/kind; and a base pairs
    dataframe containing the original peptide x target protein sequence.
    """
    abl = pd.read_csv(ablation_csv)
    rows = []
    for orig, grp in abl.groupby('orig'):
        base = df_base[df_base['pep_seq'] == orig]
        # map orig id -> pairs csv pep_id
        base2 = df_base[df_base['pep_id'].astype(str).str.startswith(orig)]
        base = base2 if len(base2) else base
        if not len(base):
            continue
        prot_seq = base.iloc[0]['prot_seq']
        base_pep = base.iloc[0]['pep_seq']
        scores = {}
        # score the parent lasso peptide itself (baseline for deltas)
        dfb = pd.DataFrame([{'prot_seq': prot_seq, 'pep_seq': base_pep, 'label': 1,
                             'pep_id': orig, 'pdb': 'ablation'}])
        lb = DataLoader(PairDataset(dfb, tok, cfg.data.max_len, 'binary'),
                        batch_size=8, shuffle=False, collate_fn=collate_fn)
        scb, _ = predict(model, lb, device)
        scores['lasso'] = float(scb[0])
        for _, r in grp.iterrows():
            dfx = pd.DataFrame([{'prot_seq': prot_seq, 'pep_seq': r['seq'], 'label': 1, 'pep_id': r['pep_id'], 'pdb': 'ablation'}])
            loader = DataLoader(PairDataset(dfx, tok, cfg.data.max_len, 'binary'),
                                batch_size=8, shuffle=False, collate_fn=collate_fn)
            sc, _ = predict(model, loader, device)
            scores[r['kind']] = float(sc[0])
        rows.append(dict(orig=orig, **scores))
    out = pd.DataFrame(rows)
    out_path = os.path.join(cfg.paths.output_dir, 'ablation_scores.csv')
    out.to_csv(out_path, index=False)
    # summary stats + paired Wilcoxon (acceptance: >=80% drop, p<0.05)
    if {'ring-open', 'distal-control'}.issubset(out.columns):
        lasso_s = out.get('lasso', np.nan)
        d_open = (out['ring-open'] - lasso_s).dropna()
        d_ctrl = (out['distal-control'] - lasso_s).dropna()
        n_drop = (d_open < 0).sum()
        ratio = n_drop / max(len(d_open), 1)
        print(f'ring-open vs lasso: {n_drop}/{len(d_open)} pairs with lower open-score '
              f'({ratio:.0%}; acceptance >=80%), mean delta {d_open.mean():.3f}')
        print(f'distal-control vs lasso: {(d_ctrl < 0).sum()}/{len(d_ctrl)} lower '
              f'(expect ~0), mean delta {d_ctrl.mean():.3f}')
        try:
            from scipy.stats import wilcoxon
            if len(d_open) >= 6:
                w, p = wilcoxon(d_open, alternative='less')
                print(f'paired Wilcoxon ring-open: W={w:.1f} p={p:.4f} -> '
                      f'{"PASS" if p < 0.05 else "FAIL"} (p<0.05)')
            if len(d_ctrl) >= 6:
                w2, p2 = wilcoxon(d_ctrl)
                print(f'paired Wilcoxon distal-control: W={w2:.1f} p={p2:.4f} '
                      f'(expect NOT significant)')
        except ImportError:
            print('scipy not installed; skipping Wilcoxon test')
    return out


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', default='config.yaml')
    parser.add_argument('--ckpt', required=True)
    parser.add_argument('--device', default='auto')
    args = parser.parse_args()
    cfg = Config.from_yaml(args.config)
    cfg._base_dir = os.path.dirname(os.path.abspath(args.config))
    device = (torch.device('cuda' if torch.cuda.is_available() else 'cpu')
              if args.device == 'auto' else torch.device(args.device))
    tok = AutoTokenizer.from_pretrained(cfg.model.esm_model)
    model = build_loaded(cfg, args.ckpt, device)
    df = load_pairs(cfg.data.pairs_csv)
    loader = DataLoader(PairDataset(df, tok, cfg.data.max_len, cfg.train.label_type),
                        batch_size=cfg.train.batch_size, shuffle=False, collate_fn=collate_fn)
    sc, keys = predict(model, loader, device)
    auc, ap = evaluate_scores(sc, df['label'].values)
    print(f'held-out: AUC {auc:.3f} AP {ap:.3f}')
    if cfg.eval.ablation_csv:
        ring_ablation(cfg, model, tok, device, df, cfg.eval.ablation_csv)


def build_loaded(cfg, ckpt, device):
    model = InteractionPredictor(
        esm_model_name=cfg.model.esm_model,
        lassoesm_model_name=cfg.model.lassoesm_model,
        h_dim=cfg.model.h_dim, n_heads=cfg.model.n_heads, k=cfg.model.k,
        dropout=cfg.model.dropout, ban_dropout=cfg.model.ban_dropout,
        freeze_esm=cfg.model.freeze_esm,
        lassoesm_unfreeze=cfg.model.lassoesm_unfreeze,
        dual_head=cfg.model.dual_head).to(device)
    return load_ckpt(model, ckpt, device)


if __name__ == '__main__':
    main()
