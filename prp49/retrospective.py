"""Retrospective validation (acceptance: >=70% known pairs rank in top 20% per target).

Scores held-out known interaction pairs with a trained checkpoint and ranks
each against a candidate pool for the same target (10 scan peptides by
default, plus any same-target negatives given via --neg_csv).
Output: retro_rank.csv (per-pair percentile) + printed pass/fail summary.
"""
import os
import argparse
import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader
from transformers import AutoTokenizer

from .config import Config
from .data import PairDataset, collate_fn
from .train import predict
from .eval import build_loaded


def score_rows(model, tok, cfg, device, rows):
    loader = DataLoader(PairDataset(rows, tok, cfg.data.max_len, 'binary'),
                        batch_size=8, shuffle=False, collate_fn=collate_fn)
    sc, _ = predict(model, loader, device)
    return sc


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', default='config.yaml')
    parser.add_argument('--ckpt', required=True)
    parser.add_argument('--retro', required=True, help='csv: pep_id,pep_seq,prot_seq')
    parser.add_argument('--peptides', default=None,
                        help='pool peptides csv (id,seq); default scan_peptides.csv next to retro')
    parser.add_argument('--neg_csv', default=None,
                        help='optional same-target negatives csv (pep_id,pep_seq,prot_seq,label)')
    parser.add_argument('--device', default='auto')
    parser.add_argument('--out', default='retro_rank.csv')
    args = parser.parse_args()
    cfg = Config.from_yaml(args.config)
    cfg._base_dir = os.path.dirname(os.path.abspath(args.config))
    device = (torch.device('cuda' if torch.cuda.is_available() else 'cpu')
              if args.device == 'auto' else torch.device(args.device))
    tok = AutoTokenizer.from_pretrained(cfg.model.esm_model)
    model = build_loaded(cfg, args.ckpt, device)

    retro = pd.read_csv(args.retro)
    pep_path = args.peptides or os.path.join(os.path.dirname(os.path.abspath(args.retro)),
                                             'scan_peptides.csv')
    pool_peps = pd.read_csv(pep_path) if os.path.exists(pep_path) else None
    neg = pd.read_csv(args.neg_csv) if args.neg_csv and os.path.exists(args.neg_csv) else None

    results = []
    for _, r in retro.iterrows():
        prot_seq = r['prot_seq']
        pool_rows = []
        if pool_peps is not None:
            for _, p in pool_peps.iterrows():
                pool_rows.append({'pep_id': p['id'], 'pep_seq': p['seq'],
                                  'prot_seq': prot_seq, 'label': 1})
        if neg is not None:
            same = neg[neg['prot_seq'] == prot_seq]
            for _, n in same.iterrows():
                pool_rows.append({'pep_id': n['pep_id'], 'pep_seq': n['pep_seq'],
                                  'prot_seq': prot_seq, 'label': float(n.get('label', 0))})
        # ensure the retro pair itself is in the pool (idempotent if duplicated)
        pool_rows.append({'pep_id': r['pep_id'], 'pep_seq': r['pep_seq'],
                          'prot_seq': prot_seq, 'label': 1})
        pool_df = pd.DataFrame(pool_rows).drop_duplicates(subset=['pep_seq', 'prot_seq'])
        sc = score_rows(model, tok, cfg, device, pool_df)
        pool_df = pool_df.copy()
        pool_df['score'] = sc
        retro_score = pool_df.loc[
            (pool_df['pep_seq'] == r['pep_seq']) & (pool_df['prot_seq'] == prot_seq),
            'score'].max()
        rank = (pool_df['score'] > retro_score).sum() + 1
        pct = rank / len(pool_df)
        results.append(dict(pep_id=r['pep_id'], prot_id=r.get('prot_id', ''),
                            score=float(retro_score), rank=int(rank),
                            pool_size=len(pool_df), percentile=float(pct)))

    out = pd.DataFrame(results)
    out.to_csv(args.out, index=False)
    n_top20 = int((out['percentile'] <= 0.20).sum())
    print(out.round(3).to_string(index=False))
    print(f'retrospective: {n_top20}/{len(out)} known pairs in top-20% '
          f'({n_top20 / max(len(out), 1):.0%}; acceptance >=70%)')
    print(f'saved {args.out}')


if __name__ == '__main__':
    main()
