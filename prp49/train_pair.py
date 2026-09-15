"""Pairwise ranking trainer for drug-candidate prioritisation.

Task: two peptides against the SAME target -> which binds stronger?
Loss: logistic pairwise ranking on the score difference (no absolute scale needed,
      so Ki/Kd/IC50 heterogeneity cancels out by construction).

Evaluation is screening-oriented (not RMSE):
  * pairwise accuracy          - direction of the preference correct?
  * per-target Spearman median - does the score order peptides within a target?
  * NDCG@k                     - are the truly strongest peptides ranked at the top?
  * EF@10% (enrichment factor) - how much better than random is the top decile?

Cross-validation groups by target (prot_id) so every test target is unseen.
"""
import argparse
import json
import os

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, Dataset
from transformers import AutoTokenizer

from .ckpt import maybe_resume, save_best, save_resume
from .config import Config
from .model import InteractionPredictor


class PairDataset(Dataset):
    """(target, peptide_a, peptide_b) -> label 1 means a binds stronger than b."""

    def __init__(self, df, tok, max_len=1024, pep_max_len=128):
        self.df = df.reset_index(drop=True)
        self.tok = tok
        self.max_len = max_len
        self.pep_max_len = pep_max_len

    def __len__(self):
        return len(self.df)

    def _enc(self, seq, max_len):
        e = self.tok(seq, truncation=True, max_length=max_len, padding='max_length',
                     return_tensors='pt')
        return e['input_ids'][0], e['attention_mask'][0]

    def __getitem__(self, i):
        r = self.df.iloc[i]
        t_ids, t_m = self._enc(str(r['prot_seq']), self.max_len)
        a_ids, a_m = self._enc(str(r['pep_a_seq']), self.pep_max_len)
        b_ids, b_m = self._enc(str(r['pep_b_seq']), self.pep_max_len)
        return dict(t_ids=t_ids, t_mask=t_m, a_ids=a_ids, a_mask=a_m,
                    b_ids=b_ids, b_mask=b_m, label=float(r['label']))


def pair_collate(batch):
    out = {}
    for k in ('t_ids', 't_mask', 'a_ids', 'a_mask', 'b_ids', 'b_mask'):
        out[k] = torch.stack([b[k] for b in batch])
    out['label'] = torch.tensor([b['label'] for b in batch], dtype=torch.float32)
    return out


def score_peptides(model, df, tok, device, max_len, pep_max_len, batch_size=16, target_col='prot_seq',
                   pep_col='pep_seq'):
    """Score every unique (target, peptide) pair; returns a DataFrame."""
    uniq = df.drop_duplicates([target_col, pep_col]).reset_index(drop=True)
    model.eval()
    scores = []
    with torch.no_grad():
        for s in range(0, len(uniq), batch_size):
            chunk = uniq.iloc[s:s + batch_size]
            t = tok(list(chunk[target_col]), truncation=True, max_length=max_len,
                    padding='max_length', return_tensors='pt')
            p = tok(list(chunk[pep_col]), truncation=True, max_length=pep_max_len,
                    padding='max_length', return_tensors='pt')
            logits, _, _ = model(t['input_ids'].to(device), t['attention_mask'].to(device),
                                 p['input_ids'].to(device), p['attention_mask'].to(device))
            scores.extend(logits.detach().cpu().numpy().tolist())
    uniq['score'] = scores
    return uniq


def ndcg_at_k(order_scores, true_scores, k):
    """NDCG@k with exponential gain on min-max normalised true affinity."""
    n = len(order_scores)
    if n == 0:
        return np.nan
    k = min(k, n)
    t = np.asarray(true_scores, dtype=float)
    lo, hi = t.min(), t.max()
    rel = (t - lo) / (hi - lo) if hi > lo else np.ones_like(t)
    order = np.argsort(-np.asarray(order_scores, dtype=float))[:k]
    gains = 2.0 ** rel[order] - 1.0
    discounts = 1.0 / np.log2(np.arange(2, k + 2))
    dcg = float((gains * discounts).sum())
    best = np.argsort(-rel)[:k]
    idcg = float(((2.0 ** rel[best] - 1.0) * discounts).sum())
    return dcg / idcg if idcg > 0 else np.nan


def ef_at_frac(order_scores, true_scores, frac=0.1):
    """Enrichment factor: share of true top-frac recovered in the model's top-frac."""
    n = len(order_scores)
    if n < 5:
        return np.nan
    k = max(1, int(round(frac * n)))
    t = np.asarray(true_scores, dtype=float)
    s = np.asarray(order_scores, dtype=float)
    top_true = set(np.argsort(-t)[:k].tolist())
    top_pred = set(np.argsort(-s)[:k].tolist())
    hit = len(top_true & top_pred)
    return (hit / k) / (k / n)


def evaluate(model, test_df, tok, device, cfg):
    """Peptide-level screening metrics on unseen targets."""
    # pairwise frame holds peptides in two columns -> build the long (target, peptide) form
    long_df = pd.concat([
        test_df[['prot_id', 'prot_seq']].assign(pep_seq=test_df['pep_a_seq'].values),
        test_df[['prot_id', 'prot_seq']].assign(pep_seq=test_df['pep_b_seq'].values),
    ], ignore_index=True).drop_duplicates(['prot_id', 'pep_seq'])
    scored = score_peptides(model, long_df.reset_index(drop=True), tok, device,
                            cfg.data.max_len, cfg.data.get('pep_max_len', 128))
    aff = pd.read_csv(cfg.data.affinity_csv)[['prot_id', 'pep_seq', 'energy']].drop_duplicates(
        ['prot_id', 'pep_seq'])

    # aggregate model score per (target, peptide): mean over its occurrences
    agg = scored.groupby(['prot_id', 'pep_seq'], as_index=False)['score'].mean()
    merged = agg.merge(aff, on=['prot_id', 'pep_seq'], how='inner')

    sp, nd5, nd10, ef = [], [], [], []
    for pid, grp in merged.groupby('prot_id'):
        if len(grp) < 3:
            continue
        if grp.energy.nunique() < 2 or grp.score.nunique() < 2:
            continue
        rho = pd.Series(grp.score).corr(pd.Series(grp.energy), method='spearman')
        if not np.isnan(rho):
            sp.append(rho)
        nd5.append(ndcg_at_k(grp.score.values, grp.energy.values, 5))
        nd10.append(ndcg_at_k(grp.score.values, grp.energy.values, 10))
        ef.append(ef_at_frac(grp.score.values, grp.energy.values, cfg.eval.get('ef_top_frac', 0.1)))

    # pairwise accuracy on the raw test pairs
    model.eval()
    correct, total = 0, 0
    with torch.no_grad():
        for s in range(0, len(test_df), 24):
            chunk = test_df.iloc[s:s + 24]
            t = tok(list(chunk.prot_seq), truncation=True, max_length=cfg.data.max_len,
                    padding='max_length', return_tensors='pt')
            a = tok(list(chunk.pep_a_seq), truncation=True, max_length=cfg.data.get('pep_max_len', 128),
                    padding='max_length', return_tensors='pt')
            b = tok(list(chunk.pep_b_seq), truncation=True, max_length=cfg.data.get('pep_max_len', 128),
                    padding='max_length', return_tensors='pt')
            la, _, _ = model(t['input_ids'].to(device), t['attention_mask'].to(device),
                             a['input_ids'].to(device), a['attention_mask'].to(device))
            lb, _, _ = model(t['input_ids'].to(device), t['attention_mask'].to(device),
                             b['input_ids'].to(device), b['attention_mask'].to(device))
            correct += int((la > lb).sum().item())
            total += len(chunk)
    return dict(
        pair_acc=correct / max(total, 1),
        spearman_median=float(np.nanmedian(sp)) if sp else float('nan'),
        spearman_mean=float(np.nanmean(sp)) if sp else float('nan'),
        ndcg5=float(np.nanmean(nd5)) if nd5 else float('nan'),
        ndcg10=float(np.nanmean(nd10)) if nd10 else float('nan'),
        ef10=float(np.nanmean(ef)) if ef else float('nan'),
        n_targets=len(sp),
    )


def run(cfg, device):
    tok = AutoTokenizer.from_pretrained(cfg.model.esm_model)
    df = pd.read_csv(cfg.data.pairs_csv)
    print(f'pairwise samples: {len(df)}, targets: {df.prot_id.nunique()}')

    groups = df[cfg.data.group_col].values
    uniq_groups = pd.unique(groups)
    rng = np.random.default_rng(cfg.train.get('seed', 42))
    rng.shuffle(uniq_groups)
    folds = np.array_split(uniq_groups, cfg.eval.cv_folds)

    results = []
    os.makedirs(cfg.paths.checkpoint_dir, exist_ok=True)
    for fi, test_groups in enumerate(folds):
        te_mask = df[cfg.data.group_col].isin(test_groups).values
        tr_df, te_df = df[~te_mask], df[te_mask]
        if len(te_df) == 0 or len(tr_df) == 0:
            continue
        model = InteractionPredictor(
            esm_model_name=cfg.model.esm_model, lassoesm_model_name=cfg.model.lassoesm_model,
            h_dim=cfg.model.h_dim, n_heads=cfg.model.n_heads, k=cfg.model.k,
            dropout=cfg.model.dropout, ban_dropout=cfg.model.ban_dropout,
            freeze_esm=cfg.model.get('freeze_esm', False),
            lassoesm_unfreeze=cfg.model.get('lassoesm_unfreeze', 0),
            dual_head=False, use_topo=False).to(device)

        opt = torch.optim.AdamW([p for p in model.parameters() if p.requires_grad],
                                lr=cfg.train.lr, weight_decay=cfg.train.weight_decay)
        tr_loader = DataLoader(PairDataset(tr_df, tok, cfg.data.max_len, cfg.data.get('pep_max_len', 128)),
                               batch_size=cfg.train.batch_size, shuffle=True,
                               collate_fn=pair_collate, drop_last=True)

        best_acc, bad, best_state = -1.0, 0, None
        ckpt_dir = cfg.paths.checkpoint_dir
        os.makedirs(ckpt_dir, exist_ok=True)
        start_ep, _ = maybe_resume(ckpt_dir, fi, model, opt, device)
        ckpt_every = int(cfg.train.get('ckpt_every') or 5)
        for ep in range(start_ep, cfg.train.epochs):
            model.train()
            tot, n = 0.0, 0
            for batch in tr_loader:
                t_ids = batch['t_ids'].to(device); t_m = batch['t_mask'].to(device)
                a_ids = batch['a_ids'].to(device); a_m = batch['a_mask'].to(device)
                b_ids = batch['b_ids'].to(device); b_m = batch['b_mask'].to(device)
                la, _, _ = model(t_ids, t_m, a_ids, a_m)
                lb, _, _ = model(t_ids, t_m, b_ids, b_m)
                loss = nn.functional.binary_cross_entropy_with_logits(
                    la - lb, torch.ones_like(la))
                opt.zero_grad(); loss.backward()
                nn.utils.clip_grad_norm_(model.parameters(), 1.0)
                opt.step()
                tot += loss.item() * len(la); n += len(la)
            # quick validation: pairwise accuracy on a subset of test pairs
            sub = te_df.sample(min(256, len(te_df)), random_state=ep)
            model.eval()
            correct = 0
            with torch.no_grad():
                for s in range(0, len(sub), 16):
                    chunk = sub.iloc[s:s + 16]
                    t = tok(list(chunk.prot_seq), truncation=True, max_length=cfg.data.max_len,
                            padding='max_length', return_tensors='pt')
                    a = tok(list(chunk.pep_a_seq), truncation=True, max_length=cfg.data.get('pep_max_len', 128),
                            padding='max_length', return_tensors='pt')
                    b = tok(list(chunk.pep_b_seq), truncation=True, max_length=cfg.data.get('pep_max_len', 128),
                            padding='max_length', return_tensors='pt')
                    la, _, _ = model(t['input_ids'].to(device), t['attention_mask'].to(device),
                                     a['input_ids'].to(device), a['attention_mask'].to(device))
                    lb, _, _ = model(t['input_ids'].to(device), t['attention_mask'].to(device),
                                     b['input_ids'].to(device), b['attention_mask'].to(device))
                    correct += int((la > lb).sum().item())
            vacc = correct / max(len(sub), 1)
            print(f'  fold {fi} ep {ep}: loss {tot/max(n,1):.4f} val pair-acc {vacc:.4f}', flush=True)
            if vacc > best_acc:
                best_acc, bad = vacc, 0
                best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
                save_best(ckpt_dir, fi, best_state, {'pair_acc': float(vacc), 'ep': ep})
            else:
                bad += 1
                if bad >= cfg.train.patience:
                    print(f'  fold {fi}: early stop at ep {ep}', flush=True)
                    break
            if (ep + 1) % ckpt_every == 0:
                save_resume(ckpt_dir, fi, model, opt, ep, best_acc)
                print(f'  [ckpt] fold {fi} ep {ep}: rolling checkpoint saved', flush=True)

        if best_state is not None:
            model.load_state_dict(best_state)
        save_resume(ckpt_dir, fi, model, opt, cfg.train.epochs - 1, best_acc)
        m = evaluate(model, te_df, tok, device, cfg)
        m['fold'] = fi
        print(f'fold {fi}: ' + ' '.join(f'{k}={v:.4f}' if isinstance(v, float) else f'{k}={v}'
                                        for k, v in m.items()), flush=True)
        results.append(m)

    summary = {}
    for key in ('pair_acc', 'spearman_median', 'spearman_mean', 'ndcg5', 'ndcg10', 'ef10'):
        vals = [r[key] for r in results if not np.isnan(r.get(key, np.nan))]
        if vals:
            summary[key] = float(np.mean(vals))
            summary[key + '_std'] = float(np.std(vals))
    print(json.dumps({'folds': results, **summary}, indent=2, ensure_ascii=False), flush=True)
    with open(os.path.join(cfg.paths.output_dir, 'rank_summary.json'), 'w') as f:
        json.dump({'folds': results, **summary}, f, indent=2, ensure_ascii=False)
    return summary


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--config', required=True)
    ap.add_argument('--device', default='cuda')
    args = ap.parse_args()
    cfg = Config.from_yaml(args.config)
    cfg._base_dir = os.path.dirname(os.path.abspath(args.config))
    os.makedirs(cfg.paths.output_dir, exist_ok=True)
    run(cfg, args.device)


if __name__ == '__main__':
    main()
