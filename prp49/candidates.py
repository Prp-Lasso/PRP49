"""Candidate prioritisation for druggability screening.

Combines every signal this project produces into one ranked table:

  binding_prob  - interaction classifier  (is this peptide-target pair a binder?)
  rank_score    - pairwise ranking model  (within-target preference among binders)
  lasso_prob    - Lasso peptide classifier (topology plausibility / lasso-likeness)
  dock_score    - AutoDock Vina score, optional (structural feasibility)

Composite = weighted z-score fusion (weights configurable).  The output is the
deliverable a screening funnel actually consumes: a shortlist, not a Kd number.

Usage:
  python -m prp49.candidates --config config_rank.yaml \
      --peptides ../mvp_cpu/scan_peptides_ext.csv --targets ../mvp_cpu/scan_targets.fasta \
      --cls_ckpt runs_improved/checkpoints/fold0_best.pt \
      --rank_ckpt runs_rank/checkpoints/fold0_best.pt \
      --lasso_ckpt ../LassoPeptideClassifier/checkpoints/best_model.pt \
      --out ../results/candidates.csv
"""
import argparse
import os

import numpy as np
import pandas as pd
import torch
from transformers import AutoTokenizer

from .config import Config
from .model import InteractionPredictor


def load_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path, errors='ignore'):
        line = line.strip()
        if not line or line[0] in ';#!':
            continue
        if line.startswith('>'):
            if name:
                seqs[name] = ''.join(buf)
            name, buf = line[1:].split()[0].split('|')[0], []   # "EDNRB|P24530|human" -> EDNRB
        else:
            buf.append(line)
    if name:
        seqs[name] = ''.join(buf)
    return seqs


def load_peptides(path):
    df = pd.read_csv(path)
    col = next(c for c in df.columns if c.lower() in ('pep_seq', 'sequence', 'seq'))
    idc = next((c for c in df.columns if c.lower() in ('pep_id', 'id', 'name')), None)
    ids = df[idc].astype(str).tolist() if idc else [f'pep{i}' for i in range(len(df))]
    return list(zip(ids, df[col].astype(str).tolist()))


@torch.no_grad()
def score_pairs(model, tok, device, pairs, max_len=1024, pep_max_len=128, batch=16):
    """pairs: list of (target_seq, peptide_seq) -> logits."""
    model.eval()
    out = []
    for s in range(0, len(pairs), batch):
        chunk = pairs[s:s + batch]
        t = tok([c[0] for c in chunk], truncation=True, max_length=max_len,
                padding='max_length', return_tensors='pt')
        p = tok([c[1] for c in chunk], truncation=True, max_length=pep_max_len,
                padding='max_length', return_tensors='pt')
        logits, _, _ = model(t['input_ids'].to(device), t['attention_mask'].to(device),
                             p['input_ids'].to(device), p['attention_mask'].to(device))
        out.extend(logits.detach().cpu().numpy().tolist())
    return out


def zscore(x):
    x = np.asarray(x, dtype=float)
    sd = x.std()
    return (x - x.mean()) / sd if sd > 1e-9 else np.zeros_like(x)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--config', required=True)
    ap.add_argument('--peptides', required=True, help='CSV with a sequence column')
    ap.add_argument('--targets', required=True, help='FASTA of target sequences')
    ap.add_argument('--cls_ckpt', required=True, help='interaction classifier checkpoint')
    ap.add_argument('--rank_ckpt', default=None, help='pairwise ranking checkpoint (optional)')
    ap.add_argument('--lasso_ckpt', default=None, help='Lasso peptide classifier (optional)')
    ap.add_argument('--dock_csv', default=None, help='CSV with peptide,target,vina score columns')
    # Default weights follow a validation finding (2026-09-13): the pairwise
    # ranking head was trained on ordinary peptides (BindingDB/ChEMBL) and
    # transfers poorly to lasso peptides - adding it moved literature-supported
    # pairs from mean rank 23 to 71. Structural docking evidence, by contrast,
    # moved them from 59 (binding only) to 23. So: binding + docking by default,
    # rank/lasso available but off unless the caller has evidence for them.
    ap.add_argument('--w_binding', type=float, default=0.60)
    ap.add_argument('--w_rank', type=float, default=0.0)
    ap.add_argument('--w_lasso', type=float, default=0.0)
    ap.add_argument('--w_dock', type=float, default=0.40)
    ap.add_argument('--out', default='candidates.csv')
    args = ap.parse_args()

    cfg = Config.from_yaml(args.config)
    cfg._base_dir = os.path.dirname(os.path.abspath(args.config))
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    tok = AutoTokenizer.from_pretrained(cfg.model.esm_model)

    peps = load_peptides(args.peptides)
    tgts = load_fasta(args.targets)
    print(f'candidates: {len(peps)} peptides x {len(tgts)} targets = {len(peps)*len(tgts)} pairs')

    rows = []
    for pid, pseq in peps:
        for tid, tseq in tgts.items():
            rows.append(dict(peptide_id=pid, peptide=pseq,
                             target_id=tid, target=tseq))
    df = pd.DataFrame(rows)

    # --- binding probability
    model = InteractionPredictor(
        esm_model_name=cfg.model.esm_model, lassoesm_model_name=cfg.model.lassoesm_model,
        h_dim=cfg.model.h_dim, n_heads=cfg.model.n_heads, k=cfg.model.k,
        dropout=cfg.model.dropout, ban_dropout=cfg.model.ban_dropout,
        freeze_esm=False, lassoesm_unfreeze=0, dual_head=False, use_topo=False).to(device)
    sd = torch.load(args.cls_ckpt, map_location=device)
    model.load_state_dict(sd.get('state_dict', sd), strict=False)
    df['binding_logit'] = score_pairs(model, tok, device,
                                      list(zip(df.target, df.peptide)),
                                      cfg.data.max_len, cfg.data.get('pep_max_len', 128))
    df['binding_prob'] = 1 / (1 + np.exp(-df.binding_logit))

    # --- pairwise ranking score (optional)
    if args.rank_ckpt and os.path.exists(args.rank_ckpt):
        rmodel = InteractionPredictor(
            esm_model_name=cfg.model.esm_model, lassoesm_model_name=cfg.model.lassoesm_model,
            h_dim=cfg.model.h_dim, n_heads=cfg.model.n_heads, k=cfg.model.k,
            dropout=cfg.model.dropout, ban_dropout=cfg.model.ban_dropout,
            freeze_esm=False, lassoesm_unfreeze=0, dual_head=False, use_topo=False).to(device)
        rsd = torch.load(args.rank_ckpt, map_location=device)
        rmodel.load_state_dict(rsd.get('state_dict', rsd), strict=False)
        df['rank_score'] = score_pairs(rmodel, tok, device, list(zip(df.target, df.peptide)),
                                       cfg.data.max_len, cfg.data.get('pep_max_len', 128))
    else:
        df['rank_score'] = np.nan

    # --- lasso probability (optional; peptide-only model)
    if args.lasso_ckpt and os.path.exists(args.lasso_ckpt):
        try:
            import sys
            lasso_dir = os.path.dirname(os.path.dirname(os.path.abspath(args.lasso_ckpt)))
            if lasso_dir not in sys.path:
                sys.path.insert(0, lasso_dir)
            from model import LassoPeptideClassifier  # type: ignore
            from utils import load_classifier_from_checkpoint  # type: ignore
            lmodel = load_classifier_from_checkpoint(args.lasso_ckpt, device=device)
            pep_list = sorted({p for _, p in peps})
            emb_tok = AutoTokenizer.from_pretrained(cfg.model.esm_model)
            from transformers import AutoModel
            esm = AutoModel.from_pretrained(cfg.model.esm_model).to(device).eval()
            probs = {}
            with torch.no_grad():
                for p in pep_list:
                    e = emb_tok(p, return_tensors='pt', truncation=True, max_length=128).to(device)
                    h = esm(**e).last_hidden_state[:, 1:-1, :]
                    probs[p] = float(torch.sigmoid(lmodel(h)).mean())
            df['lasso_prob'] = df.peptide.map(probs)
        except Exception as exc:      # keep the pipeline usable without the classifier
            print(f'[warn] lasso classifier unavailable: {exc}')
            df['lasso_prob'] = np.nan
    else:
        df['lasso_prob'] = np.nan

    # --- docking score (optional; lower is better -> invert for the composite)
    if args.dock_csv and os.path.exists(args.dock_csv):
        dk = pd.read_csv(args.dock_csv)
        pcol = next((c for c in dk.columns if c.lower() in ('peptide', 'pep_seq', 'pep_name')), None)
        tcol = next((c for c in dk.columns if c.lower() in ('target', 'prot_id', 'rec_seq')), None)
        scol = next((c for c in dk.columns if 'score' in c.lower()), None)
        if pcol and tcol and scol:
            dk = dk.rename(columns={pcol: 'peptide_id', tcol: 'target_id', scol: 'dock_score'})
            # match on IDs (the peptide/target sequence columns are not unique keys here)
            df = df.merge(dk[['peptide_id', 'target_id', 'dock_score']],
                          on=['peptide_id', 'target_id'], how='left')
        else:
            df['dock_score'] = np.nan
    else:
        df['dock_score'] = np.nan

    # --- composite score.
    # rank_score is a *within-target* preference score, so normalise it per target;
    # the other signals are pair-level and are normalised globally.
    comp = args.w_binding * zscore(df.binding_prob.fillna(df.binding_prob.mean()))
    if df.rank_score.notna().any():
        per_target = df.groupby('target_id')['rank_score'].transform(
            lambda s: (s - s.mean()) / (s.std(ddof=0) + 1e-9) if s.notna().sum() > 1 else 0.0)
        comp = comp + args.w_rank * per_target.fillna(0.0)
    if df.lasso_prob.notna().any():
        comp = comp + args.w_lasso * zscore(df.lasso_prob.fillna(df.lasso_prob.mean()))
    if df.dock_score.notna().any():
        comp = comp + args.w_dock * zscore(-df.dock_score.fillna(df.dock_score.mean()))
    df['composite'] = comp
    df['rank'] = df.composite.rank(ascending=False).astype(int)

    cols = ['rank', 'peptide_id', 'target_id', 'composite', 'binding_prob', 'binding_logit',
            'rank_score', 'lasso_prob', 'dock_score', 'peptide', 'target']
    df = df[[c for c in cols if c in df.columns]].sort_values('rank')
    df.round(4).to_csv(args.out, index=False, lineterminator='\n')

    # coverage / calibration diagnostics (absolute probabilities matter for screening)
    print(f'wrote {args.out}: {len(df)} pairs')
    print('score coverage: ' + ', '.join(f'{c}={df[c].notna().mean()*100:.0f}%'
                                         for c in ('binding_prob', 'rank_score', 'lasso_prob', 'dock_score')))
    print(f'binding_prob distribution: min {df.binding_prob.min():.3f} · median {df.binding_prob.median():.3f} '
          f'· max {df.binding_prob.max():.3f} · frac>0.9 {100*(df.binding_prob > 0.9).mean():.0f}%')
    if df.lasso_prob.notna().any():
        print(f'lasso_prob distribution: min {df.lasso_prob.min():.3f} · median {df.lasso_prob.median():.3f} '
              f'· max {df.lasso_prob.max():.3f}')
    if df.dock_score.notna().any():
        print(f'dock_score distribution: min {df.dock_score.min():.1f} · median {df.dock_score.median():.1f} '
              f'· max {df.dock_score.max():.1f}')
    print(df[['rank', 'peptide_id', 'target_id', 'composite', 'binding_prob',
              'rank_score', 'lasso_prob', 'dock_score']].head(12).to_string(index=False))


if __name__ == '__main__':
    main()
