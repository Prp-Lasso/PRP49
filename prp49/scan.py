"""Scan: score peptide x target matrix and export CSV + heatmap data."""
import os
import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader
from transformers import AutoTokenizer

from .config import Config
from .data import PairDataset, collate_fn
from .train import predict
from .eval import build_loaded


def load_fasta(path):
    seqs, cur, buf = [], None, []
    for line in open(path, encoding='utf-8'):
        line = line.strip()
        if line.startswith('>'):
            if cur:
                seqs.append((cur, ''.join(buf)))
            cur = line[1:]
            buf = []
        elif line:
            buf.append(line)
    if cur:
        seqs.append((cur, ''.join(buf)))
    return seqs


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', default='config.yaml')
    parser.add_argument('--ckpt', required=True)
    parser.add_argument('--peptides', required=True)   # csv: id,seq
    parser.add_argument('--targets', required=True)    # fasta: target proteins
    parser.add_argument('--device', default='auto')
    args = parser.parse_args()
    cfg = Config.from_yaml(args.config)
    cfg._base_dir = os.path.dirname(os.path.abspath(args.config))
    device = (torch.device('cuda' if torch.cuda.is_available() else 'cpu')
              if args.device == 'auto' else torch.device(args.device))
    tok = AutoTokenizer.from_pretrained(cfg.model.esm_model)
    model = build_loaded(cfg, args.ckpt, device)

    peps = pd.read_csv(args.peptides)
    targets = load_fasta(args.targets)
    mat = np.zeros((len(peps), len(targets)))
    for j, (tid, tseq) in enumerate(targets):
        rows = [{'prot_seq': tseq, 'pep_seq': s, 'label': 1,
                 'pep_id': f'{i}', 'pdb': 'scan'} for i, s in enumerate(peps['seq'])]
        dfx = pd.DataFrame(rows)
        loader = DataLoader(PairDataset(dfx, tok, cfg.data.max_len, 'binary'),
                            batch_size=8, shuffle=False, collate_fn=collate_fn)
        sc, _ = predict(model, loader, device)
        mat[:, j] = sc
    out = pd.DataFrame(mat, index=peps['id'], columns=[t[0] for t in targets])
    os.makedirs(cfg.scan.output and os.path.dirname(cfg.scan.output) or '.', exist_ok=True)
    out.to_csv(cfg.scan.output or 'scan_matrix.csv')
    print(out.round(3))


if __name__ == '__main__':
    main()
