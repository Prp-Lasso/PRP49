"""Route C: family-aware contrastive pretraining of the peptide encoder.

Rationale: the supervised head sees only 352 positive pairs. The LassoPred library
holds 4,749 cores in 227 similarity families; pulling same-family members together
and pushing different families apart teaches family-level structure with no labels
at all, which is exactly the kind of prior the small supervised set cannot supply.

Design
  * encoder   : LassoESM (top layers unfrozen, low LR)
  * head      : 2-layer MLP -> 128-d projection
  * loss      : NT-Xent / InfoNCE, temperature 0.07, in-batch negatives
  * sampling  : one positive per family per step, so the 1,109-member family does
                not drown out the small ones
  * output    : HuggingFace-format encoder directory, loadable straight into
                InteractionPredictor(lassoesm_model_name=...)
"""
import argparse
import json
import math
import os
import random
import time

import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, Dataset
from transformers import AutoModel, AutoTokenizer


class FamilyPairs(Dataset):
    """Yields (anchor, positive) drawn from the same family.

    Families are sampled uniformly (one item per family per epoch slot) rather than
    per sequence, otherwise the largest family dominates every batch.
    """

    def __init__(self, df, steps_per_epoch, seed=0):
        self.members = {f: g.sequence.tolist() for f, g in df.groupby('family_id') if len(g) > 1}
        self.families = list(self.members.keys())
        self.steps = steps_per_epoch
        self.rng = random.Random(seed)

    def __len__(self):
        return self.steps

    def __getitem__(self, _):
        f = self.rng.choice(self.families)
        a, b = self.rng.sample(self.members[f], 2)
        return a, b


def nt_xent(z1, z2, temperature):
    """Standard SimCLR/InfoNCE over a batch of paired embeddings."""
    n = z1.size(0)
    z = torch.cat([z1, z2], dim=0)
    z = F.normalize(z, dim=1)
    sim = z @ z.t() / temperature
    sim.fill_diagonal_(-1e9)
    # positive of row i (i<n) is i+n and vice versa
    targets = torch.cat([torch.arange(n, 2 * n), torch.arange(0, n)]).to(z.device)
    return F.cross_entropy(sim, targets)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data', required=True)
    ap.add_argument('--encoder', required=True, help='starting LassoESM directory')
    ap.add_argument('--out', required=True, help='output encoder directory (HF format)')
    ap.add_argument('--epochs', type=int, default=15)
    ap.add_argument('--batch', type=int, default=128)
    ap.add_argument('--lr', type=float, default=1e-5)
    ap.add_argument('--temperature', type=float, default=0.07)
    ap.add_argument('--unfreeze_layers', type=int, default=4)
    ap.add_argument('--max_len', type=int, default=64)
    ap.add_argument('--device', default='cuda')
    args = ap.parse_args()

    torch.manual_seed(0)
    df = pd.read_csv(args.data)
    df = df.dropna(subset=['sequence', 'family_id'])
    usable = df.groupby('family_id').filter(lambda g: len(g) > 1)
    print(f'sequences {len(df)} | families {df.family_id.nunique()} | '
          f'usable (families with >1 member): {len(usable)} '
          f'in {usable.family_id.nunique()} families', flush=True)

    tok = AutoTokenizer.from_pretrained(args.encoder)
    enc = AutoModel.from_pretrained(args.encoder)
    for p in enc.parameters():
        p.requires_grad = False
    layers = getattr(getattr(enc, 'encoder', None), 'layer', None)
    if layers is not None and args.unfreeze_layers > 0:
        for layer in layers[-args.unfreeze_layers:]:
            for p in layer.parameters():
                p.requires_grad = True
    n_train = sum(p.numel() for p in enc.parameters() if p.requires_grad)
    print(f'unfroze top {args.unfreeze_layers} layers -> {n_train/1e6:.1f}M trainable params',
          flush=True)

    class Proj(nn.Module):
        def __init__(self, d_in):
            super().__init__()
            self.net = nn.Sequential(nn.Linear(d_in, 256), nn.ReLU(),
                                     nn.Linear(256, 128))

        def forward(self, x):
            return self.net(x)

    proj = Proj(enc.config.hidden_size).to(args.device)
    enc = enc.to(args.device)
    params = [p for p in enc.parameters() if p.requires_grad] + list(proj.parameters())
    opt = torch.optim.AdamW(params, lr=args.lr, weight_decay=0.01)
    steps_per_epoch = max(20, len(usable) // 4)          # one slot per ~4 sequences
    print(f'steps/epoch = {steps_per_epoch} (batch {args.batch})', flush=True)
    ds = FamilyPairs(usable, steps_per_epoch)
    dl = DataLoader(ds, batch_size=args.batch, shuffle=False, drop_last=True)

    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=args.epochs)
    hist = []
    for ep in range(args.epochs):
        t0 = time.time()
        enc.train(); proj.train()
        tot, nb, pos_sim, neg_sim = 0.0, 0, 0.0, 0.0
        for a, b in dl:
            ta = tok(list(a), truncation=True, max_length=args.max_len,
                     padding='max_length', return_tensors='pt').to(args.device)
            tb = tok(list(b), truncation=True, max_length=args.max_len,
                     padding='max_length', return_tensors='pt').to(args.device)
            oa = enc(**ta).last_hidden_state[:, 0]
            ob = enc(**tb).last_hidden_state[:, 0]
            za, zb = proj(oa), proj(ob)
            loss = nt_xent(za, zb, args.temperature)
            opt.zero_grad(); loss.backward()
            torch.nn.utils.clip_grad_norm_(params, 1.0)
            opt.step()
            with torch.no_grad():
                zn = F.normalize(torch.cat([za, zb]), dim=1)
                s = zn @ zn.t()
                n = za.size(0)
                diag = torch.cat([s[range(n), range(n, 2 * n)], s[range(n, 2 * n), range(n)]])
                pos_sim += float(diag.mean())
                mask = ~torch.eye(2 * n, dtype=torch.bool, device=s.device)
                neg_sim += float(s[mask].mean())
            tot += float(loss); nb += 1
        sched.step()
        hist.append(dict(epoch=ep, loss=tot / max(nb, 1), pos_sim=pos_sim / max(nb, 1),
                         neg_sim=neg_sim / max(nb, 1), seconds=round(time.time() - t0)))
        print(f'  ep {ep}: loss {tot/max(nb,1):.4f} | pos-sim {pos_sim/max(nb,1):.3f} '
              f'| neg-sim {neg_sim/max(nb,1):.3f} | {time.time()-t0:.0f}s', flush=True)

    os.makedirs(args.out, exist_ok=True)
    enc.save_pretrained(args.out)
    tok.save_pretrained(args.out)
    json.dump(dict(args=vars(args), history=hist,
                   trainable_params=int(n_train),
                   separation_gain=hist[-1]['pos_sim'] - hist[-1]['neg_sim']),
              open(os.path.join(args.out, 'contrastive_report.json'), 'w'), indent=2)
    print(f'\nsaved encoder -> {args.out}', flush=True)
    print(f"pos-neg separation: {hist[0]['pos_sim']-hist[0]['neg_sim']:.3f} -> "
          f"{hist[-1]['pos_sim']-hist[-1]['neg_sim']:.3f} "
          f"(higher = better family structure)", flush=True)
    print('CONTRASTIVE_DONE')


if __name__ == '__main__':
    main()
