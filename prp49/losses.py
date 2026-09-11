"""Losses: binary BCE (+PU weighting) + attention alignment loss (L_align)."""
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F


def binary_loss(logits, labels, pos_weight=1.0, sample_weight=None):
    """BCE with optional per-sample weights (PU learning hook)."""
    loss = F.binary_cross_entropy_with_logits(
        logits, labels,
        pos_weight=torch.tensor(pos_weight, dtype=logits.dtype, device=logits.device),
        reduction='none')
    if sample_weight is not None:
        loss = loss * sample_weight
    return loss.mean()


def ranking_loss(logits, labels, margin=0.5):
    """Pairwise margin ranking (I2): positives must outscore negatives.

    Aligns the training objective with the delivery metric (per-target ranking)
    instead of pure global classification.
    """
    pos = logits[labels > 0.5]
    neg = logits[labels <= 0.5]
    if pos.numel() == 0 or neg.numel() == 0:
        return None
    diff = pos.unsqueeze(1) - neg.unsqueeze(0)      # (P, N)
    return F.relu(margin - diff).mean()


def alignment_loss(att_maps, contact_mats, pep_lens=None, prot_lens=None, mode='kl'):
    """Attention alignment loss (思路2): pull BAN attention toward true contacts.

    att_maps: (B, heads, prot_len, pep_len)
    contact_mats: list of (prot_len_true, pep_len_true) binary matrices (or None
                  for pairs without structure annotation -> excluded)
    mode: 'kl' | 'mse'
    """
    losses = []
    count = 0
    for i, M in enumerate(contact_mats):
        if M is None:
            continue
        hp, wp = att_maps[i].shape[-2], att_maps[i].shape[-1]
        Mt = torch.tensor(M, dtype=att_maps.dtype, device=att_maps.device)
        if Mt.sum() < 1:
            continue
        # pad true matrix into attention (padded) space: batch may contain
        # pairs of different lengths, and tokenizer truncation keeps the head
        # of the sequence, so truncate the matrix to (hp, wp) accordingly
        Mp = torch.zeros(hp, wp, dtype=att_maps.dtype, device=att_maps.device)
        h, w = Mt.shape
        hh, ww = min(h, hp), min(w, wp)
        Mp[:hh, :ww] = Mt[:hh, :ww]
        Mt = Mp
        # mean attention over heads
        a = att_maps[i].mean(dim=0)                    # (prot_len, pep_len)
        a = F.log_softmax(a.view(-1), dim=0).view(hp, wp)
        target = Mt / Mt.sum()
        if mode == 'kl':
            l = F.kl_div(a, target, reduction='sum')
        else:
            l = F.mse_loss(F.softmax(a.view(-1), dim=0), target.view(-1), reduction='sum')
        losses.append(l)
        count += 1
    if count == 0:
        return None
    return torch.stack(losses).mean()
