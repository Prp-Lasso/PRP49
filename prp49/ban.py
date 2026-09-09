"""BAN: Bilinear Attention Network (low-rank) — v=protein, q=peptide.

Returns fused feature and multi-head attention maps for interpretability
and alignment loss (att_maps used by losses.alignment_loss).
"""
import torch
import torch.nn as nn
from torch.nn.utils.parametrizations import weight_norm


class BANLayer(nn.Module):
    def __init__(self, v_dim, q_dim, h_dim, h_out, act='ReLU', dropout=0.2, k=3):
        super().__init__()
        self.c = 32
        self.k = k
        self.v_dim = v_dim
        self.q_dim = q_dim
        self.h_dim = h_dim
        self.h_out = h_out

        self.v_net = FCNet([v_dim, h_dim * self.k], act=act, dropout=dropout)
        self.q_net = FCNet([q_dim, h_dim * self.k], act=act, dropout=dropout)
        if 1 < k:
            self.p_net = nn.AvgPool1d(self.k, stride=self.k)

        if h_out <= self.c:
            self.h_mat = nn.Parameter(torch.Tensor(1, h_out, 1, h_dim * self.k).normal_())
            self.h_bias = nn.Parameter(torch.Tensor(1, h_out, 1, 1).normal_())
        else:
            self.h_net = weight_norm(nn.Linear(h_dim * self.k, h_out), dim=None)

        self.bn = nn.BatchNorm1d(h_dim)

    def attention_pooling(self, v, q, att_map):
        fusion_logits = torch.einsum('bvk,bvq,bqk->bk', (v, att_map, q))
        if 1 < self.k:
            fusion_logits = fusion_logits.unsqueeze(1)
            fusion_logits = self.p_net(fusion_logits).squeeze(1) * self.k
        return fusion_logits

    def forward(self, v, q, softmax=False):
        v_ = self.v_net(v)
        q_ = self.q_net(q)
        if self.h_out <= self.c:
            att_maps = torch.einsum('xhyk,bvk,bqk->bhvq', (self.h_mat, v_, q_)) + self.h_bias
        else:
            v_t = v_.transpose(1, 2).unsqueeze(3)
            q_t = q_.transpose(1, 2).unsqueeze(2)
            d_ = torch.matmul(v_t, q_t)
            d_ = d_.transpose(1, 2).transpose(2, 3)
            att_maps = self.h_net(d_)
            att_maps = att_maps.transpose(2, 3).transpose(1, 2)
        if softmax:
            v_num, q_num = v.size(1), q.size(1)
            p = torch.softmax(att_maps.view(-1, self.h_out, v_num * q_num), 2)
            att_maps = p.view(-1, self.h_out, v_num, q_num)
        logits = self.attention_pooling(v_, q_, att_maps[:, 0, :, :])
        for i in range(1, self.h_out):
            logits = logits + self.attention_pooling(v_, q_, att_maps[:, i, :, :])
        logits = self.bn(logits)
        return logits, att_maps


class FCNet(nn.Module):
    """Non-linear fully connected network (from BAN-VQA)."""

    def __init__(self, dims, act='ReLU', dropout=0):
        super().__init__()
        layers = []
        for i in range(len(dims) - 2):
            if 0 < dropout:
                layers.append(nn.Dropout(dropout))
            layers.append(weight_norm(nn.Linear(dims[i], dims[i + 1], bias=False), dim=None))
            if act:
                layers.append(getattr(nn, act)())
        if 0 < dropout:
            layers.append(nn.Dropout(dropout))
        layers.append(weight_norm(nn.Linear(dims[-2], dims[-1], bias=False), dim=None))
        if act:
            layers.append(getattr(nn, act)())
        self.main = nn.Sequential(*layers)

    def forward(self, x):
        return self.main(x)


class BCNet(nn.Module):
    """Bilinear connect network (kept for reference/ablation)."""

    def __init__(self, v_dim, q_dim, h_dim, h_out, act='ReLU', dropout=(0.2, 0.5), k=3):
        super().__init__()
        self.c = 32
        self.k = k
        self.v_dim = v_dim
        self.q_dim = q_dim
        self.h_dim = h_dim
        self.h_out = h_out
        self.v_net = FCNet([v_dim, h_dim * self.k], act=act, dropout=dropout[0])
        self.q_net = FCNet([q_dim, h_dim * self.k], act=act, dropout=dropout[0])
        self.dropout = nn.Dropout(dropout[1])
        if 1 < k:
            self.p_net = nn.AvgPool1d(self.k, stride=self.k)
        if h_out is None:
            pass
        elif h_out <= self.c:
            self.h_mat = nn.Parameter(torch.Tensor(1, h_out, 1, h_dim * self.k).normal_())
            self.h_bias = nn.Parameter(torch.Tensor(1, h_out, 1, 1).normal_())
        else:
            self.h_net = weight_norm(nn.Linear(h_dim * self.k, h_out), dim=None)

    def forward(self, v, q):
        v_ = self.dropout(self.v_net(v))
        q_ = self.q_net(q)
        if self.h_out <= self.c:
            logits = torch.einsum('xhyk,bvk,bqk->bhvq', (self.h_mat, v_, q_)) + self.h_bias
        else:
            v_t = v_.transpose(1, 2).unsqueeze(3)
            q_t = q_.transpose(1, 2).unsqueeze(2)
            d_ = torch.matmul(v_t, q_t).transpose(1, 2).transpose(2, 3)
            logits = self.h_net(d_).transpose(2, 3).transpose(1, 2)
        return logits
