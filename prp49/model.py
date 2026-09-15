"""Model: protein ESM + peptide LassoESM + BAN fusion + (optional) dual head."""
import torch
import torch.nn as nn
from torch.nn.utils.parametrizations import weight_norm
from transformers import AutoModel

from .ban import BANLayer


class InteractionPredictor(nn.Module):
    """ESM-2 (protein) x LassoESM (peptide) with bilinear attention fusion.

    forward returns (logits, att_maps, affinity)
      logits: (B,) binary binding score
      att_maps: (B, heads, prot_len, pep_len) attention maps for alignment loss
      affinity: (B,) optional regression head (only when dual_head=True)
    """

    def __init__(self, esm_model_name='facebook/esm2_t12_35M_UR50D',
                 lassoesm_model_name='./LassoESM',
                 h_dim=256, n_heads=3, k=2,
                 dropout=0.3, ban_dropout=0.3,
                 freeze_esm=False, lassoesm_unfreeze=0,
                 dual_head=False, use_topo=False, topo_dim=16,
                 affinity_bound=None, base_head=False):
        super().__init__()
        self.esm = AutoModel.from_pretrained(esm_model_name)
        self.esm_emb_dim = self.esm.config.hidden_size
        self.lassoesm = AutoModel.from_pretrained(lassoesm_model_name)
        self.lassoesm_emb_dim = self.lassoesm.config.hidden_size
        self.dual_head = dual_head
        # Regression stabilisers (see docs/PRP49_回归RMSE约束方案.md):
        #   affinity_bound=(lo,hi) squashes the prediction into the label range via
        #     tanh, so the unbounded linear head can no longer diverge (it produced
        #     RMSE 6.86 against a 1.78 mean-predictor baseline).
        #   base_head=True adds a target-only baseline term:
        #     affinity = base(target_seq) + delta(pep, target), letting the target's
        #     typical affinity level be learned from its sequence rather than
        #     extrapolated by the interaction term.
        self.affinity_bound = affinity_bound
        self.base_head = base_head
        # II2: explicit topology channel (ring/loop/tail) on the peptide side
        self.use_topo = use_topo
        self.topo_emb = nn.Embedding(5, topo_dim) if use_topo else None
        q_dim = self.lassoesm_emb_dim + (topo_dim if use_topo else 0)

        if freeze_esm:
            for p in self.esm.parameters():
                p.requires_grad = False
        for p in self.lassoesm.parameters():
            p.requires_grad = False
        if lassoesm_unfreeze == -1:
            for p in self.lassoesm.parameters():
                p.requires_grad = True
        elif lassoesm_unfreeze > 0:
            layers = getattr(getattr(self.lassoesm, 'encoder', None), 'layer', None)
            if layers is not None:
                for layer in layers[-lassoesm_unfreeze:]:
                    for p in layer.parameters():
                        p.requires_grad = True

        self.ban = weight_norm(
            BANLayer(v_dim=self.esm_emb_dim, q_dim=q_dim,
                     h_dim=h_dim, h_out=n_heads, k=k, dropout=ban_dropout),
            name='h_mat', dim=None)

        self.fc1 = nn.Linear(h_dim, 512)
        self.fc2 = nn.Linear(512, 128)
        self.out = nn.Linear(128, 1)
        if dual_head:
            self.aff_fc = nn.Linear(h_dim, 64)
            self.aff_out = nn.Linear(64, 1)
        if dual_head and base_head:
            # target-only baseline: pooled protein embedding -> scalar
            self.base_fc = nn.Linear(self.esm_emb_dim, 64)
            self.base_out = nn.Linear(64, 1)
        self.relu = nn.ReLU()
        self.dropout = nn.Dropout(dropout)

    @staticmethod
    def _extract_aa(tok_emb, input_ids):
        """Drop special tokens (<cls>/<eos>/<pad>), pad to batch max."""
        aa_mask = input_ids >= 4
        rows = []
        for i in range(tok_emb.size(0)):
            rows.append(tok_emb[i][aa_mask[i]])
        max_len = max(r.size(0) for r in rows) if rows else 1
        dim = tok_emb.size(-1)
        padded = torch.zeros(len(rows), max_len, dim, device=tok_emb.device)
        for i, r in enumerate(rows):
            padded[i, :r.size(0)] = r
        return padded

    def forward(self, seq1_ids, seq1_mask, seq2_ids, seq2_mask, seq2_topo=None):
        o1 = self.esm(input_ids=seq1_ids, attention_mask=seq1_mask)
        o2 = self.lassoesm(input_ids=seq2_ids, attention_mask=seq2_mask)
        f1 = self._extract_aa(o1.last_hidden_state, seq1_ids)   # protein
        f2 = self._extract_aa(o2.last_hidden_state, seq2_ids)   # peptide
        if self.use_topo and seq2_topo is not None:
            # map per-residue topology ids into the same residue space as f2
            aa_mask = seq2_ids >= 4
            rows = [seq2_topo[i][aa_mask[i]] for i in range(seq2_topo.size(0))]
            L2 = f2.size(1)
            topo_pad = torch.zeros(f2.size(0), L2, dtype=torch.long, device=f2.device)
            for i, r in enumerate(rows):
                n = min(r.numel(), L2)
                topo_pad[i, :n] = r[:n]
            f2 = torch.cat([f2, self.topo_emb(topo_pad)], dim=-1)
        fused, att_maps = self.ban(f1, f2)                     # (B, h_dim), (B, heads, L1, L2)
        x = self.relu(self.fc1(fused))
        x = self.dropout(x)
        x = self.relu(self.fc2(x))
        x = self.dropout(x)
        logits = self.out(x).squeeze(-1)
        affinity = None
        if self.dual_head:
            a = self.relu(self.aff_fc(fused))
            raw = self.aff_out(a).squeeze(-1)
            if self.base_head:
                # mean-pool the protein embedding (respecting the attention mask)
                m = seq1_mask.unsqueeze(-1).float()
                tgt = (o1.last_hidden_state * m).sum(1) / m.sum(1).clamp(min=1.0)
                raw = raw + self.base_out(self.relu(self.base_fc(tgt))).squeeze(-1)
            if self.affinity_bound is not None:
                lo, hi = self.affinity_bound
                mid, half = (lo + hi) / 2.0, (hi - lo) / 2.0
                affinity = mid + half * torch.tanh(raw / half)
            else:
                affinity = raw
        return logits, att_maps, affinity
