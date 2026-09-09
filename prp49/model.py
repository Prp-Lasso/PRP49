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
                 dual_head=False):
        super().__init__()
        self.esm = AutoModel.from_pretrained(esm_model_name)
        self.esm_emb_dim = self.esm.config.hidden_size
        self.lassoesm = AutoModel.from_pretrained(lassoesm_model_name)
        self.lassoesm_emb_dim = self.lassoesm.config.hidden_size
        self.dual_head = dual_head

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
            BANLayer(v_dim=self.esm_emb_dim, q_dim=self.lassoesm_emb_dim,
                     h_dim=h_dim, h_out=n_heads, k=k, dropout=ban_dropout),
            name='h_mat', dim=None)

        self.fc1 = nn.Linear(h_dim, 512)
        self.fc2 = nn.Linear(512, 128)
        self.out = nn.Linear(128, 1)
        if dual_head:
            self.aff_fc = nn.Linear(h_dim, 64)
            self.aff_out = nn.Linear(64, 1)
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

    def forward(self, seq1_ids, seq1_mask, seq2_ids, seq2_mask):
        o1 = self.esm(input_ids=seq1_ids, attention_mask=seq1_mask)
        o2 = self.lassoesm(input_ids=seq2_ids, attention_mask=seq2_mask)
        f1 = self._extract_aa(o1.last_hidden_state, seq1_ids)   # protein
        f2 = self._extract_aa(o2.last_hidden_state, seq2_ids)   # peptide
        fused, att_maps = self.ban(f1, f2)                     # (B, h_dim), (B, heads, L1, L2)
        x = self.relu(self.fc1(fused))
        x = self.dropout(x)
        x = self.relu(self.fc2(x))
        x = self.dropout(x)
        logits = self.out(x).squeeze(-1)
        affinity = None
        if self.dual_head:
            a = self.relu(self.aff_fc(fused))
            affinity = self.aff_out(a).squeeze(-1)
        return logits, att_maps, affinity
