"""Data loading: pairs CSV, alignment labels, ablation set."""
import csv
import json
import os
import numpy as np
import pandas as pd
import torch
from torch.utils.data import Dataset


class PairDataset(Dataset):
    """(protein_seq, peptide_seq, label) pairs tokenized on the fly.

    label_type 'binary' -> label column; 'energy' -> energy column (regression).
    """

    def __init__(self, df, tokenizer, max_len=1024, label_type='binary', align_map=None):
        self.df = df.reset_index(drop=True)
        self.tok = tokenizer
        self.max_len = max_len
        self.label_type = label_type
        self.align_map = align_map or {}
        self._align_cache = {}

    def __len__(self):
        return len(self.df)

    def __getitem__(self, i):
        row = self.df.iloc[i]
        prot = row['prot_seq']
        pep = row['pep_seq']
        t1 = self.tok(prot, return_tensors='pt', max_length=self.max_len,
                      truncation=True, padding='max_length')
        t2 = self.tok(pep, return_tensors='pt', max_length=self.max_len,
                      truncation=True, padding='max_length')
        if self.label_type == 'binary':
            lab = torch.tensor(float(row['label']), dtype=torch.float32)
        else:
            lab = torch.tensor(float(row['energy']), dtype=torch.float32)
        return {
            'seq1_ids': t1['input_ids'][0], 'seq1_mask': t1['attention_mask'][0],
            'seq2_ids': t2['input_ids'][0], 'seq2_mask': t2['attention_mask'][0],
            'label': lab,
            'key': f"{row.get('pep_id', i)}|{row.get('pdb', '')}",
            'align': self._match_align(row)
        }

    def _match_align(self, row):
        """Return (prot_len, pep_len) contact matrix for this pair, or None.

        Structural chains and training sequences can differ by unresolved
        residues, windowing (full chain vs domain window) or mid-peptide
        insertions, so both sides are mapped with SequenceMatcher.
        """
        if not self.align_map:
            return None
        key = (row['pep_seq'], row['prot_seq'])
        if key in self._align_cache:
            return self._align_cache[key]
        mat = self._align_match(row['pep_seq'], row['prot_seq'])
        self._align_cache[key] = mat
        return mat

    def _align_match(self, pep, prot):
        import difflib
        for (pep_ref, prot_ref), mat0 in self.align_map.items():
            # ---- peptide side: pair pep index -> structural pep index ----
            if pep == pep_ref:
                col_map = np.arange(len(pep_ref))
            else:
                smp = difflib.SequenceMatcher(None, pep, pep_ref, autojunk=False)
                pblocks = [b for b in smp.get_matching_blocks() if b.size >= 1]
                if sum(b.size for b in pblocks) < 0.9 * min(len(pep), len(pep_ref)):
                    continue
                col_map = np.full(len(pep), -1, dtype=np.int64)
                for a, b, n in pblocks:
                    col_map[a:a + n] = np.arange(b, b + n)
            # ---- protein side: pair prot index -> structural prot index ----
            smq = difflib.SequenceMatcher(None, prot, prot_ref, autojunk=False)
            qblocks = [b for b in smq.get_matching_blocks() if b.size >= 1]
            if sum(b.size for b in qblocks) < 0.9 * min(len(prot), len(prot_ref)):
                continue
            M = np.zeros((len(prot), len(pep)), dtype=np.float32)
            ok_cols = col_map >= 0
            for a, b, n in qblocks:          # prot[a:a+n] == prot_ref[b:b+n]
                sub = mat0.T[b:b + n, :]     # (n, pep_ref_len) structural contacts
                sub2 = np.zeros((n, len(pep)), dtype=np.float32)
                sub2[:, ok_cols] = sub[:, col_map[ok_cols]]
                M[a:a + n, :] = sub2
            return M
        return None


def collate_fn(batch):
    import torch
    seq1_ids = torch.stack([b['seq1_ids'] for b in batch])
    seq1_mask = torch.stack([b['seq1_mask'] for b in batch])
    seq2_ids = torch.stack([b['seq2_ids'] for b in batch])
    seq2_mask = torch.stack([b['seq2_mask'] for b in batch])
    labels = torch.stack([b['label'] for b in batch])
    keys = [b['key'] for b in batch]
    aligns = [b.get('align') for b in batch]
    return {'seq1_ids': seq1_ids, 'seq1_mask': seq1_mask,
            'seq2_ids': seq2_ids, 'seq2_mask': seq2_mask,
            'label': labels, 'keys': keys, 'align': aligns}


def load_pairs(csv_path, **kw):
    df = pd.read_csv(csv_path)
    df['prot_seq'] = df['prot_seq'].astype(str).str.replace('*', '', regex=False)
    df['pep_seq'] = df['pep_seq'].astype(str).str.replace('*', '', regex=False)
    df = df[df['prot_seq'].str.len() > 0]
    df = df[df['pep_seq'].str.len() > 0]
    if 'label' in df.columns:
        df['label'] = df['label'].astype(float)
    return df


def load_alignment_labels(json_path):
    """Load interface contact labels as {(pep_ref, prot_ref): mat(pep_len x prot_len)}.

    contacts.json holds residue-pair lists per complex; the companion npz files
    (written by interface_extract.py) store the binary contact matrix plus the
    structural chain sequences, which are used to match training pairs.
    """
    if not json_path or not os.path.exists(json_path):
        return {}
    with open(json_path, encoding='utf-8') as f:
        data = json.load(f)
    base = os.path.dirname(json_path)
    align_map = {}
    for pdb, entry in data.items():
        for pc, recs in entry['contacts'].items():
            for rc in recs:
                npz_path = os.path.join(base, f'{pdb}_{pc}_{rc}.npz')
                if not os.path.exists(npz_path):
                    continue
                z = np.load(npz_path, allow_pickle=True)
                mat = z['mat']
                pep_ref = ''.join(str(c) for c in z['pep_seq'])
                prot_ref = ''.join(str(c) for c in z['prot_seq'])
                align_map[(pep_ref, prot_ref)] = mat.astype(np.float32)
    return align_map
