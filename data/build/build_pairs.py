"""Build MVP positive/negative pairs (v2).

- Peptide seq from complex PDB (shortest/peptide chain) or literature literal.
- Protein seq: complex receptor chain if <=1022aa; else UniProt with a
  biologically-motivated domain window (annotated inline).
- Negatives: shuffled re-pairing (same pools) + shuffled-peptide hard negatives.
"""
import sys
import os
import csv
import random
from collections import defaultdict

sys.stdout.reconfigure(encoding='utf-8')

AA3 = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G',
       'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
       'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V',
       'TRP': 'W', 'TYR': 'Y', 'ASX': 'D', 'GLX': 'E', 'ILM': 'I', 'MVA': 'F',
       'ACE': '', 'NH2': '', 'HOH': '', 'HYP': 'P', 'MLY': 'K', 'MSE': 'M',
       'CSO': 'C', 'CGU': 'E', 'SEP': 'S', 'TPO': 'T', 'PTR': 'Y', 'PCA': 'Q'}

AA1 = set('ACDEFGHIKLMNPQRSTVWY')


def read_fasta(path):
    seqs = {}
    cur = None
    buf = []
    for line in open(path, encoding='utf-8'):
        line = line.strip()
        if line.startswith('>'):
            if cur:
                seqs[cur] = ''.join(buf)
            cur = line[1:].split()[0]
            buf = []
        elif line:
            buf.append(line)
    if cur:
        seqs[cur] = ''.join(buf)
    return seqs


def pdb_seqs(path):
    """{chain: one-letter seq} from SEQRES; skip non-polymer."""
    chains = defaultdict(list)
    for line in open(path, encoding='utf-8', errors='replace'):
        if line.startswith('SEQRES'):
            ch = line[11]
            chains[ch].extend(AA3.get(t, '') for t in line[19:].split())
    return {ch: ''.join(v) for ch, v in chains.items() if v}


def pdb_atom_seqs(path):
    """{chain: seq} from ATOM records (ordered, unique consecutive residues)."""
    chains = defaultdict(list)
    for line in open(path, encoding='utf-8', errors='replace'):
        if line.startswith('ATOM') and line[12:16].strip() == 'CA':
            ch = line[21]
            res = line[17:20].strip()
            aa = AA3.get(res, '')
            if not chains[ch] or chains[ch][-1] != aa:
                chains[ch].append(aa)
    return {ch: ''.join(v) for ch, v in chains.items() if v}


def load_uniprot(uid):
    d = read_fasta(f'mvp_cpu/uniprot/{uid}.fasta')
    return list(d.values())[0] if d else None


def win(seq, a, b=None):
    b = a if b is None else b
    return seq[a - 1:b]


# ---------- catalogue ----------
# (pep_id, pep_spec, prot_spec, evidence, group)
# pep_spec: ('pdb', pdbid, 'SEQRES'|'ATOM', chain) | ('lit', seq)
# prot_spec: ('pdb', pdbid, chain) | ('uni', uid) | ('uni', uid, a, b) window
CAT = [
    ('RES701-3',    ('pdb', '9KDF', 'SEQRES', 'A'), ('pdb', '9KDF', 'R'),
     'cryoEM 9KDF (ETB-RES701-3)', 'RES701'),
    ('RES701-1',    ('lit', 'GNWHGTAPDWFFNYYW'), ('uni', 'P24530'),
     'ETB antagonist biochem', 'RES701'),
    ('PB1m7',       ('lit', 'CNSNVLSWQTYSRWYC'), ('pdb', '7VF3', 'A'),
     'Xray 7VF3 PlexinB1 (peptide from ATOM chain B; disulfide-cyclized)', 'PB1'),
    ('MccJ25',      ('lit', 'GGAGHVPEYFVGIGTPISFYG'), ('uni', 'P0A8T7', 900, 1407),
     'cryoEM 6N60 RNAP secondary channel (RpoC C-term window 900-1407)', 'RNAP'),
    ('Capistruin',  ('lit', 'GTPGFQTPDARVISRFGFN'), ('uni', 'P0A8T7', 900, 1407),
     'cryoEM 6N61', 'RNAP'),
    ('Capi-var1',   ('lit', 'GTPGFQTPDNRVISRFGFN'), ('uni', 'P0A8T7', 900, 1407),
     'LP_318 family homolog', 'RNAP'),
    ('Capi-var2',   ('lit', 'GQPGYQTIDFRVVTRLGGR'), ('uni', 'P0A8T7', 900, 1407),
     'LP_12391 family homolog', 'RNAP'),
    ('Compstatin',  ('pdb', '2QKI', 'SEQRES', 'G'), ('pdb', '2QKI', 'A'),
     'Xray 2QKI C3c-beta chain', 'C3'),
    ('MCoTI-II',    ('pdb', '4GUX', 'SEQRES', None), ('uni', 'P00760'),
     'Xray 4GUX bovine trypsin', 'TRY'),
    ('MCoTI-II-h',  ('pdb', '4GUX', 'SEQRES', None), ('uni', 'P07477'),
     'homolog transfer: human trypsin', 'TRY'),
    ('p53pep',      ('pdb', '1YCR', 'SEQRES', None), ('uni', 'Q00987', 1, 110),
     'Xray 1YCR MDM2 p53-binding domain (1-110)', 'MDM2'),
    ('Cilengitide-b3', ('pdb', '1L5G', 'SEQRES', 'C'), ('pdb', '1L5G', 'B'),
     'Xray 1L5G ITGB3', 'ITG'),
    ('Cilengitide-av', ('pdb', '1L5G', 'SEQRES', 'C'), ('pdb', '1L5G', 'A'),
     'Xray 1L5G ITGAV', 'ITG'),
    ('Lassomycin',  ('pdb', '2MAI', 'SEQRES', 'A'), ('pdb', '8IBO', 'A'),
     '8IBO ClpC1-NTD (lassomycin NMR 2MAI)', 'CLPC'),
    ('Anantin',     ('lit', 'GFIGWGKDIFGHYGG'), ('uni', 'P16066', 1, 450),
     'NPR-A antagonist (ECD window 1-450)', 'ANANTIN'),
    ('Anantin-v2',  ('lit', 'GFIGWGKDIFGHYSGGF'), ('uni', 'P16066', 1, 450),
     'family homolog', 'ANANTIN'),
    ('MccJ25-FhuA', ('lit', 'GGAGHVPEYFVGIGTPISFYG'), ('uni', 'P06971'),
     'uptake transporter 4CU4', 'FhuA'),
]


def main():
    seqres = {}
    atom = {}

    def get_seqres(pdb_id, ch):
        if pdb_id not in seqres:
            p = f'mvp_cpu/pdbs/{pdb_id}.pdb'
            seqres[pdb_id] = pdb_seqs(p) if os.path.isfile(p) else {}
        return seqres[pdb_id].get(ch)

    def get_atom(pdb_id, ch):
        if pdb_id not in atom:
            p = f'mvp_cpu/pdbs/{pdb_id}.pdb'
            atom[pdb_id] = pdb_atom_seqs(p) if os.path.isfile(p) else {}
        return atom[pdb_id].get(ch)

    rows = []
    skipped = []
    for name, pep_spec, prot_spec, ev, grp in CAT:
        # --- peptide ---
        if pep_spec[0] == 'lit':
            pep_seq = pep_spec[1]
        else:
            _, pdb_id, mode, ch = pep_spec
            chains = get_seqres(pdb_id, None) if mode == 'SEQRES' else None
            if ch is None:
                # auto: shortest chain
                src = get_seqres(pdb_id, None) if mode == 'SEQRES' else get_atom(pdb_id, None)
                cands = [(c, s) for c, s in (seqres.get(pdb_id, {}) if mode == 'SEQRES' else atom.get(pdb_id, {})).items()
                         if s and len(s) <= 60 and all(a in AA1 for a in s)]
                pep_seq = min(cands, key=lambda t: len(t[1]))[1] if cands else None
            else:
                pep_seq = (get_seqres(pdb_id, ch) if mode == 'SEQRES' else get_atom(pdb_id, ch))
        # --- protein ---
        if prot_spec[0] == 'pdb':
            _, pdb_id, ch = prot_spec
            prot_seq = get_seqres(pdb_id, ch)
            prot_id = f'{pdb_id}:{ch}'
        elif prot_spec[0] == 'uni':
            uid = prot_spec[1]
            prot_seq = load_uniprot(uid)
            prot_id = uid
            if len(prot_spec) == 4:
                prot_seq = win(prot_seq, prot_spec[2], prot_spec[3])
                prot_id = f'{uid}[{prot_spec[2]}:{prot_spec[3]}]'
        if not pep_seq or not prot_seq:
            skipped.append((name, 'missing seq'))
            continue
        if any(a not in AA1 for a in pep_seq):
            skipped.append((name, f'non-std pep: {pep_seq}'))
            continue
        rows.append(dict(pep_id=name, group=grp, pep_seq=pep_seq, prot_id=prot_id,
                         prot_seq=prot_seq, label=1, evidence=ev))
        print(f'+ {name}: pep {len(pep_seq)}aa x prot {len(prot_seq)}aa ({ev})')

    print(f'\npositives: {len(rows)}, skipped: {skipped}')

    # ---------- negatives ----------
    random.seed(42)
    pos_set = {(r['pep_seq'], r['prot_seq']) for r in rows}
    neg_rows = []
    attempts = 0
    target = len(rows) * 4
    while len(neg_rows) < target and attempts < 50000:
        attempts += 1
        a = random.choice(rows)
        b = random.choice(rows)
        if a is b:
            continue
        if (a['pep_seq'], b['prot_seq']) in pos_set:
            continue
        key = (a['pep_seq'], b['prot_seq'])
        if key in {(n['pep_seq'], n['prot_seq']) for n in neg_rows}:
            continue
        neg_rows.append(dict(pep_id=f'shuf-{len(neg_rows)}', group='NEG',
                             pep_seq=a['pep_seq'], prot_id=b['prot_id'],
                             prot_seq=b['prot_seq'], label=0, evidence='shuffled-pair'))
    # hard negatives: shuffled peptide sequence, same protein
    for r in rows:
        s = list(r['pep_seq'])
        random.shuffle(s)
        s = ''.join(s)
        if s == r['pep_seq']:
            s = s[::-1]
        neg_rows.append(dict(pep_id=f'scram-{r["pep_id"]}', group='NEG',
                             pep_seq=s, prot_id=r['prot_id'], prot_seq=r['prot_seq'],
                             label=0, evidence='scrambled-peptide'))
    print(f'negatives: {len(neg_rows)} ({target} shuffled + {len(rows)} scrambled)')

    all_rows = rows + neg_rows
    with open('mvp_cpu/pairs.csv', 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=['pep_id', 'group', 'pep_seq', 'prot_id', 'prot_seq', 'label', 'evidence'])
        w.writeheader()
        for r in all_rows:
            w.writerow(r)
    print('saved mvp_cpu/pairs.csv')
    from collections import Counter
    print('labels:', Counter(r['label'] for r in all_rows))
    print('pos groups:', Counter(r['group'] for r in rows))


if __name__ == '__main__':
    main()
