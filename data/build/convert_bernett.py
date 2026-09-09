"""Convert TUnA Bernett data (Intra0/1/2 interaction tsv) to PRP49 pairs format.

Bernett rows: protA<TAB>protB<TAB>label, with sequences in
human_swissprot_oneliner.fasta (header = UniProt id). We map protA -> prot_seq
(protein encoder side) and protB -> pep_seq (peptide encoder side); for L0
warmup both encoders are ESM-2 (see config_warmup.yaml).

Usage: python scratch/convert_bernett.py [--max_pairs 50000]
"""
import argparse
import os
import random
import re

root = r'D:\deepseek_harness\prp49'
random.seed(42)


def load_fasta(path):
    seqs = {}
    cur = None
    for line in open(path, encoding='utf-8', errors='ignore'):
        line = line.strip()
        if line.startswith('>'):
            cur = line[1:].split()[0]
            seqs[cur] = []
        elif cur:
            seqs[cur].append(line)
    return {k: ''.join(v) for k, v in seqs.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--max_pairs', type=int, default=50000)
    ap.add_argument('--intras', default='0,1,2')
    ap.add_argument('--out', default=r'D:\deepseek_harness\prp49\mvp_cpu\bernett_pairs.csv')
    args = ap.parse_args()

    fasta = load_fasta(os.path.join(root, 'TUnA', 'data', 'raw', 'bernett',
                                   'human_swissprot_oneliner.fasta'))
    print(f'fasta: {len(fasta)} sequences')

    rows = []
    n_missing = 0
    for intra in args.intras.split(','):
        p = os.path.join(root, 'TUnA', 'data', 'processed', 'bernett',
                         f'Intra{intra.strip()}_interaction.tsv')
        if not os.path.exists(p):
            print(f'skip missing {p}')
            continue
        for line in open(p):
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            a, b, lab = parts[0], parts[1], parts[2]
            if a not in fasta or b not in fasta:
                n_missing += 1
                continue
            rows.append((a, b, fasta[a], fasta[b], lab))
        print(f'Intra{intra.strip()}: loaded (total rows now {len(rows)})')
    print(f'missing-sequence rows skipped: {n_missing}')

    random.shuffle(rows)
    if args.max_pairs and len(rows) > args.max_pairs:
        rows = rows[:args.max_pairs]
    pos = sum(1 for r in rows if r[4] == '1')
    print(f'writing {len(rows)} pairs (pos {pos}, neg {len(rows) - pos})')

    with open(args.out, 'w') as f:
        f.write('pep_id,pdb,pep_seq,prot_seq,prot_id,label\n')
        for a, b, sa, sb, lab in rows:
            # protB -> peptide side id; pdb='bernett'
            f.write(f'{b},bernett,{sb},{sa},{a},{lab}\n')
    print(f'saved {args.out}')


if __name__ == '__main__':
    main()
