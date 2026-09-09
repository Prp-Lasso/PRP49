"""Generate augmentation variants + ring-opening ablation controls.

Rules (based on LassoPred annotation: ring Gly1-(E/D at ring_len)):
  - conservative variants: substitutions in ring/loop/tail (never at isopeptide site)
  - ring-opening (ablation): isopeptide acceptor E/D -> S at position ring_len
  - distal control: conservative substitution far from ring (tail)
Outputs:
  mvp_cpu/augmented_positives.csv  (training augmentation)
  mvp_cpu/ablation_set.csv         (ring-open + distal control, for evaluation)
"""
import sys
import csv
from collections import defaultdict

sys.stdout.reconfigure(encoding='utf-8')

# (pep_id, seq, isopeptide_pos (1-based; None if disulfide), ring_len, note)
LASSO_PEPS = [
    ('RES701-3',    'GNWHGTSPDWFFNYYW',      9, 9, 'LP_434'),
    ('RES701-1',    'GNWHGTAPDWFFNYYW',      9, 9, 'biochem'),
    ('MccJ25',      'GGAGHVPEYFVGIGTPISFYG', 8, 8, 'LP_471'),
    ('Capistruin',  'GTPGFQTPDARVISRFGFN',   9, 9, 'LP_105'),
    ('Capi-var1',   'GTPGFQTPDNRVISRFGFN',   9, 9, 'LP_318'),
    ('Capi-var2',   'GQPGYQTIDFRVVTRLGGR',   9, 9, 'LP_12391'),
    ('Anantin',     'GFIGWGKDIFGHYGG',       8, 8, 'LP_366'),
    ('Anantin-v2',  'GFIGWGKDIFGHYSGGF',     8, 8, 'LP_388'),
    ('Lassomycin',  'GLRRLFADQLVGRNI',       8, 8, '8IBO chain C (D8 acceptor)'),
    ('PB1m7',       'CNSNVLSWQTYSRWYC',      None, None, 'disulfide C1-C16'),
]

CONS = {'E': ['D'], 'D': ['E'], 'I': ['L', 'V'], 'L': ['I', 'V'], 'V': ['I', 'L'],
        'S': ['T'], 'T': ['S'], 'F': ['Y'], 'Y': ['F'], 'N': ['Q'], 'Q': ['N'],
        'K': ['R'], 'R': ['K'], 'A': ['G'], 'M': ['L']}


def mutate(seq, pos1, new_aa):
    p = pos1 - 1
    return seq[:p] + new_aa + seq[p + 1:]


def generate():
    aug_rows = []
    abl_rows = []
    for pid, seq, iso_pos, ring_len, note in LASSO_PEPS:
        n = len(seq)
        # positions: ring = 1..ring_len, loop = ring_len+1 .. n-2, tail = n-1..n (rough)
        if ring_len:
            ring_pos = list(range(1, ring_len + 1))
            loop_pos = list(range(ring_len + 1, n + 1))
        else:
            ring_pos = [1, n]  # disulfide cysteines
            loop_pos = list(range(2, n))

        # 1) conservative variants: pick 3 positions not at iso site, apply 1-2 substitutions
        made = 0
        for pos in ring_pos + loop_pos:
            if pos == iso_pos:
                continue
            aa = seq[pos - 1]
            for naa in CONS.get(aa, []):
                if made >= 4:
                    break
                var = mutate(seq, pos, naa)
                if var != seq:
                    aug_rows.append(dict(pep_id=f'{pid}-v{pos}{aa}{naa}', orig=pid,
                                         seq=var, kind='conservative', pos=pos,
                                         mut=f'{aa}{pos}{naa}', note=note))
                    made += 1
            if made >= 4:
                break

        # 2) ring-opening ablation (isopeptide) or disulfide break
        if iso_pos:
            aa = seq[iso_pos - 1]
            if aa in 'ED':
                open_seq = mutate(seq, iso_pos, 'S')
                abl_rows.append(dict(pep_id=f'{pid}-open', orig=pid, seq=open_seq,
                                     kind='ring-open', pos=iso_pos,
                                     mut=f'{aa}{iso_pos}S', note=note))
            else:
                print(f'  [warn] {pid}: iso_pos aa={aa} not E/D, skip ring-open')
        else:
            open_seq = mutate(mutate(seq, 1, 'S'), n, 'S')
            abl_rows.append(dict(pep_id=f'{pid}-open', orig=pid, seq=open_seq,
                                 kind='disulfide-break', pos=f'1,{n}',
                                 mut='C1S,CnS', note=note))

        # 3) distal control: conservative substitution at tail (far from ring)
        tail_cand = [p for p in range(max(ring_len or 3, 3), n + 1) if p != iso_pos]
        done = False
        for pos in reversed(tail_cand):
            aa = seq[pos - 1]
            for naa in CONS.get(aa, []):
                ctrl = mutate(seq, pos, naa)
                if ctrl != seq:
                    abl_rows.append(dict(pep_id=f'{pid}-ctrl', orig=pid, seq=ctrl,
                                         kind='distal-control', pos=pos,
                                         mut=f'{aa}{pos}{naa}', note=note))
                    done = True
                    break
            if done:
                break

    with open('mvp_cpu/augmented_positives.csv', 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=['pep_id', 'orig', 'seq', 'kind', 'pos', 'mut', 'note'])
        w.writeheader()
        for r in aug_rows:
            w.writerow(r)
    with open('mvp_cpu/ablation_set.csv', 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=['pep_id', 'orig', 'seq', 'kind', 'pos', 'mut', 'note'])
        w.writeheader()
        for r in abl_rows:
            w.writerow(r)

    from collections import Counter
    print(f'augmented variants: {len(aug_rows)} | kinds: {Counter(r["kind"] for r in aug_rows)}')
    print(f'ablation/control: {len(abl_rows)} | kinds: {Counter(r["kind"] for r in abl_rows)}')
    for r in abl_rows:
        print(f'  {r["pep_id"]:22s} {r["kind"]:16s} {r["mut"]:8s} {r["seq"]}')


if __name__ == '__main__':
    generate()
