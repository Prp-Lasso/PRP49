"""Extract interface contact residue pairs from complex PDBs (for L_align labels).

Rules:
- peptide chain = shortest polymer chain (5-40 aa, standard AAs)
- receptor chain = any other polymer chain with >=1 heavy-atom contact (<5.0 A) to peptide
- output contacts.json: {pdb: {receptor_chain: [[pep_idx, prot_idx, min_dist], ...]}}
- also save a normalized per-pair contact matrix file per complex (npz) for direct use
"""
import sys
import os
import json
import numpy as np
from collections import defaultdict

sys.stdout.reconfigure(encoding='utf-8')

CUTOFF = 5.0  # heavy-atom contact distance

AA3 = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G',
       'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
       'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V',
       'TRP': 'W', 'TYR': 'Y'}
AA1 = set('ACDEFGHIKLMNPQRSTVWY')

# complexes available locally: (pdb, [peptide_chain_candidates])
COMPLEXES = {
    '9KDF': ['A'], '7VF3': ['B'], '6N60': None, '6N61': None, '8IBO': ['C'],
    '1YCR': None, '4GUX': None, '1L5G': ['C'], '2QKI': ['G', 'H'], '9KDG': [],
}


def parse_pdb(path):
    """Return {chain: {resseq: {'atoms': [(name, x,y,z)], 'aa': one-letter}}}"""
    chains = defaultdict(lambda: defaultdict(dict))
    for line in open(path, encoding='utf-8', errors='replace'):
        if line.startswith('ATOM') or line.startswith('HETATM'):
            ch = line[21]
            resn = line[17:20].strip()
            aa = AA3.get(resn, '')
            resseq = int(line[22:26])
            aname = line[12:16].strip()
            try:
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            except ValueError:
                continue
            if aname.startswith('H'):
                continue
            r = chains[ch][resseq]
            r.setdefault('atoms', []).append((aname, x, y, z))
            r['aa'] = aa
    return chains


def chain_seq(chains, ch):
    r = chains.get(ch, {})
    return ''.join(r[i]['aa'] for i in sorted(r))


def main():
    outdir = 'mvp_cpu/contacts'
    os.makedirs(outdir, exist_ok=True)
    results = {}
    for pdb, pep_cands in COMPLEXES.items():
        path = f'mvp_cpu/pdbs/{pdb}.pdb'
        if not os.path.isfile(path):
            print(f'{pdb}: no PDB file, skip')
            continue
        chains = parse_pdb(path)
        # pick peptide chains: specified or shortest polymer chain 5-40aa
        if pep_cands:
            pep_chains = [c for c in pep_cands if c in chains]
        else:
            cands = []
            for ch in chains:
                s = chain_seq(chains, ch)
                if 5 <= len(s) <= 40 and all(a in AA1 for a in s):
                    cands.append((len(s), ch))
            pep_chains = [ch for _, ch in sorted(cands)[:2]]
        # receptors: chains with contacts to peptide
        contacts = {}
        pep_seq = {}
        for pc in pep_chains:
            pep_seq[pc] = chain_seq(chains, pc)
            pres = sorted(chains[pc])
            for rc in chains:
                if rc == pc or len(chain_seq(chains, rc)) <= 40:
                    continue
                rres = sorted(chains[rc])
                pair_list = []
                for pr in pres:
                    for rr in rres:
                        mind = min(
                            ((ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2) ** 0.5
                            for _, ax, ay, az in chains[pc][pr]['atoms']
                            for _, bx, by, bz in chains[rc][rr]['atoms']
                        ) if chains[pc][pr]['atoms'] and chains[rc][rr]['atoms'] else 1e9
                        if mind <= CUTOFF:
                            pair_list.append([pr, rr, round(mind, 2)])
                if pair_list:
                    contacts.setdefault(pc, {})[rc] = pair_list
        if contacts:
            results[pdb] = dict(pep_sequences={c: s for c, s in pep_seq.items() if c in contacts},
                                contacts=contacts)
            print(f'{pdb}: peptide chains {list(contacts.keys())}, '
                  f'receptors {sorted({rc for v in contacts.values() for rc in v})}, '
                  f'total pairs {sum(len(v) for v in [vv for v in contacts.values() for vv in v.values()])}')

    out = f'{outdir}/contacts.json'
    json.dump(results, open(out, 'w'), indent=1)
    print(f'saved {out}')

    # also save per-complex contact matrices (pep_len x prot_len binary)
    for pdb, data in results.items():
        for pc, recs in data['contacts'].items():
            plen = len(data['pep_sequences'][pc])
            for rc, pairs in recs.items():
                # prot seq from pdb parse
                chains = parse_pdb(f'mvp_cpu/pdbs/{pdb}.pdb')
                pseq = chain_seq(chains, pc)
                rseq = chain_seq(chains, rc)
                mat = np.zeros((len(pseq), len(rseq)), dtype=np.float32)
                for pr, rr, d in pairs:
                    # map resseq to index via sorted position
                    pass
                # rebuild with index mapping
                pres = sorted(chains[pc])
                rres = sorted(chains[rc])
                pidx = {r: i for i, r in enumerate(pres)}
                ridx = {r: i for i, r in enumerate(rres)}
                mat = np.zeros((len(pres), len(rres)), dtype=np.float32)
                for pr, rr, d in pairs:
                    mat[pidx[pr], ridx[rr]] = 1.0
                np.savez_compressed(f'{outdir}/{pdb}_{pc}_{rc}.npz',
                                    mat=mat, pep_seq=np.array(list(pseq)), prot_seq=np.array(list(rseq)))
                print(f'saved matrix {pdb}_{pc}_{rc}: {mat.shape}, contacts {int(mat.sum())}')


if __name__ == '__main__':
    main()
