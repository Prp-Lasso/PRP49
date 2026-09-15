"""Map the experimental docking box onto AlphaFold coordinates - pure Python 3.

Replaces the PyMOL route, which failed twice over: the pymol env is Python 2 (no
f-strings, issue #21) and its licence expired 2024-05-01.

Method (no external tools):
  1. parse CA coordinates + sequence from both structures
  2. align the two SEQUENCES (difflib) - the experiment PDB may not start at
     residue 1, while the AFDB model always does
  3. Kabsch superposition over the matched CA pairs -> rotation R, translation t
  4. box centre' = R @ centre + t ; keep size, floor the edge at 38 A
  5. report RMSD and the number of matched residues so a bad fit is visible

Run with the prp49 env (Python 3 + numpy).
"""
import csv
import os
import sys
from difflib import SequenceMatcher

import numpy as np

BASE = os.path.expanduser('~/LassoPep')
AF_DIR = os.path.join(BASE, 'docking_af', 'af')
REC_DIR = os.path.join(BASE, 'docking', 'rec')
TASKS = os.path.join(BASE, 'docking', 'tasks.csv')
OUT_TASKS = os.path.join(BASE, 'docking_af', 'tasks_af.csv')
OUT_REPORT = os.path.join(BASE, 'docking_af', 'map_report.csv')

AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
       'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
       'THR':'T','TRP':'W','TYR':'Y','VAL':'V','MSE':'M','SEC':'U','PYL':'O'}


def parse_ca(path):
    """Return (sequence, {seq_index: xyz}) using CA atoms in file order."""
    seq, coords = [], []
    for line in open(path, errors='ignore'):
        if line.startswith('ATOM') and line[12:16].strip() == 'CA':
            aa = AA3.get(line[17:20].strip())
            if aa is None:
                continue
            try:
                xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                continue
            seq.append(aa)
            coords.append(xyz)
    return ''.join(seq), np.array(coords, dtype=float)


def kabsch(P, Q):
    """Rotation+translation taking P onto Q (both Nx3)."""
    Pc, Qc = P.mean(0), Q.mean(0)
    H = (P - Pc).T @ (Q - Qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    t = Qc - R @ Pc
    return R, t


# Radius around the old box centre used for the LOCAL superposition.
# A global fit fails on multi-domain proteins (PLXNB1 41.7 A, C3 27.5 A RMSD),
# because the AF monomer and the experimental complex adopt different quaternary
# arrangements. Only the residues lining the binding site need to agree, so the
# fit is restricted to atoms near the old box centre.
LOCAL_RADIUS = 15.0


def main():
    centres = {}
    with open(TASKS) as f:
        for r in csv.DictReader(f):
            centres.setdefault(r['rec'], (float(r['cx']), float(r['cy']), float(r['cz']),
                                          float(r['sx']), float(r['sy']), float(r['sz'])))
    rows = []
    for target, (cx, cy, cz, sx, sy, sz) in sorted(centres.items()):
        af_p, exp_p = os.path.join(AF_DIR, f'{target}.pdb'), os.path.join(REC_DIR, f'{target}.pdb')
        if not (os.path.exists(af_p) and os.path.exists(exp_p)):
            rows.append(dict(target=target, status='missing'))
            print(f'{target:9s} MISSING file')
            continue
        af_seq, af_xyz = parse_ca(af_p)
        ex_seq, ex_xyz = parse_ca(exp_p)
        sm = SequenceMatcher(None, ex_seq, af_seq, autojunk=False)
        pairs = [(i, j) for blk in sm.get_matching_blocks() for i, j in
                 zip(range(blk.a, blk.a + blk.size), range(blk.b, blk.b + blk.size))]
        if len(pairs) < 30:
            rows.append(dict(target=target, status=f'few matched residues ({len(pairs)})'))
            print(f'{target:9s} only {len(pairs)} matched residues - skip')
            continue

        # --- LOCAL fit: keep only matched residues within LOCAL_RADIUS of the old box
        centre = np.array([cx, cy, cz])
        local = [(i, j) for i, j in pairs
                 if np.linalg.norm(ex_xyz[i] - centre) <= LOCAL_RADIUS]
        if len(local) < 20:
            # box centre sits in a gap (missing residues) - fall back to the global fit
            local = pairs
            fit_kind = 'global(fallback)'
        else:
            fit_kind = f'local<={LOCAL_RADIUS:.0f}A'
        P = np.array([ex_xyz[i] for i, _ in local])
        Q = np.array([af_xyz[j] for _, j in local])
        R, t = kabsch(P, Q)
        fitted = (R @ P.T).T + t
        rmsd = float(np.sqrt(((fitted - Q) ** 2).sum(1).mean()))
        # displacement of the box centre is what actually matters for docking
        new_c = R @ centre + t
        shift = float(np.linalg.norm(new_c - centre))
        nsx, nsy, nsz = max(sx, 38.0), max(sy, 38.0), max(sz, 38.0)
        rows.append(dict(target=target, status='ok', rmsd=round(rmsd, 2), fit=fit_kind,
                         matched_residues=len(local), local_rmsd=round(rmsd, 2),
                         exp_len=len(ex_seq), af_len=len(af_seq),
                         cos_old=f'{cx},{cy},{cz}', cx=round(float(new_c[0]), 1),
                         cy=round(float(new_c[1]), 1), cz=round(float(new_c[2]), 1),
                         sx=round(nsx), sy=round(nsy), sz=round(nsz)))
        print(f'{target:9s} {fit_kind:16s} n={len(local):5d} RMSD {rmsd:6.2f} A | '
              f'box moves {shift:7.1f} A -> ({new_c[0]:7.1f},{new_c[1]:7.1f},{new_c[2]:7.1f})')

    with open(OUT_REPORT, 'w', newline='') as f:
        keys = sorted({k for r in rows for k in r})
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)

    ok = [r for r in rows if r.get('status') == 'ok']
    with open(TASKS) as f:
        peps = sorted({r['pep'] for r in csv.DictReader(f)})
    with open(OUT_TASKS, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['id', 'pep', 'rec', 'cx', 'cy', 'cz', 'sx', 'sy', 'sz'])
        tid = 7000
        for r in ok:
            for p in peps:
                w.writerow([tid, p, f'af_{r["target"]}', r['cx'], r['cy'], r['cz'],
                            r['sx'], r['sy'], r['sz']])
                tid += 1
    print(f'\nmapped {len(ok)}/{len(rows)} targets | {len(ok)*len(peps)} docking pairs '
          f'-> {OUT_TASKS}')
    if ok:
        rmsds = [r['rmsd'] for r in ok]
        print(f'RMSD: min {min(rmsds):.2f} median {np.median(rmsds):.2f} max {max(rmsds):.2f} A')
        big = [r['target'] for r in ok if r['rmsd'] > 8]
        if big:
            print(f'WARNING: RMSD > 8 A for {big} - box transfer less reliable there')


if __name__ == '__main__':
    main()
