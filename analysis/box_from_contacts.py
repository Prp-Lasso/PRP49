"""Define the docking box from LIGAND CONTACT RESIDUES mapped onto the AF model.

For EDNRB / MDM2 / PLXNB1 the experimental complex and the AF monomer share no
rigid frame (local-fit RMSD 10-31 A), so superposition cannot transfer the box.
Instead:
  1. take the co-crystallised ligand chain from the original complex
  2. find receptor residues within CONTACT_A of any ligand atom
  3. align the experimental receptor sequence to the AF model sequence
  4. take those residues' CA coordinates *in the AF model* - no superposition
  5. box centre = their centroid; size = their extent + padding, min edge 38 A

This is local by construction: only the binding-site residues matter, so a
different global fold no longer corrupts the box.
"""
import csv
import os
import sys
from difflib import SequenceMatcher

import numpy as np

BASE = os.path.expanduser('~/LassoPep')
AF_DIR = os.path.join(BASE, 'docking_af', 'af')
SRC_DIR = os.path.join(BASE, 'docking_inputs')
OUT = os.path.join(BASE, 'docking_af')

CONTACT_A = 8.0      # ligand-receptor contact cutoff
MIN_EDGE = 38.0      # large lasso peptides need room (issue #9)
PAD = 8.0

# target -> (complex pdb, receptor chains, ligand chains)
MAPPING = {
    'EDNRB':  ('9KDF', 'R', 'A'),
    'MDM2':   ('1YCR', 'A', 'B'),
    'PLXNB1': ('7VF3', 'A', 'B'),
    'PPIA':   ('2QKI', 'A', 'G'),
    'ITGAVB3': ('1L5G', 'A', 'C'),
    'POLR2A': ('6N61', 'C', 'I'),
}

AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
       'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
       'THR':'T','TRP':'W','TYR':'Y','VAL':'V','MSE':'M','SEC':'U','PYL':'O'}


def read_pdb(path):
    """-> (receptor_seq, receptor_ca_by_resnum, ligand_atoms)"""
    seq, ca, lig = [], {}, []
    for line in open(path, errors='ignore'):
        if line.startswith('ATOM'):
            ch, resn, resi = line[21], line[17:20].strip(), line[22:27].strip()
            try:
                xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            except ValueError:
                continue
            yield_ = (ch, resn, resi, line[12:16].strip(), xyz, 'ATOM')
        elif line.startswith('HETATM'):
            ch, resn, resi = line[21], line[17:20].strip(), line[22:27].strip()
            if resn in ('HOH', 'WAT', 'NAG', 'GOL', 'EDO', 'PEG', 'SO4', 'PO4'):
                continue
            try:
                xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            except ValueError:
                continue
            yield_ = (ch, resn, resi, line[12:16].strip(), xyz, 'HETATM')
        else:
            continue
        ch, resn, resi, atom, xyz, kind = yield_
        if kind == 'ATOM':
            if resi not in ca and atom == 'CA':
                ca[resi] = (ch, xyz)
                seq.append(AA3.get(resn, 'X'))
        else:
            lig.append((ch, resn, resi, xyz))
    # rebuild sequence in residue order
    return seq, ca, lig


def parse(path, rec_chains, lig_chains):
    rec_ca, rec_seq, rec_resids, lig_xyz = {}, [], [], []
    for line in open(path, errors='ignore'):
        if line.startswith(('ATOM', 'HETATM')):
            ch = line[21]
            resn = line[17:20].strip()
            resi = line[22:27].strip()
            try:
                xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            except ValueError:
                continue
            if line.startswith('ATOM') and ch in rec_chains:
                if line[12:16].strip() == 'CA' and resi not in rec_ca:
                    rec_ca[resi] = xyz
                    rec_seq.append(AA3.get(resn, 'X'))
                    rec_resids.append(resi)
            elif ch in lig_chains and resn not in ('HOH', 'WAT'):
                lig_xyz.append(xyz)
    return rec_seq, rec_ca, rec_resids, np.array(lig_xyz) if lig_xyz else np.zeros((0, 3))


def parse_af(path):
    seq, ca = [], []
    for line in open(path, errors='ignore'):
        if line.startswith('ATOM') and line[12:16].strip() == 'CA':
            seq.append(AA3.get(line[17:20].strip(), 'X'))
            ca.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return ''.join(seq), np.array(ca)


rows = []
for target, (pdb, rec_ch, lig_ch) in MAPPING.items():
    src = os.path.join(SRC_DIR, f'{pdb}.pdb')
    af = os.path.join(AF_DIR, f'{target}.pdb')
    if not (os.path.exists(src) and os.path.exists(af)):
        print(f'{target:9s} SKIP (missing {pdb}.pdb or AF model)')
        continue
    r_seq, r_ca, r_ids, lig = parse(src, set(rec_ch), set(lig_ch))
    af_seq, af_ca = parse_af(af)
    if len(lig) == 0:
        print(f'{target:9s} no ligand atoms found in {pdb} chains {lig_ch}')
        continue
    # contact residues
    contacts = [rid for rid in r_ids
                if rid in r_ca and np.min(np.linalg.norm(lig - r_ca[rid], axis=1)) <= CONTACT_A]
    sm = SequenceMatcher(None, ''.join(r_seq), af_seq, autojunk=False)
    # map experiment residue index -> AF index
    idx_of = {rid: i for i, rid in enumerate(r_ids)}
    af_of = {}
    for blk in sm.get_matching_blocks():
        for k in range(blk.size):
            af_of[r_ids[blk.a + k]] = blk.b + k
    ca_af = [af_ca[af_of[rid]] for rid in contacts if rid in af_of]
    if len(ca_af) < 8:
        print(f'{target:9s} only {len(ca_af)} contact residues mapped - skip')
        continue
    P = np.array(ca_af)
    centre = P.mean(0)
    ext = P.max(0) - P.min(0)
    size = np.maximum(ext + PAD, MIN_EDGE)
    rows.append(dict(target=target, method=f'contact<={CONTACT_A:.0f}A', source_pdb=pdb,
                     receptor_chain=rec_ch, ligand_chain=lig_ch,
                     n_contacts=len(contacts), n_mapped=len(ca_af),
                     cx=round(float(centre[0]), 1), cy=round(float(centre[1]), 1),
                     cz=round(float(centre[2]), 1),
                     sx=round(float(size[0])), sy=round(float(size[1])), sz=round(float(size[2]))))
    print(f'{target:9s} {pdb} contacts {len(contacts):3d} mapped {len(ca_af):3d} | '
          f'box ({centre[0]:7.1f},{centre[1]:7.1f},{centre[2]:7.1f}) size '
          f'({size[0]:.0f},{size[1]:.0f},{size[2]:.0f})')

out = os.path.join(OUT, 'contact_boxes.csv')
with open(out, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else ['target'])
    w.writeheader(); w.writerows(rows)
print(f'\nwrote {out}: {len(rows)} targets')
