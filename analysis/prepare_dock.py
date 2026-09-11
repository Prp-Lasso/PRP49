"""Prepare docking inputs: extract peptide chains, clean receptors, define boxes.

Outputs:
  docking/pep/<name>.pdb    peptide chains (ATOM only)
  docking/rec/<name>.pdb    cleaned receptors (ATOM only, optionally trimmed)
  docking/tasks.csv         one line per dock: id,pep,rec,cx,cy,cz,sx,sy,sz,note
"""
import csv
import os
import re

ROOT = r'D:\deepseek_harness\prp49'
INDIR = os.path.join(ROOT, 'docking_inputs')
OUT = os.path.join(ROOT, 'docking')
os.makedirs(os.path.join(OUT, 'pep'), exist_ok=True)
os.makedirs(os.path.join(OUT, 'rec'), exist_ok=True)


def parse_pdb(path, first_model=True):
    """Return {chain: [(resseq, resname, x, y, z, atom_line)]} for ATOM+HETATM."""
    chains = {}
    in_model = 0
    for line in open(path, errors='ignore'):
        if line.startswith('MODEL'):
            in_model += 1
            if first_model and in_model > 1:
                break
        if line.startswith(('ENDMDL',)) and first_model:
            continue
        if line.startswith(('ATOM', 'HETATM')):
            ch = line[21]
            resseq = int(line[22:26])
            resname = line[17:20].strip()
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            chains.setdefault(ch, []).append((resseq, resname, x, y, z, line))
    return chains


def chain_stats(chains):
    return {c: len(set(r[0] for r in rows)) for c, rows in chains.items()}


def write_chain(path, rows):
    with open(path, 'w') as f:
        for _, _, _, _, _, line in rows:
            f.write(line.rstrip() + '\n')


def center(rows):
    xs = [r[2] for r in rows]
    ys = [r[3] for r in rows]
    zs = [r[4] for r in rows]
    return ((min(xs) + max(xs)) / 2, (min(ys) + max(ys)) / 2, (min(zs) + max(zs)) / 2,
            max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs))


# ---- peptides: structure id -> (peptide name, expected chain or None, trim range or None) ----
PEP_SRC = {
    'RES-701-3': ('9KDF', None, None),      # chain auto-detect (shortest)
    'MccJ25': ('1Q71', None, None),
    'Capistruin': ('6N61', 'I', None),      # from RNAP complex chain I (2KRN is NOT capistruin)
    'Siamycin-I': ('1RPB', None, None),
    'Chaxapeptin': ('2N5C', None, None),
    'Sphingopyxin-I': ('5JQF', 'A', (1, 26)),  # first 26 residues = mature core
    'Lassomycin': ('2MAI', None, None),
    'Ubonodin': ('6POR', None, None),
}

# ---- receptors ----
# name -> (pdb, kind: apo/box | complex/box-from-ligand, trim (resi range or None), note)
REC_SRC = {
    'EDNRB':  ('9KDG', 'apo', None, 'ETB apo cryo-EM'),
    'PLXNB1': ('7VF3', 'complex', None, 'PlexinB1 ECD, box from grafted lasso peptide'),
    'POLR2A': ('6EXV', 'het', None, 'Pol II, box from amanitin HETATM'),
    'PPIA':   ('1RMH', 'apo', None, 'PPIA apo X-ray'),
    'FKBP1A': ('1FKJ', 'apo', None, 'FKBP12 apo X-ray'),
    'NPR1':   ('NPR1_AF', 'trim', (1, 430), 'AF2 v6 ECD (1-430)'),
    'C3':     ('2QKI', 'complex', None, 'C3c, box from compstatin'),
    'MDM2':   ('1YCR', 'complex', None, 'MDM2, box from p53 peptide'),
    'ITGAVB3': ('1L5G', 'complex', None, 'integrin aVb3, box from RGD peptide'),
    'CLPB':   ('CLPB_AF', 'trim', (1, 300), 'AF2 v6 NTD+AAA core (1-300)'),
}
# redocking positive control
CONTROL = ('9KDF', 'complex', None, 'RES-701-3 back into ETB (redocking check)')

# ---- 1) extract peptides ----
print('=== peptides ===')
for name, (pid, want_ch, trim) in PEP_SRC.items():
    chains = parse_pdb(os.path.join(INDIR, f'{pid}.pdb'))
    stats = chain_stats(chains)
    rows = None
    if want_ch and want_ch in chains:
        rows = chains[want_ch]
        if trim:
            rows = [r for r in rows if trim[0] <= r[0] <= trim[1]]
        print(f'  {name} <- {pid} chain {want_ch} ({stats[want_ch]} res'
              + (f', trimmed {trim[0]}-{trim[1]}' if trim else '') + ')')
    else:
        cands = {c: n for c, n in stats.items() if 8 <= n <= 40}
        if cands:
            pch = min(cands, key=cands.get)
            rows = chains[pch]
            print(f'  {name} <- {pid} chain {pch} ({cands[pch]} res)')
        else:
            print(f'  {name} ({pid}): NO peptide-sized chain; chains={stats}')
            continue
    out_path = os.path.join(OUT, 'pep', f'{name}.pdb')
    write_chain(out_path, rows)

# ---- 2) receptors + boxes ----
print('=== receptors ===')
tasks = []
all_recs = dict(REC_SRC)
all_recs['CTRL_9KDF'] = CONTROL
box_pad = 6.0  # padding around ligand/receptor for box
max_box = 126.0

for name, (pid, kind, trim, note) in all_recs.items():
    path = os.path.join(INDIR, f'{pid}.pdb')
    chains = parse_pdb(path)
    stats = chain_stats(chains)
    # ligand rows for box: shortest peptide-sized chain (complex) or HETATM rows (het)
    lig_rows = []
    if kind == 'complex':
        cands = {c: n for c, n in stats.items() if 3 <= n <= 40}
        if cands:
            pch = min(cands, key=cands.get)
            lig_rows = chains[pch]
        else:
            print(f'  {name}: no ligand chain found; chains={stats}')
    elif kind == 'het':
        lig_rows = [r for ch in chains.values() for r in ch
                    if r[1] not in ('HOH', 'WAT') and r[0] <= 0]  # hetero rows have resseq<=0? use resname filter
        if not lig_rows:
            # amanitin residues in 6EXV may be numbered positive; fallback: collect HETATM from raw lines
            lig_rows = []
            for ch, rows in chains.items():
                for r in rows:
                    if r[5].startswith('HETATM'):
                        lig_rows.append(r)

    # build receptor rows: ATOM only, drop water; apply trim if given
    rec_rows = []
    for ch, rows in chains.items():
        for r in rows:
            if not r[5].startswith('ATOM'):
                continue
            if r[1] in ('HOH', 'WAT'):
                continue
            if trim and not (trim[0] <= r[0] <= trim[1]):
                continue
            rec_rows.append(r)
    if not rec_rows:
        print(f'  {name}: empty receptor after cleaning; chains={stats}')
        continue
    rec_out = os.path.join(OUT, 'rec', f'{name}.pdb')
    write_chain(rec_out, rec_rows)

    if lig_rows:
        cx, cy, cz, dx, dy, dz = center(lig_rows)
        # min edge 38 A: lasso peptides (16-28 aa) need room inside the box
        sx = min(max(dx + 2 * box_pad, 38), max_box)
        sy = min(max(dy + 2 * box_pad, 38), max_box)
        sz = min(max(dz + 2 * box_pad, 38), max_box)
    else:
        # blind box around whole receptor
        cx, cy, cz, dx, dy, dz = center(rec_rows)
        sx = min(dx + 2 * box_pad, max_box)
        sy = min(dy + 2 * box_pad, max_box)
        sz = min(dz + 2 * box_pad, max_box)
    print(f'  {name} ({pid}): {len(rec_rows)} atom lines, box {sx:.0f}x{sy:.0f}x{sz:.0f} @ ({cx:.1f},{cy:.1f},{cz:.1f}) [{note}]')
    all_recs[name] = (pid, kind, trim, note)
    tasks.append(dict(rec=name, pep=None,
                      cx=round(cx, 1), cy=round(cy, 1), cz=round(cz, 1),
                      sx=round(sx), sy=round(sy), sz=round(sz)))

# ---- 3) build full task list: 8 peptides x 10 receptors + control ----
box_by_rec = {t['rec']: t for t in tasks if t['pep'] is None}
full = []
tid = 0
for pep in PEP_SRC:
    for rec in REC_SRC:
        b = box_by_rec[rec]
        tid += 1
        full.append(dict(id=tid, pep=pep, rec=rec, cx=b['cx'], cy=b['cy'], cz=b['cz'],
                         sx=b['sx'], sy=b['sy'], sz=b['sz']))
# control: RES-701-3 -> 9KDF receptor
b = box_by_rec['CTRL_9KDF']
tid += 1
full.append(dict(id=tid, pep='RES-701-3', rec='CTRL_9KDF', cx=b['cx'], cy=b['cy'], cz=b['cz'],
                 sx=b['sx'], sy=b['sy'], sz=b['sz']))

with open(os.path.join(OUT, 'tasks.csv'), 'w', newline='\n') as f:
    w = csv.writer(f)
    w.writerow(['id', 'pep', 'rec', 'cx', 'cy', 'cz', 'sx', 'sy', 'sz'])
    for t in full:
        w.writerow([t['id'], t['pep'], t['rec'], t['cx'], t['cy'], t['cz'],
                    t['sx'], t['sy'], t['sz']])
print(f'=== tasks: {len(full)} (8 pep x 10 rec + 1 control) ===')
print('DONE')
