"""Prepare extended docking set: 4 homolog/variant peptides x 10 targets (tasks 82-121).

Peptide conformers:
  RES-701-1  : 9KDF chain A coordinates (1-residue different from RES-701-3; noted)
  PB1m7      : extracted from 7VF3 complex (true experimental conformer)
  Capi-var1  : 6N61 chain I coordinates (D8N variant of capistruin; noted)
  Anantin    : deferred (no experimental structure; needs LassoPred)
"""
import csv
import os

ROOT = r'D:\deepseek_harness\prp49'
DOCK = os.path.join(ROOT, 'docking')
IN = os.path.join(ROOT, 'docking_inputs')

# 1) RES-701-1 = copy of RES-701-3 coords
src = os.path.join(DOCK, 'pep', 'RES-701-3.pdb')
dst = os.path.join(DOCK, 'pep', 'RES-701-1.pdb')
open(dst, 'w').write(open(src).read())
print('RES-701-1 <- 9KDF chain A coords')

# 2) Capi-var1 = copy of Capistruin coords
open(os.path.join(DOCK, 'pep', 'Capi-var1.pdb'), 'w').write(
    open(os.path.join(DOCK, 'pep', 'Capistruin.pdb')).read())
print('Capi-var1 <- 6N61 chain I coords')

# 3) PB1m7 from 7VF3 (find peptide-sized chain)
def parse(path):
    chains = {}
    for line in open(path, errors='ignore'):
        if line.startswith('ATOM'):
            ch = line[21]
            chains.setdefault(ch, []).append(line)
    return chains

chains = parse(os.path.join(IN, '7VF3.pdb'))
sizes = {c: len(set(l[22:26] for l in rows)) for c, rows in chains.items()}
cands = {c: n for c, n in sizes.items() if 8 <= n <= 40}
if cands:
    pch = min(cands, key=cands.get)
    with open(os.path.join(DOCK, 'pep', 'PB1m7.pdb'), 'w') as f:
        for l in chains[pch]:
            f.write(l)
    print(f'PB1m7 <- 7VF3 chain {pch} ({cands[pch]} res); all chains={sizes}')
else:
    print(f'PB1m7: no peptide chain in 7VF3; chains={sizes}')

# 4) build extended tasks.csv (append 3 peptides x 10 targets)
tasks = list(csv.DictReader(open(os.path.join(DOCK, 'tasks.csv'))))
boxes = {}
for t in tasks:
    boxes[t['rec']] = (t['cx'], t['cy'], t['cz'], t['sx'], t['sy'], t['sz'])
recs = [r for r in ['EDNRB', 'PLXNB1', 'POLR2A', 'PPIA', 'FKBP1A', 'NPR1', 'C3', 'MDM2', 'ITGAVB3', 'CLPB']]
tid = 81
new = []
for pep in ['RES-701-1', 'PB1m7', 'Capi-var1']:
    for rec in recs:
        tid += 1
        cx, cy, cz, sx, sy, sz = boxes[rec]
        new.append(dict(id=tid, pep=pep, rec=rec, cx=cx, cy=cy, cz=cz, sx=sx, sy=sy, sz=sz))

with open(os.path.join(DOCK, 'tasks_ext.csv'), 'w', newline='\n') as f:
    w = csv.writer(f)
    w.writerow(['id', 'pep', 'rec', 'cx', 'cy', 'cz', 'sx', 'sy', 'sz'])
    for t in tasks:
        w.writerow([t['id'], t['pep'], t['rec'], t['cx'], t['cy'], t['cz'], t['sx'], t['sy'], t['sz']])
    for t in new:
        w.writerow([t['id'], t['pep'], t['rec'], t['cx'], t['cy'], t['cz'], t['sx'], t['sy'], t['sz']])
print(f'tasks_ext.csv: {len(tasks)} + {len(new)} = {len(tasks)+len(new)} tasks (new ids 82-{tid})')
