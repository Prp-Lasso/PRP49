"""Merge both box-transfer methods into the final AF docking task list.

Method preference per target:
  * contact-residue boxes  - used where superposition fails (EDNRB/MDM2/PLXNB1)
  * local/global fit boxes - used where the fit is tight (RMSD < 2 A)
Each row records which method produced it and its quality metric, so the dataset
documents its own reliability.

ITGAVB3 is excluded: it is an alpha/beta heterodimer, and an AFDB single-chain
model cannot represent the dimer interface. It stays on the experimental version.
"""
import csv
import os

BASE = os.path.expanduser('~/LassoPep')
AF = os.path.join(BASE, 'docking_af')
TASKS = os.path.join(BASE, 'docking', 'tasks.csv')

GOOD_FIT_RMSD = 2.0
EXCLUDE = {'ITGAVB3', 'CTRL_9KDF'}

# --- collect candidates ---
cand = {}
# 1) contact-residue boxes (preferred where available)
cpath = os.path.join(AF, 'contact_boxes.csv')
if os.path.exists(cpath):
    for r in csv.DictReader(open(cpath)):
        cand[r['target']] = dict(target=r['target'], method=r['method'], quality='contacts',
                                 cx=float(r['cx']), cy=float(r['cy']), cz=float(r['cz']),
                                 sx=float(r['sx']), sy=float(r['sy']), sz=float(r['sz']))
# 2) superposition boxes (only if the local fit was tight)
mpath = os.path.join(AF, 'map_report.csv')
if os.path.exists(mpath):
    for r in csv.DictReader(open(mpath)):
        t = r.get('target')
        if not t or r.get('status') != 'ok' or t in cand:
            continue
        rmsd = float(r.get('rmsd') or 999)
        if rmsd > GOOD_FIT_RMSD:
            continue
        cand[t] = dict(target=t, method=r.get('fit', 'fit'), quality=f'rmsd {rmsd:.2f}A',
                       cx=float(r['cx']), cy=float(r['cy']), cz=float(r['cz']),
                       sx=float(r['sx']), sy=float(r['sy']), sz=float(r['sz']))

final = {t: v for t, v in cand.items() if t not in EXCLUDE}
dropped = {t: v for t, v in cand.items() if t in EXCLUDE}

print('=== final AF dataset targets ===')
for t, v in sorted(final.items()):
    print(f'  {t:9s} {v["method"]:18s} {v["quality"]:12s} box ({v["cx"]:7.1f},{v["cy"]:7.1f},{v["cz"]:7.1f})')
print(f'excluded: {sorted(dropped)} (dimer / control)')

peps = sorted({r['pep'] for r in csv.DictReader(open(TASKS))})
out = os.path.join(AF, 'tasks_af.csv')
with open(out, 'w', newline='') as f:
    w = csv.writer(f, lineterminator='\n')   # Windows default would emit CRLF (#24)
    w.writerow(['id', 'pep', 'rec', 'cx', 'cy', 'cz', 'sx', 'sy', 'sz'])
    tid = 7000
    for t, v in sorted(final.items()):
        for p in peps:
            w.writerow([tid, p, f'af_{t}', v['cx'], v['cy'], v['cz'],
                        round(v['sx']), round(v['sy']), round(v['sz'])])
            tid += 1

with open(os.path.join(AF, 'box_provenance.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, lineterminator='\n',
                       fieldnames=['target', 'method', 'quality', 'cx', 'cy', 'cz', 'sx', 'sy', 'sz'])
    w.writeheader()
    for t, v in sorted(final.items()):
        w.writerow(v)

print(f'\nwrote {out}: {len(final)} targets x {len(peps)} peptides = {len(final)*len(peps)} pairs')
