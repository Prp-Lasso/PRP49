"""Poll 62268571 (big-box rerun) -> aggregate -> download."""
import sys
import time
sys.path.insert(0, r'D:\deepseek_harness\prp49\hpc_deploy\ops')
from opslib import connect, run

cli = connect()
job = '62268571'
print(f'[poll] watching array {job}', flush=True)
for i in range(120):
    st, states, _ = run(cli, f'squeue -j {job} -h -o "%T" 2>/dev/null | sort | uniq -c')
    _, n, _ = run(cli, 'ls ~/LassoPep/docking/out/*.score 2>/dev/null | wc -l')
    if i % 5 == 0 or not states.strip():
        print(f'[{i*40}s] {states.strip() or "(all done)"} scores={n.strip()}/81', flush=True)
    if not states.strip():
        break
    time.sleep(40)

agg = r"""
cd ~/LassoPep/docking
python3 - <<'PYEOF'
import csv, glob, os
rows = {}
for f in glob.glob('out/*.score'):
    tid = int(os.path.basename(f).split('.')[0])
    line = open(f).read().strip()
    if not line: continue
    rows[tid] = float(line.split()[3])
taskmap = {}
with open('tasks.csv') as tf:
    for r in csv.DictReader(tf):
        taskmap[int(r['id'])] = (r['pep'], r['rec'])
missing = [t for t in taskmap if t not in rows]
print('FINAL scores:', len(rows), '/', len(taskmap), 'missing:', missing)
peps, recs = [], []
for t in sorted(taskmap):
    p, r = taskmap[t]
    if p not in peps: peps.append(p)
    if r not in recs: recs.append(r)
mat = {p: {} for p in peps}
for t in sorted(taskmap):
    p, r = taskmap[t]
    mat[p][r] = rows.get(t)
with open('dock_matrix.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['peptide'] + recs)
    for p in peps:
        w.writerow([p] + [round(mat[p][r], 1) if mat[p].get(r) is not None else '' for r in recs])
print(open('dock_matrix.csv').read())
PYEOF
"""
st, out, err = run(cli, agg, timeout=120)
print(out, flush=True)
sftp = cli.open_sftp()
sftp.get('/dssg/home/acct-clswxl/clswxl-ccmbi1/LassoPep/docking/dock_matrix.csv',
         r'D:\deepseek_harness\prp49\docking\dock_matrix.csv')
print('saved docking/dock_matrix.csv', flush=True)
sftp.close()
cli.close()
print('POLL_DONE', flush=True)
