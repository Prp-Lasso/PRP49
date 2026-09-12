"""Aggregate ALL docking scores (fixed): separate pep/rec caches + score sanity filter.

Fixes over v1:
  - peptide and receptor share the SAME task name, so one cache keyed by name
    returned the peptide sequence for both -> separate pep_cache / rec_cache.
  - vina scores are filtered to a physical range; outliers are dumped to a file
    for inspection instead of polluting the training table.
"""
SERVER_SCRIPT = r'''
import csv, glob, os, statistics

AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
       'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
       'THR':'T','TRP':'W','TYR':'Y','VAL':'V','MSE':'M','SEC':'U','PYL':'O'}

def seq_from_pdb(path, max_res=6000):
    seen, out = set(), []
    try:
        for line in open(path, errors='ignore'):
            if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                key = (line[21], line[22:27])
                if key in seen:
                    continue
                seen.add(key)
                out.append(AA3.get(line[17:20].strip(), 'X'))
                if len(out) > max_res:
                    break
    except FileNotFoundError:
        return None
    return ''.join(out) or None

rows, outliers, no_seq = [], [], 0
for base, tag, taskfile in [('docking_propedia', 'b1', 'tasks_propedia.csv'),
                            ('docking_propedia_b2', 'b2', 'tasks_propedia_b2.csv')]:
    taskmap = {}
    tp = os.path.join(base, taskfile)
    if os.path.exists(tp):
        with open(tp) as tf:
            for r in csv.DictReader(tf):
                taskmap[r['id']] = (r['pep'], r['rec'])
    pep_cache, rec_cache = {}, {}
    for f in sorted(glob.glob(os.path.join(base, 'out', '*.score'))):
        tid = os.path.basename(f).split('.')[0]
        txt = open(f).read().strip()
        if not txt or tid not in taskmap:
            continue
        try:
            score = float(txt.split()[3])
        except (IndexError, ValueError):
            outliers.append((tid, 'unparsable', txt[:60]))
            continue
        if not (-200.0 <= score <= 200.0):
            outliers.append((tid, 'range', score))
            continue
        pep_name, rec_name = taskmap[tid]
        if pep_name not in pep_cache:
            pep_cache[pep_name] = seq_from_pdb(os.path.join(base, 'pep', f'{pep_name}.pdb'))
        if rec_name not in rec_cache:
            rec_cache[rec_name] = seq_from_pdb(os.path.join(base, 'rec', f'{rec_name}.pdb'))
        ps, rs = pep_cache[pep_name], rec_cache[rec_name]
        if not ps or not rs or ps == rs:
            no_seq += 1
            continue
        rows.append(dict(id=tid, pep_name=pep_name, pep_seq=ps, rec_seq=rs,
                         vina_score=score, batch=tag, n_pep=len(ps), n_rec=len(rs)))

out = 'docking_propedia/all_scores.csv'
with open(out, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['id','pep_name','pep_seq','rec_seq','vina_score','batch','n_pep','n_rec'])
    w.writeheader()
    w.writerows(rows)

with open('docking_propedia/score_outliers.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['id', 'reason', 'value'])
    w.writerows(outliers)

print(f'wrote {out}: {len(rows)} usable pairs (dropped {len(outliers)} outliers, {no_seq} bad seqs)')
for tag in ('b1', 'b2'):
    sub = [r['vina_score'] for r in rows if r['batch'] == tag]
    if sub:
        print(f'  {tag}: n={len(sub)} mean={statistics.mean(sub):.2f} median={statistics.median(sub):.2f} '
              f'min={min(sub):.2f} max={max(sub):.2f}')
allv = [r['vina_score'] for r in rows]
if allv:
    print(f'  ALL: n={len(allv)} mean={statistics.mean(allv):.2f} sd={statistics.pstdev(allv):.2f} '
          f'min={min(allv):.2f} max={max(allv):.2f}')
    pos = sum(1 for v in allv if v > 0)
    strong = sum(1 for v in allv if v <= -8)
    print(f'  positive (clash): {pos} ({100*pos/len(allv):.1f}%) | strong binders <= -8: {strong} ({100*strong/len(allv):.1f}%)')
    print(f'  peptide aa: {min(r["n_pep"] for r in rows)}-{max(r["n_pep"] for r in rows)} | '
          f'receptor aa: {min(r["n_rec"] for r in rows)}-{max(r["n_rec"] for r in rows)}')
print(f'outliers sample: {outliers[:5]}')
'''

if __name__ == '__main__':
    import sys
    sys.path.insert(0, r'D:\deepseek_harness\prp49\hpc_deploy\ops')
    from opslib import connect, run

    cli = connect()
    st, out, err = run(cli, "cd ~/LassoPep && python3 - <<'PYEOF'\n" + SERVER_SCRIPT + "\nPYEOF", timeout=2400)
    print(out)
    if err.strip():
        print('--- stderr ---')
        print(err[-1200:])
    cli.close()
