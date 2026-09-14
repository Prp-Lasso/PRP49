"""Final evaluation done LOCALLY: binding/lasso probabilities already exist, only the
docking table changed (Anantin merged in), so the fusion can be recomputed without a GPU.
"""
import sys
import pandas as pd
import numpy as np

sys.path.insert(0, r'D:\deepseek_harness\prp49\hpc_deploy\ops')
from opslib import connect, run

# cancel the queued GPU job (it was only needed to redo this arithmetic)
cli = connect()
st, out, err = run(cli, 'scancel 62555117 2>/dev/null; sleep 2; squeue -u $USER -o "%.10i %.14j %.8T %.10M" | head -6')
print('queue after cancel:'); print(out)
sftp = cli.open_sftp()
sftp.get('/dssg/home/acct-clswxl/clswxl-ccmbi1/LassoPep/docking/dock_scores_long.csv',
         r'D:\deepseek_harness\prp49\docking\dock_scores_long.csv')
sftp.close(); cli.close()
print('downloaded updated dock table')

ROOT = r'D:\deepseek_harness\prp49'
g = pd.read_csv(f'{ROOT}\\results\\candidates_grp.csv')
dk = pd.read_csv(f'{ROOT}\\docking\\dock_scores_long.csv')
print(f'shortlist {len(g)} | dock table {len(dk)} rows, targets: {sorted(dk.target.unique())}')

# ITGAVB3 (heterodimer) -> ITGAV / ITGB3, duplicating ONLY the dimer rows
dimer = dk[dk.target == 'ITGAVB3'].copy()
dimer['target'] = 'ITGB3'
dk.loc[dk.target == 'ITGAVB3', 'target'] = 'ITGAV'
dk = pd.concat([dk, dimer], ignore_index=True)
m = g.merge(dk.rename(columns={'peptide': 'peptide_id', 'target': 'target_id'}),
            on=['peptide_id', 'target_id'], how='left', suffixes=('', '_new'))
m['dock'] = m['dock_score_new'].fillna(m['dock_score']) if 'dock_score_new' in m else m['dock_score']
if 'dock_score_new' in m.columns:
    m['dock'] = m['dock_score_new'].combine_first(m['dock_score'])
print(f'dock coverage: {m.dock.notna().mean()*100:.0f}%')


def z(x):
    x = np.asarray(x, float)
    sd = x.std()
    return (x - x.mean()) / sd if sd > 1e-9 else np.zeros_like(x)


MISSING = -0.5
bind = z(m.binding_prob.fillna(m.binding_prob.mean()))
dock_avail = m.dock.notna().values
dockz = z(-m.dock.fillna(m.dock.mean()))
comp = (0.40 * bind + 0.60 * np.where(dock_avail, dockz, MISSING))
m['composite_final'] = comp
m['rank_final'] = pd.Series(comp).rank(ascending=False).astype(int)
n = len(m)

print(f'\n=== FROZEN TRUTH SET (final: Anantin docked, 40/60 fusion, missing penalty) ===')
hp = hn = 0
for pep, tgt, tag in [('RES-701-3', 'EDNRB', 'POS-A'), ('RES-701-1', 'EDNRB', 'POS-B'),
                      ('Anantin', 'NPR1', 'POS-B'), ('MccJ25', 'ITGAV', 'NEG'),
                      ('MccJ25', 'ITGB3', 'NEG')]:
    s = m[(m.peptide_id == pep) & (m.target_id == tgt)]
    if s.empty:
        print(f'  {tag} {pep} x {tgt}: missing'); continue
    r = int(s.iloc[0]['rank_final']); pct = 100.0 * r / n
    if tag.startswith('POS'):
        ok = pct <= 20; hp += ok
    else:
        ok = pct > 50; hn += ok
    print(f'  {tag:6s} {pep:11s} x {tgt:9s} rank {r:3d}/{n} ({pct:5.1f}%)  {"OK" if ok else "MISS"}')
print(f'  => positives top-20%: {hp}/3 | negatives controlled: {hn}/2')

an = m[m.peptide_id == 'Anantin'].sort_values('rank_final')
print('\nAnantin vs all 11 targets (final ranking):')
print(an[['rank_final', 'target_id', 'composite_final', 'binding_prob', 'dock']].head(5).to_string(index=False))

EXCLUDE = {'PB1m7': 'not a lasso peptide (RaPID macrocycle graft)',
           'Capi-var1': 'computational variant, no experimental support'}
d = m[~m.peptide_id.isin(EXCLUDE)].copy()
KNOWN = {('RES-701-3', 'EDNRB'): 'A: 9KDF cryo-EM IC50 31.5 nM',
         ('RES-701-1', 'EDNRB'): 'B: IC50 10 nM',
         ('Anantin', 'NPR1'): 'B: Kd 0.6 uM (docking from PREDICTED model)',
         ('MccJ25', 'POLR2A'): 'C: E. coli RpoC homology',
         ('Capistruin', 'POLR2A'): 'C: E. coli RpoC homology',
         ('Ubonodin', 'POLR2A'): 'C: bacterial RNAP only',
         ('Lassomycin', 'CLPB'): 'C: M. tuberculosis ClpC1 homology',
         ('MccJ25', 'ITGAV'): 'NEGATIVE: WT inactive (>10 uM)',
         ('MccJ25', 'ITGB3'): 'NEGATIVE: WT inactive (>10 uM)'}
d['evidence'] = d.apply(lambda r: KNOWN.get((r.peptide_id, r.target_id), ''), axis=1)
d['struct_source'] = d.peptide_id.map(
    lambda p: 'PREDICTED (LassoPred+tleap)' if p in ('Anantin',) else 'experimental PDB')
d = d.sort_values('composite_final', ascending=False).reset_index(drop=True)
d['final_rank'] = range(1, len(d) + 1)
out_csv = f'{ROOT}\\results\\candidates_final.csv'
d[['final_rank', 'peptide_id', 'target_id', 'composite_final', 'binding_prob', 'dock',
   'evidence', 'struct_source']].round(4).to_csv(out_csv, index=False, lineterminator='\n')
print(f'\nwrote {out_csv}: {len(d)} rows | docked {d.dock.notna().sum()} ({d.dock.notna().mean()*100:.0f}%)')
print(d.head(12)[['final_rank', 'peptide_id', 'target_id', 'composite_final', 'binding_prob',
                  'dock', 'evidence']].to_string(index=False))
