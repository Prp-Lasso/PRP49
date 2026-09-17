"""Analyse the saturation-mutagenesis and ring-size scans into concrete advice tables.

Wild-type baseline comes from candidates_grp.csv's binding_logit for the same peptide
and target, so each variant is scored as a DELTA against its own parent rather than
against some global threshold.

Structural-risk variants (substituting ring-forming Cys/Gly, etc.) are reported
separately: the model has no representation of whether the macrocycle still closes,
so a favourable score there is not actionable without structural work.
"""
import os
import sys
sys.path.insert(0, r'D:\deepseek_harness\prp49\hpc_deploy\ops')
from opslib import connect, run

import pandas as pd

ROOT = r'D:\deepseek_harness\prp49'
cli = connect()
st, h, _ = run(cli, 'echo $HOME')
h = h.strip()
sftp = cli.open_sftp()
for f in ['mutscan_matrix.csv', 'ringscan_matrix.csv']:
    sftp.get(f'{h}/LassoPep/results/{f}', os.path.join(ROOT, 'results', f))
sftp.close()
cli.close()
print('downloaded both scan matrices')

mut = pd.read_csv(os.path.join(ROOT, 'results', 'mutscan_matrix.csv'))
ring = pd.read_csv(os.path.join(ROOT, 'results', 'ringscan_matrix.csv'))
meta = pd.read_csv(os.path.join(ROOT, 'mvp_cpu', 'mut_meta.csv'))
rmeta = pd.read_csv(os.path.join(ROOT, 'mvp_cpu', 'ring_meta.csv'))
cand = pd.read_csv(os.path.join(ROOT, 'results', 'candidates_grp.csv'))

print(f'mutscan: {mut.shape} | ringscan: {ring.shape}')
print('mutscan columns:', list(mut.columns))

# wild-type logit for each (peptide, target)
wt = cand.set_index(['peptide_id', 'target_id'])['binding_logit'].to_dict()

# ---------------- 1. mutation effects ----------------
long = mut.melt(id_vars='id', var_name='target_id', value_name='logit')
long = long.merge(meta[['variant_id', 'peptide_id', 'target_id', 'position', 'wt', 'mut',
                        'flag']],
                  left_on=['id', 'target_id'], right_on=['variant_id', 'target_id'],
                  how='inner')
long['wt_logit'] = [wt.get((p, t), float('nan'))
                    for p, t in zip(long.peptide_id, long.target_id)]
long['delta'] = long.logit - long.wt_logit
long = long.dropna(subset=['delta'])

print(f'\nscored variant-target pairs: {len(long)}')
print(f'  delta > 0 (predicted improvement): {int((long.delta > 0).sum())} '
      f'({100*(long.delta > 0).mean():.1f}%)')
print(f'  delta > +1.0                     : {int((long.delta > 1.0).sum())}')
print(f'  delta < -1.0                     : {int((long.delta < -1.0).sum())}')

print('\n=== per-peptide mutation sensitivity ===')
per = long.groupby('peptide_id').agg(variants=('delta', 'size'),
                                     mean_delta=('delta', 'mean'),
                                     best_delta=('delta', 'max'),
                                     worst_delta=('delta', 'min'),
                                     n_improving=('delta', lambda s: int((s > 0.5).sum())))
print(per.round(3).sort_values('best_delta', ascending=False).to_string())

print('\n=== top 20 predicted improvements (structural-risk variants excluded) ===')
safe = long[long.flag.fillna('') == '']
top = safe.sort_values('delta', ascending=False).head(20)
print(top[['id', 'peptide_id', 'target_id', 'position', 'wt', 'mut',
           'wt_logit', 'logit', 'delta']].round(3).to_string(index=False))

print('\n=== most destabilising substitutions ===')
bot = safe.sort_values('delta').head(8)
print(bot[['id', 'peptide_id', 'position', 'wt', 'mut', 'delta']].round(3).to_string(index=False))

# per-position sensitivity for the top candidate
print('\n=== position sensitivity, top candidate ===')
for pep in long.sort_values('delta', ascending=False).peptide_id.unique()[:3]:
    sub = safe[safe.peptide_id == pep]
    if sub.empty:
        continue
    pos = sub.groupby('position').agg(n=('delta', 'size'), best=('delta', 'max'),
                                      mean=('delta', 'mean')).round(3)
    hot = pos[pos.best > 0.5].sort_values('best', ascending=False)
    print(f'  {pep}: positions where some substitution helps (best delta > 0.5):')
    print('   ', hot.head(6).to_dict('index') if len(hot) else 'none')

top.to_csv(os.path.join(ROOT, 'results', 'mutation_advice.csv'), index=False,
           lineterminator='\n')
print(f'\nwrote results/mutation_advice.csv ({len(top)} recommendations)')

# ---------------- 2. ring-size trade-off ----------------
rlong = ring.melt(id_vars='id', var_name='target_id', value_name='logit')
rlong = rlong.merge(rmeta[['variant_id', 'peptide_id', 'target_id', 'ring_wt', 'ring_new',
                           'delta', 'stability_note']].rename(columns={'delta': 'ring_delta'}),
                    left_on=['id', 'target_id'], right_on=['variant_id', 'target_id'],
                    how='inner')
rlong['wt_logit'] = [wt.get((p, t), float('nan'))
                     for p, t in zip(rlong.peptide_id, rlong.target_id)]
rlong['score_delta'] = rlong.logit - rlong.wt_logit
rlong = rlong.dropna(subset=['score_delta'])
rlong.to_csv(os.path.join(ROOT, 'results', 'ring_tradeoff.csv'), index=False,
             lineterminator='\n')
print('\n=== ring-size trade-off (binding delta vs topology rule) ===')
agg = rlong.groupby(['ring_wt', 'ring_new', 'stability_note']).agg(
    n=('score_delta', 'size'), mean_delta=('score_delta', 'mean'),
    max_delta=('score_delta', 'max')).round(3).reset_index()
print(agg.to_string(index=False))
print('\nper-peptide:')
print(rlong.groupby(['peptide_id', 'ring_new']).score_delta.mean().round(3).unstack().to_string())
print('\nwrote results/ring_tradeoff.csv')
