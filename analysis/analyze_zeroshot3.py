"""Analyse the zero-shot probe: do the 13 scan peptides score meaningfully on
GUK1 / SSX1 / EXOSC1, targets the model has never seen?

The decisive diagnostic is whether the score varies with the TARGET. If a peptide scores
nearly the same on all three, the model is not reading the target at all - it is
emitting a peptide-level constant, which carries no binding information.
"""
import os
import sys
sys.path.insert(0, r'D:\deepseek_harness\prp49\hpc_deploy\ops')
from opslib import connect, run

import numpy as np
import pandas as pd

ROOT = r'D:\deepseek_harness\prp49'
cli = connect()
st, h, _ = run(cli, 'echo $HOME')
h = h.strip()
sftp = cli.open_sftp()
sftp.get(f'{h}/LassoPep/results/zeroshot_three_targets.csv',
         os.path.join(ROOT, 'results', 'zeroshot_three_targets.csv'))
sftp.close()
cli.close()

m = pd.read_csv(os.path.join(ROOT, 'results', 'zeroshot_three_targets.csv'))
pd.set_option('display.width', 220)
print('=== full matrix: 13 scan peptides x 3 unseen targets (logit) ===')
print(m.round(3).to_string(index=False))

targets = [c for c in m.columns if c != 'id']
num = m[targets]
print(f'\nglobal range: {num.values.min():.3f} .. {num.values.max():.3f}')

print('\n=== (b) does the score depend on the TARGET? ===')
spread = num.max(axis=1) - num.min(axis=1)
print('per-peptide spread across the three targets:')
for i, s in spread.items():
    print(f'  {m.id.iloc[i]:22s} spread {s:.3f}')
print(f'  median spread: {spread.median():.3f}   max spread: {spread.max():.3f}')

print('\n=== per-target means (are the three targets distinguishable?) ===')
print(num.mean().round(3).to_dict())
print(f'  spread of target means: {num.mean().max()-num.mean().min():.3f}')

print('\n=== (a) any extreme high scores (>8 logit)? ===')
print(f'  max value {num.values.max():.3f} | count >8: {int((num.values>8).sum())} '
      f'| count >5: {int((num.values>5).sum())}')

print('\n=== variance decomposition (peptide vs target) ===')
grand = num.values.mean()
ss_pep = len(targets) * ((num.mean(axis=1) - grand) ** 2).sum()
ss_tgt = len(num) * ((num.mean(axis=0) - grand) ** 2).sum()
ss_tot = ((num.values - grand) ** 2).sum()
print(f'  between-peptide: {100*ss_pep/ss_tot:5.1f}%')
print(f'  between-target : {100*ss_tgt/ss_tot:5.1f}%')
print(f'  residual       : {100*(ss_tot-ss_pep-ss_tgt)/ss_tot:5.1f}%')

print('\n=== (c/d) reference points ===')
p = os.path.join(ROOT, 'results', 'verify_real_complexes.csv')
if os.path.exists(p):
    v = pd.read_csv(p).set_index('id')
    print(f'  known real complex MccJ25 x RNAP (6N60): {v.loc["P1_MicrocinJ25"].filter(like="6N60").iloc[0]:.3f}')
    print(f'  MccJ25 x human POLR2A                 : {v.loc["P1_MicrocinJ25"]["HUMAN_POLR2A"]:.3f}')
p2 = os.path.join(ROOT, 'results', 'scan_matrix_ext.csv')
if os.path.exists(p2):
    s = pd.read_csv(p2)
    print(f'  known-target matrix range: {s.select_dtypes("number").values.min():.3f} .. '
          f'{s.select_dtypes("number").values.max():.3f}')

print('\n=== verdict ===')
mx = num.values.max()
med_spread = spread.median()
if mx < 2 and med_spread < 0.5:
    verdict = ('NO USABLE SIGNAL. Scores stay near zero and a given peptide scores almost '
               'identically on all three targets, i.e. the model is not reading the target '
               'at all on these proteins. Consistent with the earlier finding that absolute '
               'scores are peptide-identity dominated (62.8% of variance).')
elif med_spread < 0.5:
    verdict = ('WEAK/UNINFORMATIVE. Little target dependence; treat any ranking as noise.')
else:
    verdict = ('Target-dependent scores present - worth a closer look, but still pure '
               'extrapolation with no experimental reference pair.')
print(verdict)
m.to_csv(os.path.join(ROOT, 'results', 'zeroshot_three_targets.csv'), index=False, lineterminator='\n')
