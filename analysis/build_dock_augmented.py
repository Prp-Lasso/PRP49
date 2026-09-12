"""Download aggregated docking table and build a Propedia docking-augmented pair set.

Strategy (course-learning stage S3):
  positives  : experimentally solved peptide-protein complexes, kept when the
               Vina score supports the geometry (<= -12 kcal/mol)
  negatives  : sequence-shuffled peptides against the same receptors (hard
               negatives that preserve amino-acid composition)
The docking score is carried along as an auxiliary column (`dock_score`) so it
can be used as a sample weight or a multi-task target.
"""
import os
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, r'D:\deepseek_harness\prp49\hpc_deploy\ops')
from opslib import connect

ROOT = r'D:\deepseek_harness\prp49'
cli = connect()
sftp = cli.open_sftp()
local = os.path.join(ROOT, 'docking_propedia', 'all_scores.csv')
sftp.get('/dssg/home/acct-clswxl/clswxl-ccmbi1/LassoPep/docking_propedia/all_scores.csv', local)
sftp.close()
cli.close()
print('downloaded', local)

df = pd.read_csv(local)
print(f'raw pairs: {len(df)}')

# --- peptide sanity: real peptides only (3-40 aa), no X-heavy sequences
df = df[(df.n_pep >= 3) & (df.n_pep <= 40)]
df = df[~df.pep_seq.str.contains('X')]
df = df[~df.rec_seq.str.contains('X')]
df = df[(df.n_rec >= 40) & (df.n_rec <= 4000)]
print(f'after sequence filters: {len(df)}')

# --- positives: docking geometry supports the experimental complex
pos = df[df.vina_score <= -12.0].copy()
pos['label'] = 1
print(f'positives (vina <= -12): {len(pos)}  (score mean {pos.vina_score.mean():.1f})')

# --- hard negatives: shuffled peptides against the same receptor
rng = np.random.default_rng(42)
def shuffle_seq(s):
    a = list(s)
    rng.shuffle(a)
    return ''.join(a)
neg = pos.copy()
neg['pep_seq'] = neg.pep_seq.map(shuffle_seq)
neg['label'] = 0
neg['vina_score'] = np.nan
print(f'negatives (shuffled): {len(neg)}')

out = pd.concat([pos, neg], ignore_index=True).sample(frac=1.0, random_state=42).reset_index(drop=True)
out['pep_id'] = [f'prop-{i}' for i in range(len(out))]
out = out.rename(columns={'vina_score': 'dock_score'})
cols = ['pep_id', 'pep_seq', 'rec_seq', 'label', 'dock_score', 'n_pep', 'n_rec', 'batch']
out = out[cols]
out = out.rename(columns={'rec_seq': 'prot_seq'})
out['prot_id'] = out.prot_seq.str.slice(0, 12)     # group id = receptor prefix

dst = os.path.join(ROOT, 'mvp_cpu', 'propedia_dock_pairs.csv')
out.to_csv(dst, index=False, lineterminator='\n')
print(f'wrote {dst}: {len(out)} pairs (pos {int((out.label==1).sum())} / neg {int((out.label==0).sum())})')
print(f'  receptors: {out.prot_id.nunique()} | docking score: {out.dock_score.min():.1f} .. {out.dock_score.max():.1f}')
print(out.head(3)[['pep_id', 'pep_seq', 'label', 'dock_score', 'n_pep', 'n_rec']].to_string())
