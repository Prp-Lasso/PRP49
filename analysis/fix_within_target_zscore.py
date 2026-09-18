"""Fix the cross-target comparability problem in the screening composite.

Problem: docking scores are not comparable ACROSS targets - a large/deep pocket yields
more negative scores for every peptide. The first pass used global z-scores, so DPP4 and
PREP (both very negative for all 13 peptides) swept the top of the list. That is a
property of the pocket, not of peptide specificity.

Fix: standardise BOTH signals within each target, so a peptide is judged by how it
compares with the other peptides tested against that same target - which is also exactly
what the screening question is ("which peptide for this target?").

Also emits the per-target view, which the global ranking hides entirely.
"""
import os

import numpy as np
import pandas as pd

ROOT = r'D:\deepseek_harness\prp49'
d = pd.read_csv(os.path.join(ROOT, 'results', 'candidates_screen50.csv'))
pd.set_option('display.width', 220)

KNOWN = {'POLR2A', 'CLPB', 'EDNRB', 'NPR1', 'C3', 'MDM2', 'ITGAV', 'ITGB3', 'PPIA',
         'FKBP1A', 'PLXNB1'}
d['is_known_target'] = d.target_id.isin(KNOWN)

print('=== evidence that docking is not cross-target comparable ===')
per_target = d.groupby('target_id').dock_score.agg(['mean', 'std', 'count'])
per_target = per_target.sort_values('mean')
print('  the five targets with the most negative mean docking score:')
print(per_target.head(5).round(2).to_string())
print('  the five with the least negative:')
print(per_target.tail(5).round(2).to_string())
print(f'\n  spread of target means: {per_target["mean"].min():.1f} .. {per_target["mean"].max():.1f} '
      f'(a {per_target["mean"].max()-per_target["mean"].min():.1f} kcal/mol systematic offset)')
print('  gap between the two best target means vs the spread WITHIN a target:')
within = d.groupby('target_id').dock_score.std().median()
print(f'    median within-target std {within:.2f} vs between-target mean range '
      f'{per_target["mean"].max()-per_target["mean"].min():.1f}')
print('  => the between-target offset dwarfs the within-target signal the ranking should use')


def zscore(s):
    sd = s.std()
    return (s - s.mean()) / sd if sd and sd > 1e-9 else s * 0.0


# Physically impossible docking scores (positive = steric clash / bad box) must not
# enter the standardisation: they inflate the within-target spread enormously
# (e.g. DRD2 std 44.8, ACE2 std 28.8) and would make every other peptide look good.
n_clash = int((d.dock_score > 0).sum())
print(f'\nexcluding {n_clash} clash rows (dock_score > 0) from the composite')
d['dock_ok'] = d.dock_score.notna() & (d.dock_score <= 0)
d.loc[~d.dock_ok, 'dock_score_adj'] = np.nan
d['dock_score_adj'] = d.dock_score.where(d.dock_ok)

# within-target standardisation of both signals (computed on valid dock scores only)
d['binding_z_global'] = zscore(d.binding_logit)
d['dock_z_global'] = zscore(-d.dock_score_adj)
d['binding_z_target'] = d.groupby('target_id').binding_logit.transform(zscore)
d['dock_z_target'] = (d[d.dock_ok].groupby('target_id').dock_score
                      .transform(lambda s: zscore(-s)).reindex(d.index))

# missing docking -> fall back to the binding signal alone rather than to NaN
def combine(row, bz, dz, w_b=0.40, w_d=0.60):
    if pd.isna(row[dz]):
        return row[bz]
    return w_b * row[bz] + w_d * row[dz]


d['composite_global'] = d.apply(lambda r: combine(r, 'binding_z_global', 'dock_z_global'), axis=1)
d['composite_target'] = d.apply(lambda r: combine(r, 'binding_z_target', 'dock_z_target'), axis=1)
# rank without astype(int): missing values would raise IntCastingNaNError
d['rank_global'] = d.composite_global.rank(ascending=False, na_option='bottom')
d['rank_target'] = d.composite_target.rank(ascending=False, na_option='bottom')

print('\n=== TOP 12 under GLOBAL standardisation (the first pass) ===')
print(d.sort_values('composite_global').head(12)[
    ['rank_global', 'peptide_id', 'target_id', 'composite_global', 'binding_prob',
     'dock_score']].round(3).to_string(index=False))

print('\n=== TOP 12 under WITHIN-TARGET standardisation (the correction) ===')
print(d.sort_values('composite_target').head(12)[
    ['rank_target', 'peptide_id', 'target_id', 'composite_target', 'binding_prob',
     'dock_score']].round(3).to_string(index=False))

print('\n=== how the two views differ on the positive controls ===')
for pep, tgt in [('RES-701-3', 'EDNRB'), ('RES-701-1', 'EDNRB'), ('Anantin', 'NPR1'),
                 ('MccJ25', 'POLR2A'), ('Capistruin', 'POLR2A')]:
    s = d[(d.peptide_id == pep) & (d.target_id == tgt)]
    if len(s):
        r = s.iloc[0]
        print(f'  {pep:11s} x {tgt:8s} global rank {int(r.rank_global):4d} -> '
              f'within-target rank {int(r.rank_target):4d}')

print('\n=== best peptide FOR EACH TARGET (the actual screening deliverable) ===')
best = (d.sort_values('composite_target', ascending=False)
          .groupby('target_id').head(2)
          .sort_values(['is_known_target', 'target_id']))
print(best[['target_id', 'peptide_id', 'composite_target', 'binding_prob', 'dock_score',
            'is_known_target']].round(3).to_string(index=False))

out = os.path.join(ROOT, 'results', 'candidates_screen50_within_target.csv')
d.sort_values('rank_target').to_csv(out, index=False, lineterminator='\n')
print(f'\nwrote {out}')
print(f"  {len(d)} rows, columns now include rank_global / rank_target for comparison")
