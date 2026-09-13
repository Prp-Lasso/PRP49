"""Analyse the candidate shortlist for screening value."""
import pandas as pd
import numpy as np

df = pd.read_csv(r'D:\deepseek_harness\prp49\results\candidates.csv')
print(f'shortlist: {len(df)} pairs, {df.peptide_id.nunique()} peptides x {df.target_id.nunique()} targets\n')

# known positive pairs (experimentally supported in the literature)
known = [('MccJ25', 'POLR2A'), ('Capistruin', 'POLR2A'), ('RES-701-3', 'EDNRB'),
         ('Lassomycin', 'CLPB'), ('Ubonodin', 'POLR2A')]
print('=== known / literature-supported pairs ===')
for pep, tgt in known:
    row = df[(df.peptide_id == pep) & (df.target_id == tgt)]
    if row.empty:
        print(f'  {pep:12s} x {tgt:8s}  (not in matrix)')
        continue
    r = row.iloc[0]
    print(f'  {pep:12s} x {tgt:8s}  rank {int(r["rank"]):3d}/143  composite {r.composite:+.3f}  '
          f'bind {r.binding_prob:.3f}  rank_score {r.rank_score:+.2f}  '
          f'dock {r.dock_score if pd.notna(r.dock_score) else float("nan"):.1f}')

print('\n=== best target per peptide ===')
for pep, grp in df.groupby('peptide_id'):
    top = grp.nsmallest(3, 'rank')
    tops = ' | '.join(f'{r.target_id}({int(r["rank"])})' for _, r in top.iterrows())
    print(f'  {pep:14s} {tops}')

print('\n=== best peptide per target ===')
for tgt, grp in df.groupby('target_id'):
    top = grp.nsmallest(2, 'rank')
    tops = ' | '.join(f'{r.peptide_id}({int(r["rank"])})' for _, r in top.iterrows())
    print(f'  {tgt:9s} {tops}')

# how much does each signal contribute?
print('\n=== signal correlations with composite (all pairs) ===')
for c in ('binding_prob', 'rank_score', 'lasso_prob', 'dock_score'):
    sub = df[[c, 'composite']].dropna()
    if len(sub) > 5:
        rho = sub[c].corr(sub.composite, method='spearman')
        print(f'  {c:14s} rho={rho:+.3f}  (n={len(sub)})')

# ranking stability: does the top-10 survive dropping weak signals?
print('\n=== top-10 stability ===')
top10 = set(df.nsmallest(10, 'rank').peptide_id + ' x ' + df.nsmallest(10, 'rank').target_id)
alt = df.copy()
alt['alt'] = 0.55 * (alt.binding_prob - alt.binding_prob.mean()) / alt.binding_prob.std() \
             + 0.45 * alt.groupby('target_id')['rank_score'].transform(
                 lambda s: (s - s.mean()) / (s.std(ddof=0) + 1e-9))
top10_alt = set(alt.nlargest(10, 'alt').peptide_id + ' x ' + alt.nlargest(10, 'alt').target_id)
print(f'  with dock+lasso   : {len(top10)} pairs')
print(f'  model signals only: {len(top10_alt)} pairs')
print(f'  overlap: {len(top10 & top10_alt)}/10')
print('  dropped when structural/lasso terms removed:',
      ', '.join(sorted(top10 - top10_alt)) or 'none')
