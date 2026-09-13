"""Re-evaluate the shortlist against the CORRECTED validation set.

The curator's audit moved several pairs: MccJ25/Capistruin/Lassomycin/Ubonodin
against human POLR2A/CLPB are homology inference only (the experiments used
bacterial RNAP / mycobacterial ClpC1). Only four matrix cells have direct
human-target evidence, and one of those is an engineered graft.
"""
import pandas as pd

df = pd.read_csv(r'D:\deepseek_harness\prp49\results\candidates.csv')
n = len(df)

print('=== grade A/B: direct human-target evidence ===')
truth = [
    ('RES-701-3', 'EDNRB', 'A', '9KDF co-crystal; IC50 31.5 nM'),
    ('RES-701-1', 'EDNRB', 'B', 'IC50 10 nM radioligand binding'),
    ('Anantin', 'NPR1', 'B', 'Kd 0.6 uM (ANF receptor)'),
    ('PB1m7', 'PLXNB1', 'A', '7VF3 co-crystal (ENGINEERED graft, not natural lasso)'),
]
ranks = []
for pep, tgt, grade, note in truth:
    sub = df[(df.peptide_id == pep) & (df.target_id == tgt)]
    if sub.empty:
        print(f'  [{grade}] {pep} x {tgt}: not in matrix')
        continue
    row = sub.iloc[0]
    rk = int(row['rank'])
    ranks.append(rk)
    pct = 100.0 * rk / n
    flag = 'TOP-10%' if pct <= 10 else ('TOP-20%' if pct <= 20 else '')
    print(f'  [{grade}] {pep:11s} x {tgt:8s} rank {rk:3d}/{n} ({pct:4.1f}%) '
          f'comp {row.composite:+.3f} dock {row.dock_score:+.1f}  {flag}   {note}')

k20 = sum(1 for r in ranks if r <= 0.2 * n)
k10 = sum(1 for r in ranks if r <= 0.1 * n)
print(f'\n  direct human pairs in top-10%: {k10}/{len(ranks)} | top-20%: {k20}/{len(ranks)}')

print('\n=== grade C: cross-species homology inference (NOT validated on human) ===')
for pep, tgt, exp in [('MccJ25', 'POLR2A', 'E. coli RpoC (6N60)'),
                      ('Capistruin', 'POLR2A', 'E. coli RpoC (6N61)'),
                      ('Lassomycin', 'CLPB', 'M. tuberculosis ClpC1 (8IBO/8IBP)'),
                      ('Ubonodin', 'POLR2A', 'bacterial RNAP in vitro')]:
    sub = df[(df.peptide_id == pep) & (df.target_id == tgt)]
    if not sub.empty:
        row = sub.iloc[0]
        rk = int(row['rank'])
        print(f'  [C] {pep:11s} x {tgt:8s} rank {rk:3d}/{n} ({100.0*rk/n:4.1f}%) '
              f'comp {row.composite:+.3f}   experimental target was {exp}')

print('\n=== matrix cells with ANY support (curator audit) ===')
print('  supported cells: 8 of 143 (4 direct human, 4 homology inference)')
print('  NO support for: Capi-var1 (computational variant), Siamycin-I, Chaxapeptin,')
print('                  Sphingopyxin-I, Lariocidin; targets PPIA/FKBP1A/C3/MDM2')
print('  TRAP avoided: MccJ25 x ITGAV/ITGB3 requires the RGD-grafted variant,')
print('                not the natural sequence scored here')
