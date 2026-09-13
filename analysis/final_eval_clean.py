"""Final shortlist evaluation against the curator's frozen 25-pair audit.

Acceptance cells (direct human-target evidence, inside the 13x11 matrix):
  RES-701-3 x EDNRB (A) | RES-701-1 x EDNRB (B) | Anantin x NPR1 (B)
Hard negative: wild-type MccJ25 on ITGAV/ITGB3 (>10,000 nM, inactive)
Reference only (grade C, cross-species homology): MccJ25/Capistruin/Ubonodin x POLR2A,
  Lassomycin x CLPB
"""
import pandas as pd

df = pd.read_csv(r'D:\deepseek_harness\prp49\results\candidates.csv')
n = len(df)
print(f'shortlist rows: {n} | columns: {list(df.columns)[:8]}')
print(f'dock coverage: {df.dock_score.notna().mean()*100:.0f}%\n')


def show(pep, tgt):
    sub = df[(df.peptide_id == pep) & (df.target_id == tgt)]
    if sub.empty:
        return None
    row = sub.iloc[0]
    return dict(rank=int(row['rank']), pct=100.0 * int(row['rank']) / n,
                comp=float(row['composite']), bind=float(row['binding_prob']),
                dock=float(row['dock_score']) if pd.notna(row['dock_score']) else float('nan'))


print('=== ACCEPTANCE: hard positives (direct human-target evidence) ===')
hard = [('RES-701-3', 'EDNRB', 'A', '9KDF cryo-EM, IC50 31.5 nM'),
        ('RES-701-1', 'EDNRB', 'B', 'IC50 10 nM radioligand'),
        ('Anantin', 'NPR1', 'B', 'Kd 0.6 uM')]
h10 = h20 = 0
for pep, tgt, g, note in hard:
    s = show(pep, tgt)
    if not s:
        print(f'  [{g}] {pep} x {tgt}: NOT IN MATRIX')
        continue
    h10 += s['pct'] <= 10
    h20 += s['pct'] <= 20
    flag = 'TOP-10%' if s['pct'] <= 10 else ('TOP-20%' if s['pct'] <= 20 else '')
    print(f"  [{g}] {pep:11s} x {tgt:9s} rank {s['rank']:3d}/{n} ({s['pct']:5.1f}%)  "
          f"comp {s['comp']:+.3f}  bind {s['bind']:.3f}  dock {s['dock']:+.1f}  {flag}  ({note})")
print(f'\n  >>> top-10% hit rate {h10}/3 ({100*h10//3 if h10 else 0}%) | '
      f'top-20% hit rate {h20}/3   [acceptance threshold: >=70% in top-20%]')

print('\n=== HARD NEGATIVE (wild-type MccJ25 inactive on integrins, >10 uM) ===')
for tgt in ('ITGAV', 'ITGB3'):
    s = show('MccJ25', tgt)
    if s:
        verdict = 'OK (ranked in bottom half)' if s['pct'] > 50 else 'FAIL (ranked high)'
        print(f"  MccJ25 x {tgt:9s} rank {s['rank']:3d}/{n} ({s['pct']:5.1f}%)  "
              f"comp {s['comp']:+.3f}  {verdict}")

print('\n=== REFERENCE ONLY: grade-C homology cells (bacterial targets inferred) ===')
for pep, tgt in [('MccJ25', 'POLR2A'), ('Capistruin', 'POLR2A'),
                 ('Ubonodin', 'POLR2A'), ('Lassomycin', 'CLPB')]:
    s = show(pep, tgt)
    if s:
        print(f"  {pep:11s} x {tgt:9s} rank {s['rank']:3d}/{n} ({s['pct']:5.1f}%)  comp {s['comp']:+.3f}")

print('\n=== final top 12 ===')
print(df.head(12)[['rank', 'peptide_id', 'target_id', 'composite', 'binding_prob', 'dock_score']]
      .to_string(index=False))
