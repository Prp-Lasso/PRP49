"""Corrected fusion simulation.

Root cause of the false positives: the docking table names the integrin dimer
'ITGAVB3' while the candidate matrix uses 'ITGAV' and 'ITGB3' separately, so both
columns silently had NO docking score. The contested pairs (wild-type MccJ25 on
integrins) therefore competed on binding probability alone.

Fix: map ITGAVB3 -> {ITGAV, ITGB3}, then re-fuse.
"""
import numpy as np
import pandas as pd

ROOT = r'D:\deepseek_harness\prp49'
cand = pd.read_csv(f'{ROOT}\\results\\candidates_xneg.csv')
dock = pd.read_csv(f'{ROOT}\\docking\\dock_scores_long.csv')

# --- name normalisation: the dimer score applies to BOTH chains
# (only the ITGAVB3 rows may be duplicated - duplicating the whole table would
#  fabricate an ITGB3 score from unrelated targets, which is what went wrong before)
dimer = dock[dock.target == 'ITGAVB3'].copy()
dimer['target'] = 'ITGB3'
dock['target'] = dock['target'].replace({'ITGAVB3': 'ITGAV'})
dock = pd.concat([dock, dimer], ignore_index=True)
dock = dock[dock.target.isin(set(cand.target_id))]
print(f'dock rows after ITGAVB3 -> ITGAV/ITGB3 mapping: {len(dock)} '
      f'({len(dimer)} duplicated for ITGB3)')

m = cand.merge(dock[['peptide', 'target', 'dock_score']]
               .rename(columns={'peptide': 'peptide_id', 'target': 'target_id',
                                'dock_score': 'dock_fix'}),
               on=['peptide_id', 'target_id'], how='left')
print(f'dock coverage: before {cand.dock_score.notna().mean()*100:.0f}%  '
      f'after {m.dock_fix.notna().mean()*100:.0f}%')


def z(x):
    x = np.asarray(x, dtype=float)
    sd = x.std()
    return (x - x.mean()) / sd if sd > 1e-9 else np.zeros_like(x)


def compose(w_bind, w_dock, penalty, dock_col):
    sig = [(w_bind, z(m.binding_prob.fillna(m.binding_prob.mean())), m.binding_prob.notna().values)]
    if w_dock:
        sig.append((w_dock, z(-m[dock_col].fillna(m[dock_col].mean())), m[dock_col].notna().values))
    comp = np.zeros(len(m)); tw = sum(w for w, _, _ in sig)
    for w, v, a in sig:
        comp += w * np.where(a, v, penalty)
    return comp / tw


def report(comp, label):
    mm = m.copy(); mm['composite'] = comp
    mm['rank'] = mm.composite.rank(ascending=False).astype(int)
    n = len(mm)

    def rk(p, t):
        s = mm[(mm.peptide_id == p) & (mm.target_id == t)]
        return int(s.iloc[0]['rank']) if not s.empty else None

    pos = [('RES-701-3', 'EDNRB'), ('RES-701-1', 'EDNRB'), ('Anantin', 'NPR1')]
    neg = [('MccJ25', 'ITGAV'), ('MccJ25', 'ITGB3')]
    hp = [rk(*p) for p in pos]; hn = [rk(*q) for q in neg]
    pok = sum(1 for r in hp if r and r <= 0.2 * n)
    nok = sum(1 for r in hn if r and r > 0.5 * n)
    print(f'\n--- {label} ---')
    for (p, t), r in zip(pos, hp):
        print(f'  POS {p:11s} x {t:9s} rank {r:3d}/{n} {"OK " if r and r<=0.2*n else "MISS"}')
    for (p, t), r in zip(neg, hn):
        print(f'  NEG {p:11s} x {t:9s} rank {r:3d}/{n} {"OK " if r and r>0.5*n else "HIGH"}')
    print(f'  => positives {pok}/3 in top-20% | negatives controlled {nok}/2')


print('\n' + '=' * 74)
report(compose(0.6, 0.4, -0.5, 'dock_fix'), 'FIXED dock mapping + missing penalty -0.5')
report(compose(0.5, 0.5, -0.5, 'dock_fix'), '50/50 weights')
report(compose(0.4, 0.6, -0.5, 'dock_fix'), 'docking-dominant 40/60')
report(compose(0.6, 0.4, -0.5, 'dock_score'), 'old (broken) mapping, for reference')

print('\n' + '=' * 74)
print('docking scores for the contested integrin pairs (after mapping):')
for pep in ('MccJ25', 'RES-701-3', 'RES-701-1', 'Capistruin'):
    for tgt in ('ITGAV', 'ITGB3'):
        row = m[(m.peptide_id == pep) & (m.target_id == tgt)]
        if not row.empty:
            v = row.iloc[0]['dock_fix']
            print(f'  {pep:11s} x {tgt:6s} dock {"missing" if pd.isna(v) else f"{v:+.1f}"}')
