"""Trust-weighted docking fusion, compared against the conservative consensus.

The conservative rule (worse-of-two) failed for a diagnosable reason: the two
versions are NOT equally reliable per target. AF-vs-exp deltas range from -1.96
(PLXNB1) to +8.54 (EDNRB) kcal/mol, and for CLPB/NPR1 the "experimental" matrix
was itself built from an AlphaFold model (box rmsd 0.00), so neither version is
experimental there.

Rules compared:
  P0 exp-only        : what the shipped shortlist used
  P1 max(exp, af)    : conservative consensus (cannot be inflated by a lucky structure)
  P2 trust-weighted  : w*exp + (1-w)*af, with w set by how trustworthy the
                       "experimental" receptor really is for that target
"""
import json
import os

import numpy as np
import pandas as pd

ROOT = r'D:\deepseek_harness\prp49'
HARD_POS = [('RES-701-3', 'EDNRB'), ('RES-701-1', 'EDNRB'), ('Anantin', 'NPR1')]
HARD_NEG = [('MccJ25', 'ITGAV'), ('MccJ25', 'ITGB3')]
W_BIND, W_DOCK, MISSING = 0.40, 0.60, -0.5

# --- how much should the experimental receptor be trusted, per target? ---------
# grounded in box_provenance.csv + the AF/exp comparison:
#   0.85 - solved structure complex (or near-identical local superposition) backs the box
#   0.65 - box from contact residues of a solved complex (weaker than local fit)
#   0.50 - the "experimental" entry is itself an AlphaFold model (CLPB, NPR1),
#          so both versions are predictions and deserve equal weight
TRUST_EXP = {'EDNRB': 0.85, 'POLR2A': 0.80, 'PPIA': 0.90, 'FKBP1A': 0.90, 'C3': 0.85,
             'MDM2': 0.65, 'PLXNB1': 0.65, 'CLPB': 0.50, 'NPR1': 0.50}
DEFAULT_TRUST = 0.70

prov = pd.read_csv(os.path.join(ROOT, 'docking_af', 'box_provenance.csv'))
per_t = pd.read_csv(os.path.join(ROOT, 'docking_af', 'af_vs_exp_per_target.csv'))
af = pd.read_csv(os.path.join(ROOT, 'docking_af', 'af_scores.csv'))
ex = pd.read_csv(os.path.join(ROOT, 'docking', 'dock_scores_long.csv'))
cand = pd.read_csv(os.path.join(ROOT, 'results', 'candidates_grp.csv'))

print('=== receptor trust assignment ===')
for _, r in prov.iterrows():
    t = r.target
    print(f'  {t:9s} box={r.method:16s} quality={r.quality:12s} -> exp trust '
          f'{TRUST_EXP.get(t, DEFAULT_TRUST):.2f}')

# heterodimer split, mirroring the pipeline
dimer = ex[ex.target == 'ITGAVB3'].copy()
ex = ex[ex.target != 'ITGAVB3'].copy()
da, db = dimer.copy(), dimer.copy()
da['target'], db['target'] = 'ITGAV', 'ITGB3'
ex = pd.concat([ex, da, db], ignore_index=True)
ex = ex[ex.target != 'CTRL_9KDF']


def z(x):
    x = np.asarray(x, dtype=float)
    sd = np.nanstd(x)
    return (x - np.nanmean(x)) / sd if sd > 1e-9 else np.zeros_like(x)


rows = []
for _, c in cand.iterrows():
    pep, tgt = c.peptide_id, c.target_id
    a = af[(af.peptide == pep) & (af.target == tgt)].af_score
    e = ex[(ex.peptide == pep) & (ex.target == tgt)].dock_score
    av = float(a.iloc[0]) if len(a) else None
    ev = float(e.iloc[0]) if len(e) else None
    w = TRUST_EXP.get(tgt, DEFAULT_TRUST)
    if av is not None and ev is not None:
        worst = max(av, ev)
        wtd = w * ev + (1 - w) * av
        src = 'both'
    elif ev is not None:
        worst = wtd = ev; src = 'exp_only'
    elif av is not None:
        worst = wtd = av; src = 'af_only'
    else:
        worst = wtd = None; src = 'none'
    rows.append(dict(peptide_id=pep, target_id=tgt, af_score=av, exp_score=ev,
                     dock_worst=worst, dock_weighted=wtd, dock_source=src,
                     trust_exp=w, binding_prob=float(c.binding_prob),
                     lasso_prob=float(c.lasso_prob) if 'lasso_prob' in c else None))
d = pd.DataFrame(rows)
d['unreliable'] = d.apply(lambda r: bool(
    (r.af_score is not None and r.af_score > 0) or (r.exp_score is not None and r.exp_score > 0)), axis=1)

n = len(d)
res = {}
for label, col in [('P0_exp_only', 'exp_score'), ('P1_worst_of_two', 'dock_worst'),
                   ('P2_trust_weighted', 'dock_weighted')]:
    vals = d[col].values.astype(float)
    have = ~np.isnan(vals)
    comp = W_BIND * z(d.binding_prob.values) \
        + W_DOCK * np.where(have, z(-np.nan_to_num(vals, nan=np.nanmean(vals))), MISSING)
    d[f'composite_{label}'] = comp
    rk = pd.Series(comp).rank(ascending=False).astype(int).values
    d[f'rank_{label}'] = rk
    hp = sum(1 for p, t in HARD_POS
             if ((d.peptide_id == p) & (d.target_id == t)).any() and
             int(d[(d.peptide_id == p) & (d.target_id == t)][f'rank_{label}'].iloc[0]) <= 0.2 * n)
    hn = sum(1 for p, t in HARD_NEG
             if ((d.peptide_id == p) & (d.target_id == t)).any() and
             int(d[(d.peptide_id == p) & (d.target_id == t)][f'rank_{label}'].iloc[0]) > 0.5 * n)
    ranks = {f'{p}x{t}': int(d[(d.peptide_id == p) & (d.target_id == t)][f'rank_{label}'].iloc[0])
             for p, t in HARD_POS + HARD_NEG if ((d.peptide_id == p) & (d.target_id == t)).any()}
    res[label] = dict(positives_top20=f'{hp}/3', negatives_controlled=f'{hn}/2', ranks=ranks,
                      dock_coverage_pct=round(100 * have.mean(), 1))
    print(f'\n--- {label} ---')
    print(f'  dock coverage {100*have.mean():.0f}% | positives top-20% {hp}/3 | negatives controlled {hn}/2')
    for k, v in ranks.items():
        mark = ''
        if k in [f'{p}x{t}' for p, t in HARD_POS]:
            mark = 'OK' if v <= 0.2 * n else 'MISS'
        else:
            mark = 'OK' if v > 0.5 * n else 'MISS'
        print(f'    {k:24s} rank {v:4d}/{n} ({100*v/n:5.1f}%)  {mark}')

out = d[['peptide_id', 'target_id', 'binding_prob', 'af_score', 'exp_score', 'dock_worst',
         'dock_weighted', 'trust_exp', 'dock_source', 'unreliable',
         'rank_P0_exp_only', 'rank_P1_worst_of_two', 'rank_P2_trust_weighted']] \
    .sort_values('rank_P2_trust_weighted')
out_csv = os.path.join(ROOT, 'results', 'candidates_consensus.csv')
out.round(4).to_csv(out_csv, index=False, lineterminator='\n')
json.dump({'protocols': res, 'trust_exp': TRUST_EXP, 'fusion': {'w_binding': W_BIND, 'w_dock': W_DOCK}},
          open(os.path.join(ROOT, 'results', 'consensus_evaluation.json'), 'w'), indent=2, ensure_ascii=False)
print(f'\nwrote {out_csv}')

print('\n=== summary ===')
print(pd.DataFrame({k: {'positives': v['positives_top20'], 'negatives': v['negatives_controlled'],
                        'dock_cov': v['dock_coverage_pct']} for k, v in res.items()}).T.to_string())
print('\ntop 10 (P2 trust-weighted):')
print(out.head(10)[['rank_P2_trust_weighted', 'peptide_id', 'target_id', 'binding_prob',
                    'af_score', 'exp_score', 'dock_weighted']].round(3).to_string(index=False))
