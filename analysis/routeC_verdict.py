"""Route C verdict (runs ON THE CLUSTER inside PRP49/).

Compares the 5-fold contrastive fine-tune against the family-grouped baseline.
Both runs use config with the same grouped CV and seed, so the folds correspond and a
PAIRED comparison is the right one - it removes the between-fold variance that
dominates raw means (baseline folds span 0.755-0.867).

Prints PROGRESS if folds are missing; prints the VERDICT once all five exist.
"""
import json
import os
import sys

import numpy as np

BASE_JSON = ['runs_improved_grouped/cv_summary.json', 'runs_family/cv_summary.json']
# fallback: the family-grouped baseline recorded in the handover document
BASE_FALLBACK = [0.852, 0.867, 0.777, 0.819, 0.755]

folds, missing = {}, []
for i in range(5):
    p = f'runs_ctrft_fold{i}/cv_summary.json'
    if os.path.exists(p):
        d = json.load(open(p))
        folds[i] = d['folds'][0]['auc'] if d.get('folds') else d.get('auc_mean')
    else:
        missing.append(i)

print(f'route C folds available: {sorted(folds)} | missing: {missing}')

base, base_src = None, 'fallback constants'
for p in BASE_JSON:
    if os.path.exists(p):
        d = json.load(open(p))
        if d.get('folds'):
            base = [x['auc'] for x in d['folds']]
            base_src = p
            break
if base is None:
    base = BASE_FALLBACK
print(f'baseline: {[round(x, 4) for x in base]}  (source: {base_src})')

if missing:
    print('\nPROGRESS ONLY - folds still running:')
    for i in sorted(folds):
        print(f'  fold {i}: AUC {folds[i]:.4f}  (baseline {base[i]:.4f}, '
              f'delta {folds[i]-base[i]:+.4f})')
    done = [folds[i] - base[i] for i in sorted(folds)]
    print(f'  mean delta over {len(done)} folds: {np.mean(done):+.4f}')
    sys.exit(0)

cur = [folds[i] for i in range(5)]
deltas = [c - b for c, b in zip(cur, base)]
print('\n=== per-fold (paired) ===')
print(f'{"fold":>5} {"baseline":>9} {"routeC":>9} {"delta":>8}')
for i in range(5):
    print(f'{i:>5} {base[i]:>9.4f} {cur[i]:>9.4f} {deltas[i]:>+8.4f}')
print(f'{"mean":>5} {np.mean(base):>9.4f} {np.mean(cur):>9.4f} {np.mean(deltas):>+8.4f}')
print(f'{"std":>5} {np.std(base):>9.4f} {np.std(cur):>9.4f}')

n = len(deltas)
se = np.std(deltas, ddof=1) / np.sqrt(n) if n > 1 else float('nan')
t = np.mean(deltas) / se if se and se == se and se > 0 else float('nan')
wins = sum(1 for d in deltas if d > 0)
print(f'\npaired mean delta {np.mean(deltas):+.4f} | sd {np.std(deltas, ddof=1):.4f} | '
      f'se {se:.4f} | t = {t:+.2f} (df={n-1})')
print(f'folds improved: {wins}/5')

if np.mean(deltas) > 0 and wins >= 4:
    verdict = 'ROUTE C HELPS (consistent gain across folds)'
elif np.mean(deltas) > 0:
    verdict = 'ROUTE C MARGINAL (positive mean but inconsistent across folds)'
elif abs(np.mean(deltas)) < 0.01:
    verdict = 'ROUTE C NEUTRAL (no material difference)'
else:
    verdict = 'ROUTE C HURTS'
print(f'VERDICT: {verdict}')

out = dict(baseline=base, routeC=cur, deltas=deltas,
           mean_delta=float(np.mean(deltas)), paired_t=float(t),
           folds_improved=int(wins), verdict=verdict, baseline_source=base_src)
json.dump(out, open('../results/routeC_verdict.json', 'w'), indent=2)
print('wrote ../results/routeC_verdict.json')
print('ROUTEC_VERDICT_DONE')
