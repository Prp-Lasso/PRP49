"""Merge per-fold RMSE JSONs produced by the parallel array into one summary."""
import glob
import json
import os
import sys

import numpy as np

pattern = sys.argv[1] if len(sys.argv) > 1 else 'results/reg_bound_A_fold*.json'
out_path = sys.argv[2] if len(sys.argv) > 2 else 'results/reg_bound_A.json'

folds = []
for p in sorted(glob.glob(pattern)):
    try:
        d = json.load(open(p))
    except Exception as exc:
        print(f'  skip {p}: {exc}')
        continue
    fs = d.get('folds', [])
    if fs:
        folds.append(fs[0])
        print(f'  {os.path.basename(p)}: fold {fs[0]["fold"]} RMSE {fs[0]["rmse"]:.3f} '
              f'baseline {fs[0]["rmse_mean_baseline"]:.3f} r={fs[0].get("pearson")}')

folds.sort(key=lambda r: r['fold'])
if not folds:
    print('no folds found'); sys.exit(1)

ok = [r for r in folds if r['rmse'] == r['rmse']]
summary = {
    'folds': folds,
    'n_folds': len(folds),
    'rmse_mean': float(np.mean([r['rmse'] for r in ok])),
    'rmse_std': float(np.std([r['rmse'] for r in ok])),
    'baseline_rmse_mean': float(np.mean([r['rmse_mean_baseline'] for r in ok])),
    'improvement_pct': float(100 * (1 - np.mean([r['rmse'] for r in ok]) /
                                    np.mean([r['rmse_mean_baseline'] for r in ok]))),
    'pearson_mean': float(np.mean([r['pearson'] for r in ok if r.get('pearson') is not None])),
    'spearman_mean': float(np.mean([r['spearman'] for r in ok if r.get('spearman') is not None])),
    'per_target_rmse_median': float(np.median([r['rmse_per_target_median'] for r in ok
                                               if r.get('rmse_per_target_median') is not None])),
    'folds_better_than_baseline': int(sum(1 for r in ok
                                          if r['rmse'] < r['rmse_mean_baseline'])),
    'merged_from_parallel_array': True,
}
summary['verdict'] = ('USABLE (<1.0)' if summary['rmse_mean'] < 1.0 else
                      'better than baseline' if summary['rmse_mean'] < summary['baseline_rmse_mean']
                      else 'NO BETTER THAN BASELINE')
json.dump(summary, open(out_path, 'w'), indent=2)
print(f'\nmerged {len(folds)} folds -> {out_path}')
print(f"  RMSE {summary['rmse_mean']:.3f} +- {summary['rmse_std']:.3f} "
      f"vs baseline {summary['baseline_rmse_mean']:.3f} "
      f"({summary['improvement_pct']:+.1f}%)")
print(f"  folds better than baseline: {summary['folds_better_than_baseline']}/{len(folds)}")
print(f"  pearson {summary['pearson_mean']:.3f} | spearman {summary['spearman_mean']:.3f}")
print(f"  verdict: {summary['verdict']}")
