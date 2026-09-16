"""Assemble release v1 (2026-09-16).

What changed vs v0:
  * shortlist upgraded from single-source docking to DUAL-STRUCTURE with
    per-candidate rank-stability annotations
  * route A (family label propagation) evaluated and REJECTED: family-grouped
    CV 0.7481 +- 0.0506 vs the 0.8142 baseline -> the model is unchanged
  * three fusion protocols compared and scored on the frozen truth set
  * RMSE-constrained regression partially evaluated (schemes still running)

The model itself is deliberately the same as v0: no candidate improvement
cleared the bar this round, and shipping an unvalidated model would be worse
than shipping a better-described shortlist.
"""
import json
import os
import shutil
from datetime import datetime

import pandas as pd

ROOT = r'D:\deepseek_harness\prp49'
TAG = 'v1_2026-09-16_dualstructure'
OUT = os.path.join(ROOT, 'releases', TAG)
os.makedirs(os.path.join(OUT, 'config'), exist_ok=True)

# ---- copy the deliverables ----
files = {
    'candidates.csv': os.path.join(ROOT, 'results', 'candidates_consensus.csv'),
    'consensus_evaluation.json': os.path.join(ROOT, 'results', 'consensus_evaluation.json'),
    'rank_stability.json': os.path.join(ROOT, 'results', 'rank_stability.json'),
    'routeA_cv_summary.json': os.path.join(ROOT, 'PRP49', 'runs_prop', 'cv_summary.json'),
}
for dst, src in files.items():
    if os.path.exists(src):
        shutil.copy2(src, os.path.join(OUT, dst))
        print(f'copied {dst}')
    else:
        print(f'MISSING {src}')

for cfg in ['config_improved_grouped.yaml', 'config_prop.yaml', 'config_reg_a.yaml']:
    src = os.path.join(ROOT, 'PRP49', cfg)
    if os.path.exists(src):
        shutil.copy2(src, os.path.join(OUT, 'config', cfg))

cand = pd.read_csv(os.path.join(OUT, 'candidates.csv'))
stab = json.load(open(os.path.join(OUT, 'rank_stability.json')))
cons = json.load(open(os.path.join(OUT, 'consensus_evaluation.json')))
n = len(cand)

model_ref = dict(
    tag=TAG, created=datetime.now().isoformat(timespec='seconds'),
    cluster_path='~/LassoPep/PRP49/runs_improved_grouped/checkpoints',
    fold0_md5='d0c8bcd8c6a51cab2f138ef5a4eda057',
    unchanged_from_v0=True,
    reason=('route A (family label propagation) scored 0.7481 +- 0.0506 on family-grouped CV, '
            'below the 0.8142 baseline, so it was rejected; no other candidate model cleared '
            'the bar, so v1 ships the same weights as v0 with a materially better-described shortlist'))
json.dump(model_ref, open(os.path.join(OUT, 'model_ref.json'), 'w'), indent=2, ensure_ascii=False)

notes = f"""# PRP49 Screening Release — {TAG}

Created: {datetime.now():%Y-%m-%d %H:%M} · shortlist rows: **{n}**
Weights: `{model_ref['cluster_path']}` (fold0 md5 `{model_ref['fold0_md5']}`) — **same as v0**

## 1. What changed vs v0

| aspect | v0 | v1 |
|---|---|---|
| receptor structures | single source per target | **dual structure** (experimental + AFDB v6), per-candidate stability |
| fusion rule | exp-only | three protocols compared, exp-only retained |
| shortlist columns | rank, composite | + af_score, exp_score, ranks under **3 rules**, rank_spread, stable flags |
| route A evaluated | — | ✅ evaluated and **rejected** (0.7481 < 0.8142) |

## 2. Route A (family label propagation) — rejected, honestly

| fold | 0 | 1 | 2 | 3 | 4 | **mean** |
|---|---|---|---|---|---|---|
| family-grouped AUC | 0.726 | 0.726 | 0.677 | 0.805 | 0.807 | **0.7481 ± 0.0506** |

Baseline (family-grouped, no propagation): **0.8142** → route A is **0.066 worse**.
Interpretation: propagated soft labels (352 → 759 positives, graded by family identity)
inject more noise than signal at this data scale.

> Method note: a 2-fold reading (both 0.726) suggested the gap was 0.09; the true
> 5-fold gap is 0.066 with fold-level sd 0.05. Partial-fold verdicts were wrong and
> are logged as mistake #35.

## 3. Dual-structure finding (the substantive improvement in v1)

Same (peptide, target) pair scored against an AlphaFold model vs the experimental
receptor correlates at only **Pearson +0.195**, and the top-1 peptide agrees in just
**1 of 9** targets. The worst offender is EDNRB, whose AF model sits **+8.54 kcal/mol**
above its crystal structure — exactly the target of the grade-A validation pair.

Three fusion rules, scored on the frozen truth set (3 hard positives / 2 hard negatives):

| rule | positives top-20% | negatives controlled |
|---|---|---|
| **P0 exp-only (retained)** | **3/3** | **2/2** |
| P1 worst-of-two (the "conservative" rule originally proposed) | 1/3 | 0/2 |
| P2 trust-weighted | 3/3 | 0/2 |

**The originally proposed conservative rule is harmful** and has been dropped.

## 4. Robust candidates (new deliverable)

**{stab['stable_top20']} candidates rank inside the top-20% under ALL three protocols**;
{stab['volatile']} of {n} swing by 30+ places (structure choice dominates for those).
Promote the stable set; treat volatile pairs as structure-sensitive.

Top of the stable list:
{chr(10).join('  - ' + x for x in stab['robust'][:10])}

## 5. Status of other threads

* **RMSE-constrained regression**: schemes A/B running, C re-queued after a crash caused
  by three schemes sharing one checkpoint dir (mistake #34). fold-0 numbers so far:
  A RMSE 1.584 vs baseline 2.295 (r=0.720), B 1.510 vs 1.610 (r=0.497), C 1.725 vs 2.295 (r=0.675).
  All beat their baselines; A has the best correlation. **Not yet a 5-fold conclusion.**
* **Multitask cross-domain**: in-domain 0.7429 ± 0.0444, cross-domain **0.6421 ± 0.0335**
  (single-task 0.470) — stable 5-fold gain.
* **Route C (family contrastive pretraining)**: not yet started.

## 6. Contents
`candidates.csv` (shortlist with 3-rule ranks) · `consensus_evaluation.json` ·
`rank_stability.json` · `routeA_cv_summary.json` · `model_ref.json` · `config/`
"""
open(os.path.join(OUT, 'RELEASE_NOTES.md'), 'w', encoding='utf-8').write(notes)
print(f'\nRELEASE {TAG} FROZEN -> {OUT}')
for f in sorted(os.listdir(OUT)):
    p = os.path.join(OUT, f)
    if os.path.isfile(p):
        print(f'  {f}  {os.path.getsize(p)/1024:.1f} KB')
