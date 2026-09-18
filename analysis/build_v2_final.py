"""D4: build release v2 (full) and pick FINAL under the corrected dual-track criteria.

Criteria note: the D4 instruction specified family-grouped CV as the primary metric, but
that was superseded on 09-17 (docs/PRP49_版本对齐.md): v0/v1 ship the grouped-by-peptide
checkpoint (0.8333 +- 0.0135), while 0.8142 belongs to a different run (runs_family,
+- 0.0477 - 3.5x noisier). Mixing them is mistake #43. So: grouped is primary, family is
a falsification check, and any quoted AUC carries its protocol.
"""
import json
import os
import shutil
from datetime import datetime

import pandas as pd

ROOT = r'D:\deepseek_harness\prp49'
REL = os.path.join(ROOT, 'releases')
TAG = 'v2_2026-09-18_full'
OUT = os.path.join(REL, TAG)
os.makedirs(os.path.join(OUT, 'config'), exist_ok=True)
os.makedirs(os.path.join(OUT, 'assets'), exist_ok=True)

# ---------------- copy deliverables ----------------
FILES = {
    'candidates.csv': 'results/candidates_consensus.csv',
    'candidates_single_structure.csv': 'results/candidates_final.csv',
    'consensus_evaluation.json': 'results/consensus_evaluation.json',
    'rank_stability.json': 'results/rank_stability.json',
    'routeC_verdict.json': 'results/routeC_verdict.json',
    'mutation_advice.csv': 'results/mutation_advice.csv',
    'ring_tradeoff.csv': 'results/ring_tradeoff.csv',
    'tail_modification_advice.csv': 'results/tail_modification_advice.csv',
    'verify_real_complexes.csv': 'results/verify_real_complexes.csv',
}
for dst, src in FILES.items():
    p = os.path.join(ROOT, src)
    if os.path.exists(p):
        shutil.copy2(p, os.path.join(OUT, dst))
        print(f'  copied {dst}')
    else:
        print(f'  MISSING {src}')

ASSETS = {
    'lasso_target_db': 'lasso_target_db',
    'homolog_domains.csv': 'mvp_cpu/homolog_domains.csv',
    'db_aux_pairs.csv': 'mvp_cpu/db_aux_pairs.csv',
    'mut_meta.csv': 'mvp_cpu/mut_meta.csv',
    'ring_meta.csv': 'mvp_cpu/ring_meta.csv',
    'contrastive_pairs.csv': 'mvp_cpu/contrastive_pairs.csv',
}
for dst, src in ASSETS.items():
    p = os.path.join(ROOT, src)
    d = os.path.join(OUT, 'assets', dst)
    if os.path.isdir(p):
        shutil.copytree(p, d, dirs_exist_ok=True)
        print(f'  copied dir {dst}/')
    elif os.path.exists(p):
        shutil.copy2(p, d)
        print(f'  copied {dst}')

for cfg in ['config_improved_grouped.yaml', 'config_contrastive.yaml', 'config_mtl.yaml']:
    p = os.path.join(ROOT, 'PRP49', cfg)
    if os.path.exists(p):
        shutil.copy2(p, os.path.join(OUT, 'config', cfg))
for doc in ['PRP49_推进总览.md', 'PRP49_版本对齐.md', 'D3_结构优化建议.md',
            'target_analysis_GUK1_SSX1_EXOSC1.md']:
    p = os.path.join(ROOT, 'docs', doc)
    if os.path.exists(p):
        shutil.copy2(p, os.path.join(OUT, doc))

# ---------------- version comparison ----------------
VERSIONS = [
    dict(tag='v0_2026-09-15_baseline', weights='runs_improved_grouped', grouped=0.8333,
         family=0.8142, hard='3/3 + 2/2', shortlist='single-structure docking',
         extra='baseline screening model'),
    dict(tag='v1_2026-09-16_dualstructure', weights='runs_improved_grouped', grouped=0.8333,
         family=0.8142, hard='3/3 + 2/2', shortlist='DUAL-structure, 3 protocols, stability flags',
         extra='route A evaluated and rejected'),
    dict(tag=TAG, weights='runs_improved_grouped', grouped=0.8333, family=0.8142,
         hard='3/3 + 2/2', shortlist='dual-structure + stability (same as v1)',
         extra='lasso_target_db (6 real complexes incl. RNAP), 3 engineering tables, '
               '4-route methodology verdicts, RNase-free verification matrix'),
]
print('\n=== version comparison (corrected dual-track criteria) ===')
print(f'{"version":30s} {"weights":24s} {"grouped":>8s} {"family":>8s} {"hard":>10s}')
for v in VERSIONS:
    print(f'{v["tag"]:30s} {v["weights"]:24s} {v["grouped"]:8.4f} {v["family"]:8.4f} {v["hard"]:>10s}')

verdict = ('All three ship the SAME checkpoint (no candidate model exceeded the grouped '
           'baseline this round), so the primary criterion does not separate them. '
           'Selection therefore falls to deliverable completeness and evidence quality, '
           'where v2 dominates: it adds the lasso-target structure database (including the '
           'RNAP complexes that give POLR2A real structural support), three engineering '
           'tables, the four-route methodology verdicts, and a real-complex readout check. '
           'FINAL = v2.')
print('\n' + verdict)

# ---------------- FINAL ----------------
FINAL = os.path.join(REL, 'FINAL_2026-09-18_v2_full')
if os.path.isdir(FINAL):
    shutil.rmtree(FINAL)
shutil.copytree(OUT, FINAL)
print(f'\ncreated {FINAL}')

notes = f"""# PRP49 FINAL release — {TAG}

Created {datetime.now():%Y-%m-%d %H:%M} · **FINAL = v2** (chosen 2026-09-18)

## 1. Why this version, and how the others differ

| version | weights | grouped CV (primary) | family CV (robustness) | hard gates | shortlist |
|---|---|---|---|---|---|
| v0 (09-15) | `runs_improved_grouped` | 0.8333 ± 0.0135 | 0.8142 ± 0.0477 | 3/3 + 2/2 | single-structure docking |
| v1 (09-16) | same | 0.8333 | 0.8142 | 3/3 + 2/2 | **dual-structure + 3-protocol ranks + stability flags** |
| **v2 (09-18)** | **same** | **0.8333** | **0.8142** | **3/3 + 2/2** | same as v1 |

**All three ship the same checkpoint.** No candidate model beat the grouped baseline this
round:

| candidate improvement | result |
|---|---|
| route A (family label propagation) | 0.7481 (family CV) — **rejected** |
| route C (family contrastive pretraining) | 0.8248 grouped, **5/5 folds below baseline** — rejected |
| route D' (domain-homolog auxiliary) | in-domain +0.061 but cross-domain −0.072 — a trade-off |
| route D v2 (real complexes added) | 0.7915 in-domain, **below D'** — not adopted |
| multitask cross-domain | 0.470 → 0.6421 — **the one clear win**, but a different task |

Since the model is identical, selection falls to **deliverable completeness** — and v2 wins:

* `lasso_target_db/` — 98 PDB entries, **6 genuine lasso/protein complexes**, incl.
  6N60/6N61/6N62 (MccJ25 & capistruin × E. coli RNAP β′ → the POLR2A counterpart)
* three engineering tables (saturation mutagenesis, ring-size trade-off, C-terminal modifications)
* `verify_real_complexes.csv` — real-complex readout check
* four-route methodology verdicts with quantitative boundaries

## 2. Primary metrics (protocol always quoted)

```
grouped-by-peptide CV : 0.8333 +- 0.0135   (runs_improved_grouped)  <- PRIMARY
family-grouped CV     : 0.8142 +- 0.0477   (runs_family)            <- robustness / falsification
stratified CV         : 0.8145             (runs_improved, LEAKY)   <- reference only
by-target CV          : 0.6430             (affinity classification, new targets)
cross-domain (MTL)    : 0.6421 vs 0.470 single-task
```

**Frozen validation gates (not relaxed):** 3/3 hard positives inside top-20%;
2/2 hard negatives controlled.

**Weights**: `~/LassoPep/PRP49/runs_improved_grouped/checkpoints`, fold0 md5 `d0c8bcd8c6a51cab2f138ef5a4eda057`.

## 3. Contained deliverables
`candidates.csv` (dual-structure shortlist, 17 candidates robust across all three fusion
protocols) · `candidates_single_structure.csv` (the earlier single-source list) ·
`mutation_advice.csv` · `ring_tradeoff.csv` · `tail_modification_advice.csv` ·
`verify_real_complexes.csv` · `routeC_verdict.json` · `assets/` · `config/` · 4 documents.

## 4. Honest limitations
1. Affinity regression RMSE 1.47–1.53 — **not usable** (threshold <1.0)
2. Absolute model scores are peptide-identity dominated (62.8% of variance); temperature
   scaling does not transfer out of domain
3. Docking is not a robust signal: changing the receptor structure changes the ranking
   (top-1 agreement 1/9 targets)
4. Engineering suggestions are **model extrapolation** — no engineered variant exists in
   the training data, and the model has no representation of whether a macrocycle closes
5. Fusion weights rest on only 5 validation points
6. EDNRB/PLXNB1/ITGAV/ITGB3 have no bacterial counterpart, so homology transfer cannot
   cover them
"""
open(os.path.join(OUT, 'RELEASE_NOTES.md'), 'w', encoding='utf-8').write(notes)
shutil.copy2(os.path.join(OUT, 'RELEASE_NOTES.md'),
             os.path.join(FINAL, 'RELEASE_NOTES.md'))

ref = json.load(open(os.path.join(OUT, 'model_ref.json'))) if os.path.exists(
    os.path.join(OUT, 'model_ref.json')) else {}
ref.update(dict(tag=TAG, created=datetime.now().isoformat(timespec='seconds'),
                cluster_path='~/LassoPep/PRP49/runs_improved_grouped/checkpoints',
                fold0_md5='d0c8bcd8c6a51cab2f138ef5a4eda057',
                unchanged_from_v0=True,
                cv_protocol='grouped by peptide',
                cv_auc=0.8333, cv_std=0.0135, family_grouped_auc=0.8142,
                reason='No candidate model exceeded the grouped baseline; v2 selected for '
                       'deliverable completeness (structure DB, engineering tables, methodology).'))
json.dump(ref, open(os.path.join(OUT, 'model_ref.json'), 'w'), indent=2, ensure_ascii=False)
shutil.copy2(os.path.join(OUT, 'model_ref.json'), os.path.join(FINAL, 'model_ref.json'))

print('\nreleases now:')
for d in sorted(os.listdir(REL)):
    p = os.path.join(REL, d)
    if os.path.isdir(p):
        n = sum(len(f) for _, _, f in os.walk(p))
        print(f'  {d:34s} {n:4d} files')
