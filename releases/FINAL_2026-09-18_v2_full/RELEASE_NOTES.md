# PRP49 FINAL release — v2_2026-09-18_full

Created 2026-09-18 08:00 · **FINAL = v2** (chosen 2026-09-18)

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
