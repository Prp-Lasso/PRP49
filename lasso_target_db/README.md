# lasso_target_db — a target-annotated lasso peptide structure database

Built 2026-09-16 to fix a specific gap: **Propedia contains no RNA polymerase receptor**,
so route D (homology-weighted transfer from bacterial targets to human targets) could
not cover POLR2A — the target of the grade-C reference pairs MccJ25 and capistruin.

## Why the first attempt failed, and what fixed it

| step | result |
|---|---|
| full-text search `"lasso peptide"` | 98 PDB entries, but only **5** paired a lasso peptide with a larger protein, and **zero** were RNAP |
| target-driven search (`RNA polymerase microcin`) | **6N60 / 6N61 / 6N62** — three RNAP-lasso complexes the full-text search missed entirely, because their titles never contain the word "lasso" |

> Lesson: a single keyword strategy silently under-collects. Query by **target** as well
> as by **peptide class**.

## What the database contains

| PDB | peptide | target | human counterpart | confidence |
|---|---|---|---|---|
| **6N60** | MccJ25 (21 aa) + 24/29 aa | RNAP beta-prime (rpoC, 1409 aa) | **POLR2A** | high |
| **6N61** | capistruin (19/24/29 aa) | RNAP beta-prime | **POLR2A** | high |
| **6N62** | 24/29 aa lasso | RNAP beta-prime | **POLR2A** | high |
| **8IBO / 8IBP** | lassomycin (16 aa) | ClpC1 (142 aa, *M. tuberculosis*) | **CLPB** | high |
| 4CU4 | MccJ25 | FhuA receptor (706 aa) | — | low (uptake, not target) |
| 9KDF | RES-701-3 (16 aa) | EDNRB (908 aa, human) | EDNRB | experimental (grade A) |

Each complex is a full RNAP/ClpC1 assembly, so a peptide appears paired with every
subunit; the biologically relevant partner is the one named above (beta-prime for RNAP,
ClpC1 for lassomycin).

## Files
- `pdb_lasso_entries.csv` — 98 entries from the peptide-class search
- `lasso_target_pairs.csv` — all (peptide, protein) combinations extracted from them
- `rnap_lasso_entities.csv` / `rnap_lasso_pairs.csv` — the target-driven harvest
- `lasso_target_pairs_curated.csv` — **the usable pairs, with human counterparts**
- `target_functional_map.csv` — human target → counterpart, with evidence and confidence
- `search_hits.json`, `target_driven_hits.json` — raw query provenance

## What it still cannot supply
- EDNRB, PLXNB1, ITGAV, ITGB3 have **no bacterial counterpart** (metazoan-specific folds)
- PPIA/FKBP1A/NPR1/MDM2/C3 counterparts exist as **sequences** (via Pfam domains) but
  have no lasso-peptide complex structures
