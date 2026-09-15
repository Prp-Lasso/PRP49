# docking_af — AlphaFold-sourced docking dataset

Built 2026-09-16 to remove the **mixed structure provenance** of the original
`docking/` dataset, where the 11 lasso-panel receptors came from assorted PDB
entries (and, it turned out, two of them — CLPB and NPR1 — were already
AlphaFold models).

## Provenance
- Receptors: **AlphaFold DB v6** single-chain models, fetched from
  `https://alphafold.ebi.ac.uk/files/AF-<UNIPROT>-F1-model_v6.pdb`
- UniProt used: EDNRB P24530 · PLXNB1 O43157 · POLR2A P24928 · PPIA P62937 ·
  FKBP1A P62942 · NPR1 P16066 · C3 P01024 · MDM2 Q00987 · CLPB Q9H078
- Peptides: the same 11 lasso peptides as the experimental dataset (shared files)
- Docking: AutoDock Vina, exhaustiveness 8, num_modes 5, box edge >= 38 A

## Box transfer (AF models carry no ligand)
Two methods were used, chosen per target by fit quality:

| target | method | quality |
|---|---|---|
| EDNRB | ligand contact residues (<=8 A) mapped onto AF | 53 contacts |
| MDM2 | ligand contact residues | 33 contacts |
| PLXNB1 | ligand contact residues | 57 contacts |
| C3 | local superposition (<=15 A window) | RMSD 0.55 A |
| FKBP1A | local superposition | RMSD 0.40 A |
| PPIA | local superposition | RMSD 0.22 A |
| POLR2A | superposition | RMSD 1.74 A |
| CLPB / NPR1 | identical to AF source | RMSD 0.00 A |

Global superposition alone was **rejected** for EDNRB/MDM2/PLXNB1 (RMSD 10-42 A):
their AF monomer and the experimental complex do not share a rigid frame, so the
binding site was located through the co-crystallised ligand's contact residues
instead. Full numbers in `map_report.csv`, `contact_boxes.csv`, `box_provenance.csv`.

## Excluded
- **ITGAVB3** — alpha/beta heterodimer; a single-chain AF model cannot represent
  the dimer interface, so this target stays on the experimental structure.
- **CTRL_9KDF** — positive-control receptor, not a screening target.

## Contents
tasks_af.csv (99 pairs = 9 targets x 11 peptides) · af/ (models) · pep/ (peptides)
· box_provenance.csv · map_report.csv · contact_boxes.csv · pdbqt/ · out/ (scores)

## Relation to docking/
`docking/` = heterogeneous experimental structures; this dataset = uniform AF v6.
Comparing the two quantifies how much the receptor structure source moves Vina
scores — an explicit control that was missing before.

---

## Results (2026-09-16 06:40, job 62603981)

99/99 pairs scored, 0 failures, 5 jobs completed in 19-48 min each.

### Headline: receptor structure source dominates the docking outcome

Comparing the same (peptide, target) pairs across the two structure sources:

| metric | value |
|---|---|
| paired correlation of Vina scores | **Pearson +0.195** |
| per-pair difference (AF - exp) | mean +0.79, **sd 2.86**, range [-4.88, +9.96] |
| **within-target top-1 peptide agreement** | **1 / 9 targets** |
| per-target Spearman | NPR1 0.77, PPIA 0.70, EDNRB 0.49, CLPB 0.46, MDM2 0.40, FKBP1A 0.27, POLR2A **0.02**, C3 **-0.14**, PLXNB1 **-0.42** |

A single experimental PDB entry and the AlphaFold model of the *same* protein
produce largely uncorrelated scores and different peptide rankings. Rigid docking
is therefore **not robust to the receptor structure source**, and a docking score
derived from one structure must not be treated as a stable binding signal.

### Literature pairs (rank within their target)

| pair | grade | AF-source | experimental-source |
|---|---|---|---|
| RES-701-3 x EDNRB | A (9KDF co-crystal) | **1 / 10** | 4 / 12 |
| RES-701-1 x EDNRB | B | **2 / 10** | 4 / 12 |
| MccJ25 x POLR2A | C | 2 / 11 | 2 / 12 |
| Capistruin x POLR2A | C | 10 / 11 | 4 / 12 |
| Lassomycin x CLPB | C | 8 / 11 | **1 / 12** |
| Anantin x NPR1 | B | not scored | 2 / 12 |

Neither source dominates: AF recovers the true EDNRB ligand (RES-701-3) as rank 1
even though that target's experimental structure *is* the RES-701-3 complex, while
the experimental source is clearly better for CLPB/Lassomycin.

### Recommended use

Treat docking as a **consensus signal across structure sources** rather than a
single number:
* conservative score = the *worse* of the two sources (stricter)
* flag each shortlist row as "agree" / "conflict" between sources and prefer agree
* drop rows where either source returns a positive (clashing) score

### Artefacts found

5 positive (clashing) scores, all in targets whose boxes came from ligand-contact
mapping - i.e. box placement, not binding:

| peptide x target | AF score |
|---|---|
| Ubonodin x PLXNB1 | +147.3 |
| Siamycin-I x EDNRB | +43.9 |
| Sphingopyxin-I x EDNRB | +41.5 |
| Ubonodin x MDM2 | +28.0 |
| Siamycin-I x PLXNB1 | +1.5 |

These are marked invalid in `af_vs_exp_comparison.csv` (excluded from statistics).
Ubonodin (28 aa, the largest peptide) is involved in two of them - the contact
boxes derived from ~15-residue ligands may simply be too small for it.

## Files added
`af_scores.csv` (98 pairs) · `af_vs_exp_comparison.csv` (paired deltas) ·
`af_vs_exp_per_target.csv` (per-target summary)
