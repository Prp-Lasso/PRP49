# Lasso peptide × target protein — curated experimentally supported pairs

Companion note to `lasso_known_pairs.csv` (25 rows). Compiled for use as a validation set
for a lasso-peptide / human-protein interaction scorer.

Verification basis: RCSB PDB REST + search APIs (entry titles, polymer entities, UniProt
cross-references, and partner chains recomputed directly from deposited coordinates),
UniProt REST, PubMed/Europe PMC abstracts and full texts, the local LassoPred database
(`lassopred_database.csv`, 4,749 entries), and independent verification threads. Every row
in the CSV carries a DOI or PMID.

---

## 1. Counts per evidence grade

| Grade | Meaning | Rows |
|---|---|---|
| **A** | Co-crystal / cryo-EM structure **of the complex** exists | **6** |
| **B** | Direct biochemical binding or inhibition assay (radioligand, ELISA, in vitro transcription, autophosphorylation, enzyme inhibition) | **13** |
| **C** | Genetic / functional / homology-based inference only | **6** |
| | **Total** | **25** |

Grade A (6): MccJ25×RpoC (6N60), MccJ25×FhuA (4CU4), Capistruin×RpoC (6N61),
RES-701-3×EDNRB (9KDF), Lassomycin×ClpC1-NTD (8IBO/8IBP), Lariocidin×30S/16S rRNA
(9DFC/9DFD/9DFE).

Grade C (6): MccJ25×POLR2A, Capistruin×POLR2A, Lassomycin×CLPB, Ubonodin×POLR2A
(homology inferences onto the *human* target), Siamycin-I×HIV Env, RP-71955×HIV Env.
**No grade C row was promoted to B.** Where a bacterial complex exists but the human
ortholog was never tested, the row stays C while the bacterial counterpart is a separate
grade A/B row. The siamycin×HIV row is a deliberately *strong* grade C — a real cell-fusion
ED50 with resistance mapping, but no binding assay and no complex.

### Corrections to earlier notes in this workspace
`docs/lasso_known_targets.md` needs updating:
* It lists **6N62** as the MccJ25–RNAP complex. Verified: **6N60** is MccJ25, **6N61** is
  capistruin, **6N62** is RNAP alone (the apo reference in the same paper).
* It states lassomycin has no complex structure. **8IBO / 8IBP are MtClpC1-NTD–lassomycin
  co-crystals** — lassomycin is grade A.
* It omits lariocidin (Nature 2025), ubonodin/RNAP, FhuA–MccJ25, BI-32169×GCGR, and the
  engineered MccJ25-RGD–integrin pairs.
* Anantin's sequence is `GFIGWGNDIFGHYSGDF`, not `GFIGW...` (Gly6/Asn7).
* Siamycin I and RP 71955 are **different peptides** and must not be merged.
* Its chaxapeptin and sphingopyxin I entries overstate the evidence (see §4).

### Coordinate-level checks made in this pass
Partner chains and binding subunits were recomputed from deposited coordinates rather than
trusted from titles:
* **6N60** (MccJ25, chain M): 62 peptide–RNAP atom pairs within 4 Å to **β′ / RpoC** (min
  2.90 Å) and 22 to β / RpoB (min 2.46 Å); **zero** contacts to α, ω, σ70 or DNA.
* **6N61** (capistruin, chain I): 78 pairs within 4 Å to **β′ / RpoC** (min 2.46 Å) and only
  1 to β / RpoB — essentially exclusively RpoC.
* **9KDF**: partner chains are the calcineurin-fusion ETB and FKBP1B; FKBP1B is a fusion
  chaperone, not a target (see §2).
* **8IBO / 8IBP**: partner is the MtClpC1 N-terminal domain (P9WPC9) — a real complex.
* **McjD entries 4PL0 / 5OFP / 5OFR / 5EG1 / 8PX9** are all **apo** — there is no
  McjD–MccJ25 co-crystal, so MccJ25×McjD is not a grade A pair.
* **1RPB / 1RPC** (RP 71955): disulfides Cys1–Cys13 and Cys7–Cys19 plus isopeptide
  Cys1 N → Asp9 CG, i.e. a 9-residue ring — confirmed in the deposited mmCIF.

---

## 2. Human targets in the set

Directly measured human targets:

| Target | UniProt | Peptide | Grade |
|---|---|---|---|
| Endothelin receptor type B (EDNRB) | P24530 | RES-701-3 | A |
| Endothelin receptor type B (EDNRB) | P24530 | RES-701-1 | B |
| Atrial natriuretic peptide receptor 1 (NPR1 / GC-A) | P16066 | Anantin | B |
| Glucagon receptor (GCGR) | P47871 | BI-32169 | B |
| Integrin αvβ3 (ITGAV / ITGB3) | P06756;P05106 | MccJ25(RGD) 12 nM, MccJ25(RGDF) 4.1 nM — **engineered grafts** | B |
| Integrin αvβ6 / αvβ8 (ITGAV / ITGB6 / ITGB8) | P06756;P18564;P26012 | Lassotide-44/45/46 — **engineered lasso variants** | B |

The lassotides (J Am Chem Soc 2025;147:32522, PMID 40891737) are MccJ25-derived non-natural
lasso peptides with **IC50 1–2 nM** against αvβ6/αvβ8 measured by LAP-competition
ELISA/AlphaLISA — the strongest modern example of an engineered lasso peptide with a
validated human target, and the only peer-reviewed one.

Homology-inferred human targets (grade C, never measured): POLR2A (P24928) from *E. coli*
RpoC for MccJ25 / Capistruin / Ubonodin; CLPB (Q9H078) from *M. tuberculosis* ClpC1 for
Lassomycin.

**Not a target:** the FKBP1B chain (P68106) in cryo-EM entry **9KDF** is a
calcineurin/FK506 fusion chaperone used to solve the GPCR structure. It is a construct
artefact and must **not** become an FKBP1A/FKBP1B positive label.

Only **four natural lasso peptides** have a defined human protein target at all
(RES-701-1, RES-701-3, anantin, BI-32169). **None** of them was found by pull-down,
chemical proteomics or phage display — all four came from classical pharmacology.
Proprietary display-derived "lassotides" exist but are not peer-reviewed.

---

## 3. Usability in the 13-peptide × 11-target scoring matrix

Matrix peptides: RES-701-3, RES-701-1, MccJ25, Capistruin, Capi-var1, Siamycin-I,
Chaxapeptin, Sphingopyxin-I, Lassomycin, Ubonodin, Anantin, PB1m7, Lariocidin.
Matrix targets: EDNRB, PLXNB1, POLR2A, PPIA, FKBP1A, NPR1, C3, MDM2, ITGAV, ITGB3, CLPB.

**7 of the 143 cells have any literature support; only 3 are direct (grade A/B).**

| Peptide | Target | Grade | Basis |
|---|---|---|---|
| RES-701-3 | EDNRB | **A** | 9KDF cryo-EM complex |
| RES-701-1 | EDNRB | **B** | IC50 10 nM, radioligand binding |
| Anantin | NPR1 | **B** | Kd 0.6 µM |
| MccJ25 | POLR2A | C | homology from *E. coli* RpoC |
| Capistruin | POLR2A | C | homology from *E. coli* RNAP |
| Lassomycin | CLPB | C | homology from Mtb ClpC1 |
| Ubonodin | POLR2A | C | homology; weakest |

### Rows with no usable positive
* **PB1m7** — **removed from the lasso set entirely; it is not a lasso peptide.** See §4.
* **Capi-var1** — `GTPGFQTPDNRVISRFGFN` is a *computational* D8N variant of capistruin
  (LassoPred LP_318), used in this workspace as a 6N61-derived conformer. **Zero**
  experimental evidence. Exclude from validation; its matrix row cannot supply positives.
* **Siamycin-I** — verified targets are lipid II, FsrC/VanS and MLCK; none is in the matrix.
* **Chaxapeptin** — no molecular target (see §4).
* **Sphingopyxin-I** — no molecular target **and no reported bioactivity at all** (see §4).
* **Lariocidin** — targets 16S rRNA / the 30S subunit, not a matrix target.
* **PPIA, FKBP1A, C3, MDM2, PLXNB1** — no validated lasso peptide partner exists. Legitimate
  *negative* candidates, but absence of evidence is not evidence of absence.

### Traps worth guarding
1. **ITGAV/ITGB3.** The only experimental lasso-scaffold binders of αvβ3 are the
   **engineered grafts** `GGAGHVPEYFVRGDTPISFYG` and `GGAGHVPEYFVRGDFPISFYG`. In the same
   ELISA, **wild-type MccJ25 gives >10,000 nM** — the natural peptide is inactive, so a
   model rewarded for scoring natural MccJ25 high on ITGAV/ITGB3 is being taught a wrong
   label. Note further that **MccJ25(RGD) is potent but NON-selective**: the inventors'
   patent WO2012146729A1 (Table 7) gives αvβ3 17 ± 9, αvβ5 170 ± 37, α5β1 855 ± 191 and
   platelet **αIIbβ3 29.7 ± 2.9 nM** — an αIIbβ3/αvβ3 ratio of only ~1.7, with the patent
   itself noting the "absence of selectivity for the platelet receptor". If the model
   treats each peptide as having one target, label MccJ25(RGD) as αvβ3 but flag it as
   multi-integrin/promiscuous, and use **MccJ25(RGDF)** (4.1 nM αvβ3, 725 nM αIIbβ3, a
   177-fold window) as the clean αvβ3-selective entry.
2. **POLR2A rows.** All are grade C. Down-weight them against the EDNRB/NPR1 cells and do
   not let them be the only thing separating model weights.
3. **PLXNB1 is now empty.** It lost its only matrix positive when PB1m7 was excluded — do
   not re-add it without a genuine lasso peptide binder.

### Recommended use
* **Hard validation cells (3):** RES-701-3×EDNRB, RES-701-1×EDNRB, Anantin×NPR1. These are
  human GPCR pairs with genuinely different pharmacology (inverse agonist vs antagonist vs
  enzyme-coupled receptor), which is a useful spread.
* **Soft positives (4):** the grade C homology cells — ranking sanity checks only.
* **Off-matrix positives (17 rows in the CSV):** MccJ25×RpoC, MccJ25×FhuA,
  Capistruin×RpoC, Klebsidin×RNAP, Acinetodin×RNAP, Ubonodin×RNAP, Lassomycin×ClpC1 (×2),
  Lariocidin×30S, BI-32169×GCGR, Siamycin×lipid II, Siamycin×FsrC/VanS, Siamycin×MLCK,
  MccJ25(RGD/RGDF)×αvβ3, Siamycin×HIV Env, RP-71955×HIV Env. Enough distinct positives to
  fit a handful of weights without touching the 3 in-matrix validation cells.
* **Hard negatives available:** wild-type MccJ25 vs αvβ3/αvβ5 (>10 µM, same ELISA as the
  17.2 nM graft); RES-701-1 at rat ETB (~100× weaker, IC50 1.2 µM, than human ETB 10 nM).

---

## 4. Excluded and unverified — **do not use as validation data**

### Excluded because they are not lasso peptides
* **PB1m7 and PB1m6A9 (PDB 7VF3, 7VG7).** The PDB titles read "lasso-grafted", but
  *lasso-grafting* is a protein-engineering term from Mihara et al., Nat Commun
  2021;12:1543, doi:10.1038/s41467-021-21875-0 (PMID 33750839) — implanting **RaPID
  (mRNA-display) thioether-macrocyclic peptides** into surface-exposed loops of carrier
  proteins such as IgG Fc, serum albumin and AAV capsids. Verbatim, the thioether
  ring-closure moiety is *replaced by a well-folded natural protein domain*. These are
  display-derived macrocycles, not lasso peptides, and not derived from any natural lasso
  peptide. Composition of the three related peptides (lower-case = D-amino acid; each
  N-chloroacetylated and thioether-cyclised to its C-terminal Cys):

  | Peptide | Free-peptide sequence | KD vs human PlxnB1 (SPR) |
  |---|---|---|
  | PB1m6 | Ac-**w**RPRVARWTGQIIYC | 3.5 nM |
  | PB1m6A9 (affinity-matured PB1m6) | Ac-**w**RPYIERWTGRLIVC | 0.28 nM (44 nM mouse) |
  | PB1m7 | Ac-**w**NSNVLSWQTYSWYC | ~275–300 nM |

  The 7VF3/7VG7 crystals contain a single-chain tandem **uteroglobin (SCGB1A1, P11684)
  with the peptide inserted between the two uteroglobin units** — uteroglobin was chosen
  only as a soluble crystallization carrier. 7VF3 insert = `CNSNVLSWQTYSWYC` (its only
  disulfide is the two *added* flanking Cys, so this graft is a disulfide-closed ring);
  7VG7 insert = `WRPYIERWTGRLIV` (no Cys, uncyclised). The therapeutic format is the
  **IgG1 Fc** fusion, where the two siblings behave oppositely in the PlxnB1 cell-collapse
  assay: **PB1m7-Fc is the agonist (EC50 0.30 nM)** and **PB1m6A9-Fc the antagonist
  (IC50 ≈ 10.3 nM)**. No 1:1 Kd is reported for the divalent Fc fusions. The complexes are
  real (2.29 Å and 2.50 Å) but they are not lasso-peptide positives. (Potency figures come
  from the author's open dissertation, doi:10.18910/82338, because the Structure 2022 paper
  is closed access.) PDB 5B4W is the parent PB1m6–PlxnB1 complex — also a RaPID macrocycle,
  not a lasso peptide.

### Unverified / needs checking
| Item | Status |
|---|---|
| **Capi-var1** `GTPGFQTPDNRVISRFGFN` | Computational variant (capistruin D8N, LassoPred LP_318). No experimental evidence of any kind. |
| **Chaxapeptin** `GFGSKPLDSFGLNFF` | Genuine lasso peptide — PDB **2N5C**, isopeptide Gly1 α-NH2 → Asp8 β-COOH (8-residue ring + 7-residue tail), verified from coordinates. **No molecular target.** Reject three circulating claims: the "MIC 30–35 µg/mL" is a DRAMP assertion absent from the primary text; a supposedly separate lung-cancer paper does not exist (43 Europe PMC records checked); and an INTEDE entry naming *S. aureus* nfrA as a target is a text-mining artefact with zero supporting records. Discovery paper author order is Elsayed SS et al. (Rateb ME last), J Org Chem 2015;80:10252, doi:10.1021/acs.joc.5b01878. |
| **Sphingopyxin-I** `GIEPLGPVDEDQGEHYLFAGG` | Verified lasso sequence (isopeptide Gly1→Asp9; 9-residue ring, 12-residue tail; UniProt A0A1D5B387; 5JQF is classified "unknown function"). **Target unknown and no bioactivity ever reported** — not even an MIC. Discovery: Hegemann et al., Biopolymers 2013;100(5):527-542, doi:10.1002/bip.22326. Its isopeptidase SpI-IsoP (5JRK/5JRL) is its own catabolism/self-immunity enzyme, **not** a target. |
| **Lassomycin affinity** | Gavrish 2014 used **no direct binding assay** — no SPR, ITC, MST, NMR titration, pull-down or thermal shift. The 0.41 µM figure is an *apparent* constant from a Hill fit of ATPase **activation**, i.e. functional, and must not be reported as a measured Kd. No lassomycin–ClpC1 Kd exists in the accessible literature. Resistance maps to *clpC1* (Gln17→Arg/His, Arg21→Ser, Pro79→Thr). Structures 8IBO (1.83 Å) and 8IBP (1.45 Å) report no Kd. |
| **Lassomycin sequence and fold** | The three depositions disagree: 2MAI = `GLRRLFANQLVGRRN` + C-terminal methyl ester (ILM = methyl L-isoleucinate), 8IBO = `GLRRLFADQLVGRRN`, 8IBP = `GLRRLFADQAVGRRN`. 2MAI records the macrolactam as Gly1–**Asn8** while the paper and 8IBO/8IBP give **Asp8**; whether this is an uncorrected deposition error is unknown, and 8IBP carries an unexplained Ala10 vs Leu10. Gavrish 2014 states lassomycin "lacks the characteristic knot … the C-terminal end packs tightly against the N-terminal ring instead of passing through the macrolactam" — an *unthreaded* lasso-fold peptide — although total synthesis work argues it may in fact be threaded. Its status as a canonical lasso peptide is therefore arguable even though the ClpC1 complex is solid. |
| **Siamycin family identity** | Three distinct 21-mers: siamycin I = MS-271 = BMY-29304 = NP-06, `CLGVGSCNDFAGCGYAIVCFW` (V4/I17, UniProt P85078); siamycin II, `CLGIGSCNDFAGCGYAIVCFW` (I4/I17); RP 71955 = aborycin, `CLGIGSCNDFAGCGYAVVCFW` (I4/V17, UniProt P37046). Do not merge. Only RP 71955 has a structure (1RPB/1RPC); the "MS-271" hits 8ITG/8ITH/8GQA are the epimerase MslH, not the peptide. |
| **Siamycin × HIV gp41 mechanism** | Fusion blockade is established (Lin 1996, ED50 0.08 µM, resistance maps to gp160, gp120-CD4 binding refuted). The specific **gp41 six-helix-bundle / HR1-HR2** mechanism is untested — zero experimental support, labelled "(HYPOTHESIS)" in a 2021 review. Kept at grade C. |
| **Siamycin is promiscuous** | It inhibits nine other *E. faecalis* sensor kinases, *R. sphaeroides* PrrB, a porcine Na⁺-ATPase and bovine PKA. A quoted 10.9 µM IC50 belongs to **sviceucin**, not siamycin I. |
| **MccJ25(RGD/RGDF) IC50 provenance** | Both studies report only nmolar values, but in **opposite ELISA orientations** (2011: immobilised integrin + biotinylated vitronectin; 2014: immobilised ligand + soluble integrin), which explains 17.2 vs 12 nM (αvβ3) and 169.4 vs 38 nM (αvβ5) — quote each with its study. Compound **identity** is resolved: the three ChEMBL entries correspond to real grafted MccJ25 variants, not free peptides (Hegemann thesis SI: isolated yields 0.2 mg/L for MccJ25(FRGD), 4.6 mg/L for MccJ25(RGDF)). A **unit discrepancy** was reported in ChEMBL3297736 (displayed as µM against the paper's own "[nM]" footnote) — take nM values from the paper and check units before using ChEMBL. |
| **Lassotide IC50 tables** | The verified numbers (lassotides 44–46 IC50 1–2 nM vs αvβ6/αvβ8; lassotide 14 <100 nM; lassotide 11 >10 µM) are **preprint running-text** statements. The per-compound IC50 **tables** could not be read — bioRxiv, the ACS SI, the Wayback Machine and two proxies all returned Cloudflare 403/429/522 — so both the numbers *and the units* inside those tables are unread. Do not use any per-compound lassotide number. 9NY3/9NWZ are 20-residue, 15-model NMR structures of the peptides alone, annotated "designed, de novo protein"; their mapping onto commercial lassotide numbers (19/36/47) is unverified. |
| **"Anantin C"** `GFIGWGDDIFGHYSGDF` | Could not be verified against any primary source. Plain anantin **is** verified. |
| **RES-701-3 IC50** | Conflict: 31.5 nM (2025 Nat Commun) vs a 4 nM literature value cited there; the 1995 series paper gives only a 5–20 nM range for RES-701-2/-3/-4 collectively. |
| **Propeptin × PREP (P48147)** | Claimed in `docs/lasso_known_targets.md` as grade B; not independently verified in this pass. |
| **Lassomycin × human CLPB** | Homology inference only — retained at grade C, never upgraded. |

---

## 5. Files

* `docs/lasso_known_pairs.csv` — 25 verified pairs, columns: `peptide, pep_seq, target,
  target_uniprot, organism, evidence_grade, pdb_ids, assay, affinity_value,
  affinity_units, reference, notes`.
* `docs/lasso_known_pairs.md` — this note.
