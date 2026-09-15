# PRP49 Screening Model Release — v0_2026-09-15_baseline

Created: 2026-09-15 14:10
Weights: `~/LassoPep/PRP49/runs_improved_grouped/checkpoints` (fold0 md5 `d0c8bcd8c6a51cab2f138ef5a4eda057`)

## Frozen validation set (3 hard positives / 2 hard negatives)
| pair | rank | % | verdict |
|---|---|---|---|
| RES-701-3xEDNRB | 11/121 | 9.1% | OK |
| RES-701-1xEDNRB | 10/121 | 8.3% | OK |
| AnantinxNPR1 | 20/121 | 16.5% | OK |
| MccJ25xITGAV | 68/121 | 56.2% | OK |
| MccJ25xITGB3 | 69/121 | 57.0% | OK |

**Score card: positives 3/3 in top-20% · negatives 2/2 controlled**
Docking coverage in shortlist: 90.9%

## CV protocol ladder
| protocol | question | AUC |
|---|---|---|
| peptide-grouped | new peptide x known target | 0.8333 |
| **family-grouped** | **new family x known target** | **0.8142** |
| target-grouped | known peptide x new target | 0.6312 |
| cross-domain | ordinary peptide x protein | 0.47 -> 0.64 (multitask) |

## Contents
- `candidates.csv` — screening shortlist (121 rows)
- `evaluation.json` — frozen truth-set result
- `cv_metrics.json` — CV ladder
- `model_ref.json` — weight provenance (weights stay on the cluster)
- `config/` — exact configs used

## Note
Baseline release: current working screening model. Peptide-grouped CV 0.8333, family-grouped CV 0.8142 (no family leakage), frozen validation 3/3 positives top-20% and 2/2 negatives controlled. Weights stay on the cluster.
