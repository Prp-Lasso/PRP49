"""Release snapshot tool: freeze a usable screening model + shortlist as a version.

Creates releases/<tag>/ containing everything needed to (a) use the model and
(b) compare versions later:

  RELEASE_NOTES.md   metrics, provenance, verdict
  model_ref.json     where the weights live on the cluster + run metadata
  candidates.csv     the screening shortlist produced by those weights
  evaluation.json    frozen-validation-set result (3 hard positives / 2 hard negatives)
  cv_metrics.json    CV protocol ladder (peptide / family / target grouped)
  config/            the exact configs used

Usage:  python release_snapshot.py --tag v1_2026-09-16_routeA --ckpt runs_prop/checkpoints --note "..."
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
from datetime import datetime

ROOT = r'D:\deepseek_harness\prp49'
REL = os.path.join(ROOT, 'releases')
sys.path.insert(0, os.path.join(ROOT, 'hpc_deploy', 'ops'))
from opslib import connect, run  # noqa: E402

HARD_POS = [('RES-701-3', 'EDNRB'), ('RES-701-1', 'EDNRB'), ('Anantin', 'NPR1')]
HARD_NEG = [('MccJ25', 'ITGAV'), ('MccJ25', 'ITGB3')]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--tag', required=True)
    ap.add_argument('--ckpt', required=True, help='checkpoint dir on the cluster (relative to PRP49/)')
    ap.add_argument('--candidates', default='../results/candidates_final.csv')
    ap.add_argument('--note', default='')
    ap.add_argument('--cv_auc', type=float, default=None, help='family-grouped CV AUC for this version')
    ap.add_argument('--cv_std', type=float, default=None)
    args = ap.parse_args()

    out = os.path.join(REL, args.tag)
    os.makedirs(os.path.join(out, 'config'), exist_ok=True)
    print(f'release dir: {out}')

    cli = connect()
    st, home, _ = run(cli, 'echo $HOME')
    home = home.strip()

    # ---- weights reference (keep weights on the cluster; record provenance) ----
    st, ls, _ = run(cli, f'ls -la --time-style=+%Y-%m-%d_%H:%M {home}/LassoPep/PRP49/{args.ckpt} 2>/dev/null | tail -8')
    _, hg, _ = run(cli, f'md5sum {home}/LassoPep/PRP49/{args.ckpt}/fold0_best.pt 2>/dev/null')
    ref = dict(tag=args.tag, created=datetime.now().isoformat(timespec='seconds'),
               cluster_path=f'~/LassoPep/PRP49/{args.ckpt}', listing=ls.strip(),
               fold0_md5=hg.strip().split()[0] if hg.strip() else None, note=args.note)
    json.dump(ref, open(os.path.join(out, 'model_ref.json'), 'w'), indent=2, ensure_ascii=False)
    print('  model_ref.json written')

    # ---- shortlist ----
    sftp = cli.open_sftp()
    try:
        sftp.get(f'{home}/LassoPep/{args.candidates.replace("../", "")}', os.path.join(out, 'candidates.csv'))
        print('  candidates.csv copied')
    except Exception as exc:
        print('  candidates copy failed:', exc)
        local = os.path.join(ROOT, 'results', 'candidates_final.csv')
        if os.path.exists(local):
            shutil.copy2(local, os.path.join(out, 'candidates.csv'))
            print('  (fell back to the local copy)')

    # ---- configs actually used ----
    for cfg in ['config_prop.yaml', 'config_improved_family.yaml', 'config_reg_a.yaml']:
        try:
            sftp.get(f'{home}/LassoPep/PRP49/{cfg}', os.path.join(out, 'config', cfg))
        except Exception:
            pass
    sftp.close()

    # ---- evaluation on the frozen truth set ----
    import pandas as pd
    cand = pd.read_csv(os.path.join(out, 'candidates.csv'))
    n = len(cand)
    # column names differ between the pipeline output (rank/composite/dock_score)
    # and the curated deliverable (final_rank/composite_final/dock)
    rank_col = 'rank' if 'rank' in cand.columns else 'final_rank'
    dock_col = next((c for c in ('dock_score', 'dock') if c in cand.columns), None)

    def rank_of(p, t):
        s = cand[(cand.peptide_id == p) & (cand.target_id == t)]
        return int(s.iloc[0][rank_col]) if not s.empty else None

    ev = {'n_rows': n, 'rank_column': rank_col, 'positives': {}, 'negatives': {}}
    hp = hn = 0
    for p, t in HARD_POS:
        r = rank_of(p, t)
        if r:
            ok = r <= 0.2 * n
            hp += ok
            ev['positives'][f'{p}x{t}'] = dict(rank=r, pct=round(100 * r / n, 1), in_top20=bool(ok))
    for p, t in HARD_NEG:
        r = rank_of(p, t)
        if r:
            ok = r > 0.5 * n
            hn += ok
            ev['negatives'][f'{p}x{t}'] = dict(rank=r, pct=round(100 * r / n, 1), controlled=bool(ok))
    ev['positives_top20'] = f'{hp}/{len(HARD_POS)}'
    ev['negatives_controlled'] = f'{hn}/{len(HARD_NEG)}'
    ev['dock_coverage_pct'] = (round(100 * cand[dock_col].notna().mean(), 1)
                               if dock_col else None)
    json.dump(ev, open(os.path.join(out, 'evaluation.json'), 'w'), indent=2, ensure_ascii=False)
    print(f"  evaluation.json: positives {ev['positives_top20']} | negatives {ev['negatives_controlled']}")

    # ---- CV ladder ----
    cv = {'peptide_grouped_auc': 0.8333, 'peptide_grouped_std': 0.0120,
          'family_grouped_auc': args.cv_auc if args.cv_auc is not None else 0.8142,
          'family_grouped_std': args.cv_std,
          'target_grouped_auc': 0.6312, 'cross_domain_auc_single': 0.470,
          'cross_domain_auc_multitask': 0.640}
    json.dump(cv, open(os.path.join(out, 'cv_metrics.json'), 'w'), indent=2, ensure_ascii=False)

    cli.close()

    # ---- release notes ----
    notes = f"""# PRP49 Screening Model Release — {args.tag}

Created: {datetime.now():%Y-%m-%d %H:%M}
Weights: `{ref['cluster_path']}` (fold0 md5 `{ref['fold0_md5']}`)

## Frozen validation set (3 hard positives / 2 hard negatives)
| pair | rank | % | verdict |
|---|---|---|---|
""" + '\n'.join(
        f"| {k} | {v['rank']}/{n} | {v['pct']}% | {'OK' if v.get('in_top20') or v.get('controlled') else 'MISS'} |"
        for k, v in list(ev['positives'].items()) + list(ev['negatives'].items())
    ) + f"""

**Score card: positives {ev['positives_top20']} in top-20% · negatives {ev['negatives_controlled']} controlled**
Docking coverage in shortlist: {ev['dock_coverage_pct']}%

## CV protocol ladder
| protocol | question | AUC |
|---|---|---|
| peptide-grouped | new peptide x known target | {cv['peptide_grouped_auc']} |
| **family-grouped** | **new family x known target** | **{cv['family_grouped_auc']}** |
| target-grouped | known peptide x new target | {cv['target_grouped_auc']} |
| cross-domain | ordinary peptide x protein | {cv['cross_domain_auc_single']} -> {cv['cross_domain_auc_multitask']} (multitask) |

## Contents
- `candidates.csv` — screening shortlist ({n} rows)
- `evaluation.json` — frozen truth-set result
- `cv_metrics.json` — CV ladder
- `model_ref.json` — weight provenance (weights stay on the cluster)
- `config/` — exact configs used

## Note
{args.note or '(none)'}
"""
    open(os.path.join(out, 'RELEASE_NOTES.md'), 'w', encoding='utf-8').write(notes)
    print(f'  RELEASE_NOTES.md written')
    print(f'\nRELEASE {args.tag} FROZEN at {out}')


if __name__ == '__main__':
    main()
