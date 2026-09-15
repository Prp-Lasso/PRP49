"""Preflight checks: run BEFORE submitting any long job.

Turns the recurring mistakes from docs/PRP49_问题总结与教训.md into automatic
checks, so prevention does not depend on remembering 31 lessons.

Usage examples
--------------
  python scratch/preflight.py --data mvp_cpu/train_pairs_v2.csv
  python scratch/preflight.py --script hpc_deploy/scripts/job_reg_bound2.sh \
        --data mvp_cpu/affinity_pairs.csv --epochs 30 --folds 5 --min-per-fold 45
  python scratch/preflight.py --trainer PRP49/prp49/train.py --data mvp_cpu/train_pairs_prop2.csv
"""
import argparse
import os
import re
import sys

ROOT = r'D:\deepseek_harness\prp49'
OK, WARN, FAIL = 'PASS', 'WARN', 'FAIL'
results = []


def check(name, status, detail=''):
    results.append((status, name, detail))
    icon = {OK: '[ok]  ', WARN: '[warn]', FAIL: '[FAIL]'}[status]
    print(f'{icon} {name}' + (f'\n        {detail}' if detail else ''))


# ---------------------------------------------------------------- data checks
def check_data(path):
    """CRLF, empty values, required columns, target-key hygiene."""
    full = path if os.path.isabs(path) else os.path.join(ROOT, path)
    if not os.path.exists(full):
        check(f'data exists: {path}', FAIL, 'file not found')
        return
    raw = open(full, 'rb').read()
    size_mb = len(raw) / 1e6

    n_cr = raw.count(b'\r\n')
    if n_cr:
        check(f'data line endings: {path}', FAIL,
              f'{n_cr} CRLF line breaks -> sanitize with: sed -i "s/\\r$//" {path}')
    else:
        check(f'data line endings: {path}', OK, f'LF only ({size_mb:.1f} MB)')

    if len(raw) == 0:
        check(f'data non-empty: {path}', FAIL, 'file is 0 bytes')
        return

    header = raw.split(b'\n', 1)[0].decode('utf-8', 'replace').strip()
    cols = [c.strip() for c in header.split(',')]
    check(f'data header: {path}', OK, f'{len(cols)} cols: {cols[:8]}{"..." if len(cols) > 8 else ""}')

    # target-key hygiene: a chain descriptor is not a target id (#28)
    if 'prot_id' in cols:
        idx = cols.index('prot_id')
        bad = 0
        total = 0
        for line in raw.decode('utf-8', 'replace').splitlines()[1:]:
            parts = line.split(',')
            if len(parts) > idx:
                total += 1
                if re.match(r'^(chain|chains)\b', parts[idx].strip(), re.I):
                    bad += 1
        if total and bad / total > 0.2:
            check('target-key hygiene (prot_id)', FAIL,
                  f'{bad}/{total} ({100*bad/total:.0f}%) look like PDB chain descriptors '
                  f'"Chains A, B" -> key targets by SEQUENCE instead')
        elif bad:
            check('target-key hygiene (prot_id)', WARN, f'{bad}/{total} chain-descriptor values')
        else:
            check('target-key hygiene (prot_id)', OK)

    if 'propagated' in cols:
        check('evaluation isolation (propagated column)', WARN,
              'dataset contains inferred labels: the training script MUST filter them '
              'out of every test fold (train.py has eval_mask) - verify it is active')


# --------------------------------------------------------------- script checks
def check_script(path, epochs=None, folds=None, min_per_fold=None):
    full = path if os.path.isabs(path) else os.path.join(ROOT, path)
    if not os.path.exists(full):
        check(f'script exists: {path}', FAIL, 'not found')
        return
    s = open(full, encoding='utf-8', errors='replace').read()

    m = re.search(r'#SBATCH --time=(\d+)-(\d+):(\d+):(\d+)', s) or \
        re.search(r'#SBATCH --time=(\d+):(\d+):(\d+)', s)
    limit_h = None
    if m:
        g = [int(x) for x in m.groups()]
        limit_h = (g[0] * 24 + g[1] + g[2] / 60) if len(g) == 4 else (g[0] + g[1] / 60)
        check(f'wall-clock limit: {os.path.basename(path)}', OK, f'{limit_h:.1f} h')

    if folds and min_per_fold and limit_h:
        # --min-per-fold means the measured wall-clock for ONE fold (all epochs
        # included). An earlier version multiplied by epochs as well, double
        # counting and reporting nonsense estimates - the checker itself had a bug,
        # which is exactly the "verification is broken" failure mode (#17/#18/#26).
        need = folds * min_per_fold / 60.0
        if need > limit_h:
            check('time-limit feasibility', FAIL,
                  f'estimated need {need:.1f} h (folds {folds} x {min_per_fold:.0f} min/fold) '
                  f'> limit {limit_h:.1f} h -> raise --time or cut epochs/folds')
        elif need * 1.4 > limit_h:
            check('time-limit feasibility', WARN,
                  f'need ~{need:.1f} h vs limit {limit_h:.1f} h - margin under 40%, '
                  f'a slow fold or a queue restart can still overflow it')
        else:
            check('time-limit feasibility', OK, f'need ~{need:.1f} h vs limit {limit_h:.1f} h')
    elif folds or min_per_fold:
        check('time-limit feasibility', WARN, 'pass --folds and --min-per-fold for the estimate')

    # any array job that runs several heavy configs must be split (#27)
    n_cases = len(re.findall(r'^\s*python -m prp49\.', s, re.M))
    has_array = '--array' in s
    if n_cases > 1 and not has_array:
        check('multi-config job split', FAIL,
              f'{n_cases} training runs in ONE job without --array -> a timeout kills all of them (#27)')
    elif n_cases > 1:
        check('multi-config job split', OK, f'{n_cases} runs, array job')
    else:
        check('multi-config job split', OK, f'{n_cases} training invocation')


def check_trainer(path):
    """Checkpoint coverage + eval isolation + baseline reporting."""
    full = path if os.path.isabs(path) else os.path.join(ROOT, path)
    if not os.path.exists(full):
        check(f'trainer exists: {path}', FAIL, 'not found')
        return
    s = open(full, encoding='utf-8', errors='replace').read()

    if 'save_resume(' in s and 'maybe_resume(' in s:
        check('periodic checkpoint + resume', OK,
              f'save_resume x{s.count("save_resume(")}, maybe_resume x{s.count("maybe_resume(")}')
    else:
        check('periodic checkpoint + resume', FAIL,
              'no save_resume/maybe_resume -> a wall-clock kill discards completed epochs (#30)')

    if 'drop_last' in s:
        check('BatchNorm batch>1 guard', OK, 'drop_last present')
    else:
        check('BatchNorm batch>1 guard', WARN, 'no drop_last: a final batch of size 1 will crash BN (#12)')

    if path.endswith('train_reg_bound.py'):
        for token, label in [('baseline', 'mean-predictor baseline'), ('pearsonr', 'Pearson reporting')]:
            check(f'regression {label}', OK if token in s else FAIL,
                  '' if token in s else f'missing "{token}" -> RMSE alone can hide a useless model (#D3)')
    if 'all_scores[eval_te]' in s or 'eval_mask' in s:
        check('scatter/overall index consistency', OK, 'uses the same index for scatter and overall')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data', action='append', default=[])
    ap.add_argument('--script', action='append', default=[])
    ap.add_argument('--trainer', action='append', default=[])
    ap.add_argument('--epochs', type=int)
    ap.add_argument('--folds', type=int)
    ap.add_argument('--min-per-fold', type=float,
                    help='measured wall-clock minutes for ONE fold (all epochs included)')
    args = ap.parse_args()

    print('=' * 72)
    print('PRP49 preflight - run this BEFORE submitting a long job')
    print('=' * 72)
    for d in args.data:
        check_data(d)
    for t in args.trainer:
        check_trainer(t)
    for s in args.script:
        check_script(s, args.epochs, args.folds, args.min_per_fold)

    n_fail = sum(1 for st, _, _ in results if st == FAIL)
    n_warn = sum(1 for st, _, _ in results if st == WARN)
    print('=' * 72)
    print(f'SUMMARY: {len(results)} checks | {n_fail} FAIL | {n_warn} WARN')
    if n_fail:
        print('=> fix the FAIL items before submitting; each maps to a mistake already paid for.')
    print('=' * 72)
    return 1 if n_fail else 0


if __name__ == '__main__':
    sys.exit(main())
