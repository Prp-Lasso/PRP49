"""Version alignment: pin down which weights each release actually ships, under which
CV protocol, so the 9-18 comparison is apples-to-apples.

The problem found: the acceptance criterion names family-grouped CV as the primary
metric (runs_family = 0.8142), but the v0 snapshot records the checkpoint md5 of
runs_improved_grouped (0.8333, grouped by peptide). Same model architecture, different
protocol - so a naive comparison would silently mix two protocols.
"""
import json
import os

import pandas as pd

ROOT = r'D:\deepseek_harness\prp49'
R = os.path.join(ROOT, 'releases')

RUNS = {
    'runs_family': dict(auc=0.8142, protocol='family-grouped (held-out FAMILIES)',
                        folds=[0.8518, 0.8673, 0.7770, 0.8192, 0.7554],
                        md5='1659dc6027e2b0e1acf92c8ee5803b1e'),
    'runs_improved_grouped': dict(auc=0.8333, protocol='grouped by PEPTIDE (held-out peptides)',
                                  folds=[0.8123, 0.8490, 0.8369, 0.8375, 0.8308],
                                  md5='d0c8bcd8c6a51cab2f138ef5a4eda057'),
    'runs_improved': dict(auc=0.8145, protocol='stratified (LEAKY - reference only)',
                          folds=[0.8029, 0.7974, 0.8123, 0.8124, 0.8472], md5='-'),
    'runs_ctrft_fold*': dict(auc=0.8248, protocol='grouped by peptide (route C encoder)',
                             folds=[0.8116, 0.8443, 0.8177, 0.8313, 0.8192], md5='-'),
    'runs_prop': dict(auc=0.7481, protocol='family-grouped (route A propagation)',
                      folds=[0.7257, 0.7264, 0.6767, 0.8052, 0.8068], md5='-'),
    'runs_xneg': dict(auc=0.7940, protocol='grouped (xneg specificity training)',
                      folds=[0.7731, 0.8124, 0.7742, 0.7972, 0.8131], md5='-'),
    'runs_s2': dict(auc=0.7539, protocol='grouped (Propedia cross-topology transfer)',
                    folds=[0.7408, 0.7550, 0.7488, 0.7701, 0.7545], md5='-'),
    'runs_reg_cls': dict(auc=0.6430, protocol='by TARGET (affinity classification A+C)',
                         folds=[0.6166, 0.6517, 0.6428, 0.6754, 0.6284], md5='-'),
    'runs_warmup': dict(auc=0.6391, protocol='grouped (S1 curriculum warmup)',
                        folds=[0.6482, 0.6353, 0.6397, 0.6106, 0.6617], md5='-'),
}

print('=== all trained runs, by protocol ===')
rows = []
for name, d in RUNS.items():
    rows.append(dict(run=name, auc=d['auc'], protocol=d['protocol'],
                     std=round(float(pd.Series(d['folds']).std()), 4)))
t = pd.DataFrame(rows).sort_values('auc', ascending=False)
print(t.to_string(index=False))

print('\n=== what each release actually ships ===')
for rel in sorted(os.listdir(R)):
    p = os.path.join(R, rel)
    if not os.path.isdir(p):
        continue
    mr = os.path.join(p, 'model_ref.json')
    info = json.load(open(mr)) if os.path.exists(mr) else {}
    print(f'\n  {rel}')
    print(f"    weights      : {info.get('cluster_path', info.get('reason', '?'))[:96]}")
    print(f"    md5          : {info.get('fold0_md5', '-')}")
    print(f"    cv recorded  : {info.get('cv_auc', '-')} +- {info.get('cv_std', '-')}")
    cands = [f for f in os.listdir(p) if f.endswith('.csv')]
    print(f'    candidate csv: {cands}')

print('\n=== alignment verdict ===')
print('  v0/v1 ship runs_improved_grouped  -> grouped-by-peptide protocol -> 0.8333')
print('  but the acceptance criterion names family-grouped (0.8142) as primary.')
print('  => either relabel the releases to the grouped protocol, or add the family')
print('     checkpoint as the primary deliverable. They must not be quoted interchangeably.')
