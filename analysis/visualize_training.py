"""Visualise the current training stage from the downloaded logs.

Produces figures/*.png:
  1. rmse_schemes.png       - RMSE vs epoch for schemes A/B/C, with per-fold baselines
  2. cv_ladder.png          - CV protocol ladder + all model versions
  3. routeA_progress.png    - route A (propagation) fold-by-fold AUC vs baseline
  4. mtl_domains.png        - multitask in-domain vs cross-domain per fold
  5. loss_curves.png        - training loss curves
"""
import os
import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = r'D:\deepseek_harness\prp49'
LOGS = os.path.join(ROOT, 'logs')
FIG = os.path.join(ROOT, 'figures')
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({'figure.dpi': 130, 'font.size': 10, 'axes.grid': True,
                     'grid.alpha': 0.3, 'axes.spines.top': False, 'axes.spines.right': False})


def read(name):
    p = os.path.join(LOGS, name)
    return open(p, encoding='utf-8', errors='replace').read() if os.path.exists(p) else ''


# ---------------- parse helpers ----------------
# classification logs:  fold 0 ep 0: loss 1.7259 val AUC 0.650 AP 0.242 (align batches 31)
EP = re.compile(r'fold (\d+) ep (\d+): loss ([\d.]+) val (\w+) ([-\d.]+) (\w+) ([-\d.]+)')
# regression logs:      fold 0 ep 40: loss 0.0403 RMSE 1.593 (mean-baseline 2.295)
EP_REG = re.compile(r'fold (\d+) ep (\d+): loss ([\d.]+) RMSE ([\d.]+) \(mean-baseline ([\d.]+)\)')
FOLD_CLS = re.compile(r'^fold (\d+): AUC ([\d.]+) AP ([\d.]+)', re.M)
FOLD_REG = re.compile(r'^fold (\d+): RMSE ([\d.]+) vs mean-baseline ([\d.]+) \| per-target median '
                      r'([\d.]+) \| pred range \[([-\d.]+), ([-\d.]+)\] \| r=([-\d.]+)', re.M)
FOLD_MTL = re.compile(r'^fold (\d+): in-domain ([\d.]+) \| cross-domain ([\d.]+)', re.M)


def parse_epochs(txt, kind='cls'):
    """kind='cls' -> val metric + loss;  kind='reg' -> RMSE + baseline + loss."""
    if kind == 'reg':
        return [dict(fold=int(m.group(1)), ep=int(m.group(2)), loss=float(m.group(3)),
                     value=float(m.group(4)), baseline=float(m.group(5)))
                for m in EP_REG.finditer(txt)]
    return [dict(fold=int(m.group(1)), ep=int(m.group(2)), loss=float(m.group(3)),
                 metric=m.group(4), value=float(m.group(5)))
            for m in EP.finditer(txt)]


# ---------------- 1. RMSE schemes ----------------
fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.0))
schemes = {'A: full + tanh + base_head': 'regb2_62573246_1.out',
           'B: Kd+Ki only': 'regb2_62573246_2.out',
           'C: full + tanh (no base_head)': 'regb2_62573246_3.out'}
colors = ['#1f77b4', '#2ca02c', '#d62728']
for ax, (label, fname), col in zip(axes, schemes.items(), colors):
    txt = read(fname)
    eps = parse_epochs(txt, kind='reg')
    folds = sorted({e['fold'] for e in eps})
    for i, fold in enumerate(folds):
        sub = sorted([e for e in eps if e['fold'] == fold], key=lambda e: e['ep'])
        ax.plot([e['ep'] for e in sub], [e['value'] for e in sub],
                marker='o', ms=3, lw=1.6, color=col, alpha=0.45 + 0.25 * i,
                label=f'fold {fold}')
        if sub and 'baseline' in sub[0]:
            ax.axhline(sub[0]['baseline'], ls=':', lw=1.2, color=col, alpha=0.6)
    if folds:
        ax.set_xlim(0, max(e['ep'] for e in eps) + 2)
    ax.set_title(label, fontsize=9.5)
    ax.set_xlabel('epoch'); ax.set_ylabel('val RMSE (pAffinity units)')
    ax.legend(fontsize=7.5)
axes[0].text(0.02, 0.03, 'dotted = per-fold mean-predictor baseline',
             transform=axes[0].transAxes, fontsize=7.5, color='gray')
fig.suptitle('RMSE-constrained regression: bounded head stops the runaway (was RMSE 6.86)',
             fontsize=11.5, y=1.02)
fig.tight_layout(); fig.savefig(os.path.join(FIG, 'rmse_schemes.png'), bbox_inches='tight')
plt.close(fig)
print('  1/5 rmse_schemes.png')

# ---------------- 2. CV ladder / model versions ----------------
versions = [('8M MVP\n(CPU)', 0.663, '#bbbbbb'),
            ('650M random\ninit baseline', 0.594, '#bbbbbb'),
            ('I1 curriculum\nonly', 0.566, '#bbbbbb'),
            ('I1+I2+I5\n(improved)', 0.8145, '#4c72b0'),
            ('peptide-grouped\n(headline)', 0.8333, '#4c72b0'),
            ('family-grouped\n(strictest)', 0.8142, '#55a868'),
            ('route A\npropagation', 0.726, '#c44e52'),
            ('A+C affinity\n(new target)', 0.6312, '#8172b2'),
            ('cross-domain\n(single task)', 0.470, '#937860'),
            ('cross-domain\n(multitask)', 0.6421, '#da8bc3')]
fig, ax = plt.subplots(figsize=(11, 4.2))
xs = np.arange(len(versions))
vals = [v[1] for v in versions]
ax.bar(xs, vals, color=[v[2] for v in versions], width=0.68)
ax.axhline(0.75, ls='--', color='green', lw=1.4)
ax.text(len(versions) - 0.4, 0.762, 'acceptance 0.75', color='green', fontsize=8.5, ha='right')
ax.axhline(0.5, ls=':', color='gray', lw=1.2)
ax.text(len(versions) - 0.4, 0.512, 'random', color='gray', fontsize=8.5, ha='right')
for x, v in zip(xs, vals):
    ax.text(x, v + 0.012, f'{v:.3f}', ha='center', fontsize=8.5)
ax.set_xticks(xs); ax.set_xticklabels([v[0] for v in versions], fontsize=8)
ax.set_ylabel('AUC'); ax.set_ylim(0, 0.95)
ax.set_title('All model versions and evaluation protocols (route A underperforms its baseline)',
             fontsize=11)
fig.tight_layout(); fig.savefig(os.path.join(FIG, 'cv_ladder.png'), bbox_inches='tight')
plt.close(fig)
print('  2/5 cv_ladder.png')

# ---------------- 3. route A progress ----------------
txt = read('prp49_prop-62593230.out')
folds_a = [(int(m.group(1)), float(m.group(2))) for m in FOLD_CLS.finditer(txt)]
eps_a = parse_epochs(txt)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.0))
if folds_a:
    x, y = zip(*folds_a)
    ax1.bar([f'fold {i}' for i in x], y, color='#c44e52', width=0.55)
    ax1.axhline(0.8142, ls='--', color='#55a868', lw=1.6, label='family-CV baseline 0.8142')
    ax1.axhline(0.75, ls=':', color='green', lw=1.2, label='acceptance 0.75')
    for xi, yi in zip(range(len(x)), y):
        ax1.text(xi, yi + 0.012, f'{yi:.3f}', ha='center', fontsize=9)
    ax1.set_ylim(0, 1.0); ax1.set_ylabel('val AUC')
    ax1.set_title('Route A (propagation) finished folds', fontsize=10)
    ax1.legend(fontsize=8)
for fold in sorted({e['fold'] for e in eps_a}):
    sub = sorted([e for e in eps_a if e['fold'] == fold], key=lambda e: e['ep'])
    if sub:
        ax2.plot([e['ep'] for e in sub], [e['value'] for e in sub], marker='o', ms=4,
                 label=f"fold {fold} AUC")
        ax2.plot([e['ep'] for e in sub], [e['loss'] for e in sub], marker='s', ms=3,
                 ls='--', alpha=0.5, label=f"fold {fold} loss")
ax2.set_xlabel('epoch'); ax2.set_title('Route A: validation AUC and training loss', fontsize=10)
ax2.legend(fontsize=7.5)
fig.suptitle('Route A: family label propagation (352 -> 759 positives) does NOT beat the baseline',
             fontsize=11, y=1.02)
fig.tight_layout(); fig.savefig(os.path.join(FIG, 'routeA_progress.png'), bbox_inches='tight')
plt.close(fig)
print('  3/5 routeA_progress.png')

# ---------------- 4. MTL domains ----------------
txt = read('prp49_mtl-62558477.out')
mtl = [(int(m.group(1)), float(m.group(2)), float(m.group(3))) for m in FOLD_MTL.finditer(txt)]
fig, ax = plt.subplots(figsize=(8, 4.2))
if mtl:
    x = np.arange(len(mtl)); w = 0.36
    ax.bar(x - w / 2, [m[1] for m in mtl], w, label='in-domain (lasso peptides)', color='#4c72b0')
    ax.bar(x + w / 2, [m[2] for m in mtl], w, label='cross-domain (ordinary peptides)', color='#da8bc3')
    ax.axhline(0.8333, ls='--', color='#4c72b0', lw=1.2, alpha=0.7)
    ax.axhline(0.470, ls='--', color='#da8bc3', lw=1.2, alpha=0.7)
    ax.text(len(mtl) - 0.5, 0.845, 'single-task in-domain 0.833', fontsize=8, ha='right', color='#4c72b0')
    ax.text(len(mtl) - 0.5, 0.482, 'single-task cross-domain 0.470', fontsize=8, ha='right', color='#da8bc3')
    ax.set_xticks(x); ax.set_xticklabels([f'fold {m[0]}' for m in mtl])
    ax.set_ylabel('AUC'); ax.set_ylim(0, 1.0)
    ax.set_title(f'Multitask training: cross-domain 0.470 -> '
                 f'{np.mean([m[2] for m in mtl]):.3f} (mean of {len(mtl)} folds)')
    ax.legend(fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(FIG, 'mtl_domains.png'), bbox_inches='tight')
plt.close(fig)
print('  4/5 mtl_domains.png')

# ---------------- 5. loss curves ----------------
fig, ax = plt.subplots(figsize=(8.5, 4.2))
txt = read('prp49_prop-62593230.out')
for fold in sorted({e['fold'] for e in parse_epochs(txt)}):
    sub = sorted([e for e in parse_epochs(txt) if e['fold'] == fold], key=lambda e: e['ep'])
    ax.plot([e['ep'] for e in sub], [e['loss'] for e in sub], marker='o', ms=4,
            label=f'route A fold {fold}')
for label, fname, col in [('reg A fold 0', 'regb2_62573246_1.out', '#1f77b4'),
                          ('reg C fold 0', 'regb2_62573246_3.out', '#d62728')]:
    eps = parse_epochs(read(fname), kind='reg')
    sub = sorted([e for e in eps if e['fold'] == 0], key=lambda e: e['ep'])
    if sub:
        ax.plot([e['ep'] for e in sub], [e['loss'] for e in sub], marker='s', ms=3,
                ls='--', color=col, label=label)
ax.set_xlabel('epoch'); ax.set_ylabel('training loss'); ax.set_yscale('log')
ax.set_title('Training loss curves (log scale)')
ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(FIG, 'loss_curves.png'), bbox_inches='tight')
plt.close(fig)
print('  5/5 loss_curves.png')

print(f'\nfigures written to {FIG}')
