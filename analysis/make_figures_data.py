"""Generate the data-driven figures (3, 7, 8, 9, 10) with real project numbers.

Labels are in English so the figures can be used in the paper directly; the outline
document carries the Chinese captions.
"""
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch, Polygon

ROOT = r'D:\deepseek_harness\prp49'
FIG = os.path.join(ROOT, 'docs', 'figures')
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 10,
    'axes.spines.top': False, 'axes.spines.right': False,
    'figure.dpi': 200, 'savefig.bbox': 'tight',
})
C_MAIN, C_ACC, C_WARN, C_GREY = '#2b6cb0', '#dd6b20', '#c53030', '#718096'

# ================= fig 8: CV protocol ladder =================
fig, ax = plt.subplots(figsize=(7.2, 3.6))
labels = ['Random init\n(baseline)', 'Stratified\n(leaky)', 'Grouped by\npeptide',
          'Grouped by\nfamily (strict)', 'Grouped by\ntarget (same dist.)']
vals = [0.594, 0.8145, 0.8333, 0.8142, 0.8738]
errs = [np.nan, 0.0194, 0.0120, 0.0477, 0.0571]
cols = [C_GREY, C_WARN, C_MAIN, C_MAIN, C_ACC]
bars = ax.bar(range(len(vals)), vals, color=cols, width=0.62,
              yerr=[0 if np.isnan(e) else e for e in errs], capsize=4,
              error_kw=dict(ecolor='#4a5568', lw=1))
for i, (b, v) in enumerate(zip(bars, vals)):
    ax.text(b.get_x() + b.get_width() / 2, v + 0.022, f'{v:.4f}', ha='center',
            fontsize=9, fontweight='bold')
ax.axhline(0.8142, ls='--', lw=1, color='#a0aec0')
ax.text(4.42, 0.8142, 'family-grouped\nbaseline', fontsize=7.5, color='#4a5568', va='center')
ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels, fontsize=8.5)
ax.set_ylabel('Cross-validated AUC', fontsize=10)
ax.set_ylim(0.5, 0.97)
ax.text(1, 0.62, 'leakage\ninflates', fontsize=8, color=C_WARN, ha='center', style='italic')
ax.set_title('Leakage-controlled evaluation ladder: the same model family, five protocols',
             fontsize=10.5, pad=10)
fig.savefig(os.path.join(FIG, 'fig08_cv_ladder.png'))
plt.close(fig)
print('fig08 done')

# ================= fig 9: screening funnel =================
fig, ax = plt.subplots(figsize=(7.6, 3.8))
stages = ['All pairs\n13 peptides x 50 targets', 'After hard filters\n(drug target, non-high bias,\nusable docking, top-3 in target)',
          'After diversity\n(<=2 peptides per target)']
counts = [650, 43, 10]
widths = [1.0, 0.42, 0.16]
colors = [C_GREY, C_MAIN, C_ACC]
y = 0
for i, (s, c, w, col) in enumerate(zip(stages, counts, widths, colors)):
    ytop = -i * 1.15
    ax.add_patch(Polygon([[-w, ytop], [w, ytop], [w * 0.72, ytop - 0.9], [-w * 0.72, ytop - 0.9]],
                         closed=True, facecolor=col, alpha=0.85, edgecolor='white', lw=1.5))
    ax.text(0, ytop - 0.45, f'{c}', ha='center', va='center', color='white',
            fontsize=15 if i == 0 else 13, fontweight='bold')
    ax.text(w + 0.06, ytop - 0.45, s, ha='left', va='center', fontsize=8.5)
# exclusion annotations
ann = ['-143 known Lasso targets (positive controls)\n-200 peptide bias = high\n-104 no usable docking score\n-~550 within-target rank > 3',
       '-33 kept as eligible pool,\n  only 10 carried forward']
ax.text(-1.06, -0.45, ann[0], ha='right', va='center', fontsize=7.5, color='#4a5568')
ax.text(-0.86, -1.6, ann[1], ha='right', va='center', fontsize=7.5, color='#4a5568')
ax.text(0, -3.65, 'Selection rules: drug-panel target; peptide_bias != high; docking score usable;\n'
                  'within-target rank <= 3; two signals agree; <= 2 peptides per target',
        ha='center', fontsize=8, color='#2d3748',
        bbox=dict(boxstyle='round,pad=0.45', fc='#f7fafc', ec='#cbd5e0'))
ax.set_xlim(-2.5, 2.5)
ax.set_ylim(-4.3, 0.5)
ax.axis('off')
ax.set_title('Candidate screening funnel', fontsize=11, pad=6)
fig.savefig(os.path.join(FIG, 'fig09_funnel.png'))
plt.close(fig)
print('fig09 done')

# ================= fig 10: final 10 candidates =================
top = pd.read_csv(os.path.join(ROOT, 'results', 'final_candidates_top10.csv'))
top = top.sort_values('composite_target', ascending=True)
fig, ax = plt.subplots(figsize=(7.6, 4.2))
labels = [f'{r.target_id} x {r.peptide_id}' for _, r in top.iterrows()]
vals = top.composite_target.values
rob = top.n_weightings_selected.values
cols = [C_ACC if r == 5 else (C_MAIN if r >= 4 else C_WARN) for r in rob]
bars = ax.barh(range(len(vals)), vals, color=cols, height=0.68)
for i, (b, v, r) in enumerate(zip(bars, vals, rob)):
    ax.text(v + 0.03, b.get_y() + b.get_height() / 2, f'{v:.3f}   [{int(r)}/5]',
            va='center', fontsize=8)
ax.set_yticks(range(len(labels)))
ax.set_yticklabels(labels, fontsize=8.5)
ax.set_xlabel('Fusion score (within-target standardised)', fontsize=10)
ax.set_xlim(0, 2.0)
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=C_ACC, label='robust (selected in all 5 weightings)'),
                   Patch(color=C_MAIN, label='sensitive (4/5)'),
                   Patch(color=C_WARN, label='weighting-dependent (<=3/5)')],
          fontsize=7.5, loc='lower right', frameon=False)
ax.set_title('Final 10 candidate pairs and their fusion-weight robustness', fontsize=10.5, pad=8)
fig.savefig(os.path.join(FIG, 'fig10_final_candidates.png'))
plt.close(fig)
print('fig10 done')

# ================= fig 7: family augmentation =================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.8, 3.3))
cats = ['Reported\npositives', 'After family\nhomolog augmentation']
n = [352, 759]
b = ax1.bar(cats, n, color=[C_WARN, C_MAIN], width=0.5)
for bb, v in zip(b, n):
    ax1.text(bb.get_x() + bb.get_width() / 2, v + 12, str(v), ha='center', fontweight='bold')
ax1.set_ylabel('Training positives', fontsize=9.5)
ax1.set_title('Sample augmentation', fontsize=10)
ax1.set_ylim(0, 880)

prot = ['grouped by\npeptide', 'grouped by\nfamily']
base = [0.8333, 0.8142]
aug = [np.nan, 0.7481]     # route A (propagation) measured under family-grouped CV
x = np.arange(len(prot))
ax2.bar(x - 0.19, base, width=0.36, color=C_MAIN, label='baseline')
ax2.bar(x + 0.19, [np.nan, 0.7481], width=0.36, color=C_WARN, label='with propagation (route A)')
for i, (bb, v) in enumerate(zip(x - 0.19, base)):
    ax2.text(bb, v + 0.012, f'{v:.4f}', ha='center', fontsize=8)
ax2.text(x[1] + 0.19, 0.7481 + 0.012, '0.7481', ha='center', fontsize=8, color=C_WARN,
         fontweight='bold')
ax2.set_xticks(x); ax2.set_xticklabels(prot, fontsize=8.5)
ax2.set_ylim(0.6, 0.92)
ax2.set_ylabel('CV AUC', fontsize=9.5)
ax2.legend(fontsize=7.5, frameon=False, loc='lower left')
ax2.set_title('Augmentation did not help (family CV)', fontsize=10)
fig.suptitle('Positive-sample scarcity and the family-homolog augmentation strategy',
             fontsize=10.5, y=1.02)
fig.savefig(os.path.join(FIG, 'fig07_family_augmentation.png'))
plt.close(fig)
print('fig07 done')

# ================= fig 3: target spectrum =================
fig, ax = plt.subplots(figsize=(7.4, 3.6))
groups = ['Bacterial targets', 'Human targets']
peps_b = ['MccJ25', 'Capistruin', 'Lassomycin', 'Lariocidin', 'Ubonodin', 'Klebsidin']
targ_b = ['RNAP beta-prime', 'RNAP beta-prime', 'ClpC1 (Mtb)', '30S ribosome', 'RNAP', 'RNAP']
grade_b = ['A', 'A', 'A', 'B', 'B', 'B']
peps_h = ['RES-701-3', 'RES-701-1', 'Anantin', 'BI-32169']
targ_h = ['ETB receptor', 'ETB receptor', 'NPR1', 'GCGR']
grade_h = ['A', 'B', 'B', 'B']
gmap = {'A': C_ACC, 'B': C_MAIN}

ypos, ylabs = [], []
k = 0
for pep, tg, g in zip(peps_b, targ_b, grade_b):
    ax.scatter(0, k, s=210, color=gmap[g], zorder=3)
    ax.text(0, k, g, ha='center', va='center', color='white', fontsize=8, fontweight='bold', zorder=4)
    ax.text(0.12, k, f'{pep}  →  {tg}', va='center', fontsize=8.5)
    ylabs.append(k); k += 1
k += 0.6
for pep, tg, g in zip(peps_h, targ_h, grade_h):
    ax.scatter(0, k, s=210, color=gmap[g], zorder=3)
    ax.text(0, k, g, ha='center', va='center', color='white', fontsize=8, fontweight='bold', zorder=4)
    ax.text(0.12, k, f'{pep}  →  {tg}', va='center', fontsize=8.5)
    ylabs.append(k); k += 1
ax.axhline(5.5, ls=':', color='#a0aec0')
ax.text(-0.55, 2.5, 'Bacterial', rotation=90, va='center', ha='center', fontsize=10, color=C_MAIN)
ax.text(-0.55, 9.0, 'Human', rotation=90, va='center', ha='center', fontsize=10, color=C_ACC)
ax.set_xlim(-0.75, 2.0); ax.set_ylim(-0.8, k - 0.4)
ax.set_yticks([]); ax.set_xticks([])
for s in ax.spines.values():
    s.set_visible(False)
ax.legend(handles=[Patch(color=C_ACC, label='evidence grade A (co-crystal / cryo-EM)'),
                   Patch(color=C_MAIN, label='evidence grade B (binding assay)')],
          fontsize=7.5, loc='lower right', frameon=False)
ax.set_title('Characterised lasso peptide–target pairs: bacterial vs human targets',
             fontsize=10.5, pad=8)
fig.savefig(os.path.join(FIG, 'fig03_target_spectrum.png'))
plt.close(fig)
print('fig03 done')
print(f'\nfigures written to {FIG}')
