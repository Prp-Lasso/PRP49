"""fig09 caption fix + fig10 regrouped by robustness.

fig09: the annotation read "-33 kept as eligible pool", which inverts the meaning - those
       33 pairs were DROPPED by the diversity cap (they remain available as backups in the
       audited table, but they are not in the shortlist).
fig10: sorted by score, the robustness colours interleaved (orange/blue/red/orange...), so
       the one thing the figure exists to show was hard to read. Now sorted by robustness
       group first, score second, and the legend sits outside the bars.
"""
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch, Polygon

ROOT = r'D:\deepseek_harness\prp49'
FIG = os.path.join(ROOT, 'docs', 'figures')
MAIN, ACC, WARN, GREY = '#2b6cb0', '#dd6b20', '#c53030', '#718096'

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'figure.dpi': 200, 'savefig.bbox': 'tight'})

# ---------------- fig09 ----------------
fig, ax = plt.subplots(figsize=(7.8, 3.9))
stages = ['All pairs\n13 peptides \u00d7 50 targets',
          'Hard filters\n(drug target, non-high bias,\nusable docking, top-3 in target)',
          'Diversity cap\n(\u2264 2 peptides per target)']
counts = [650, 43, 10]
widths = [1.0, 0.42, 0.16]
cols = [GREY, MAIN, ACC]
for i, (s_, c, wd, col) in enumerate(zip(stages, counts, widths, cols)):
    ytop = -i * 1.25
    ax.add_patch(Polygon([[-wd, ytop], [wd, ytop], [wd * 0.72, ytop - 0.95], [-wd * 0.72, ytop - 0.95]],
                         closed=True, facecolor=col, alpha=0.9, edgecolor='white', lw=1.6))
    ax.text(0, ytop - 0.48, f'{c}', ha='center', va='center', color='white',
            fontsize=16 if i == 0 else 13.5, fontweight='bold')
    ax.text(wd + 0.07, ytop - 0.48, s_, ha='left', va='center', fontsize=8.5)
ax.text(-1.08, -0.48,
        '\u2212143 known Lasso targets (positive controls)\n'
        '\u2212200 peptide bias = high\n'
        '\u2212104 no usable docking score\n'
        '\u2212\u2248550 within-target rank > 3',
        ha='right', va='center', fontsize=7.5, color='#4a5568')
ax.text(-0.70, -1.73,
        '\u221233 removed by the diversity cap\n(they remain in the audited table\nas backups, not in the shortlist)',
        ha='right', va='center', fontsize=7.5, color='#4a5568')
ax.text(0, -4.05,
        'Selection rules: drug-panel target; peptide bias \u2260 high; docking score usable;\n'
        'within-target rank \u2264 3; the two signals agree; \u2264 2 peptides per target',
        ha='center', fontsize=8, color='#2d3748',
        bbox=dict(boxstyle='round,pad=0.45', fc='#f7fafc', ec='#cbd5e0'))
ax.set_xlim(-2.6, 2.6)
ax.set_ylim(-4.75, 0.55)
ax.axis('off')
ax.set_title('Candidate screening funnel', fontsize=11, pad=6)
fig.savefig(os.path.join(FIG, 'fig09_funnel.png'))
plt.close(fig)
print('fig09 rebuilt (diversity-cap annotation corrected)')

# ---------------- fig10 ----------------
top = pd.read_csv(os.path.join(ROOT, 'results', 'final_candidates_top10.csv'))
# sort by robustness group first (desc), then by score (desc), then invert for barh
top['grp'] = top.n_weightings_selected
top = top.sort_values(['grp', 'composite_target'], ascending=[False, False]).reset_index(drop=True)
plot = top.iloc[::-1].reset_index(drop=True)     # barh draws bottom-up

fig, ax = plt.subplots(figsize=(8.4, 4.4))
labels = [f'{r.target_id} \u00d7 {r.peptide_id}' for _, r in plot.iterrows()]
vals = plot.composite_target.values
rob = plot.n_weightings_selected.values
cols = [ACC if r == 5 else (MAIN if r >= 4 else WARN) for r in rob]
bars = ax.barh(range(len(vals)), vals, color=cols, height=0.66)
for i, (b, v, r) in enumerate(zip(bars, vals, rob)):
    ax.text(v + 0.03, b.get_y() + b.get_height() / 2, f'{v:.3f}   [{int(r)}/5]',
            va='center', fontsize=8.5)
# group separators
groups = plot.grp.values
for i in range(1, len(groups)):
    if groups[i] != groups[i - 1]:
        ax.axhline(i - 0.5, ls=':', color='#cbd5e0', lw=1)
ax.set_yticks(range(len(labels)))
ax.set_yticklabels(labels, fontsize=8.5)
ax.set_xlabel('Fusion score (within-target standardised)', fontsize=10)
ax.set_xlim(0, 2.15)
ax.legend(handles=[Patch(color=ACC, label='robust: selected under all 5 weightings'),
                   Patch(color=MAIN, label='sensitive: 4/5 weightings'),
                   Patch(color=WARN, label='weighting-dependent: \u2264 3/5')],
          fontsize=7.5, loc='lower right', frameon=False)
ax.set_title('Final 10 candidate pairs, grouped by fusion-weight robustness', fontsize=10.5, pad=8)
fig.savefig(os.path.join(FIG, 'fig10_final_candidates.png'))
plt.close(fig)
print('fig10 rebuilt (sorted by robustness group)')
print(top[['target_id', 'peptide_id', 'composite_target', 'n_weightings_selected']].to_string(index=False))
