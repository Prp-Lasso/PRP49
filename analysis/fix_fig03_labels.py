"""fig03 v2: fix the group labels.

The v1 bug: rows were listed bacterial-first but plotted with y = len-1-k, which puts
k=0 at the TOP - so the bacteria landed in the upper half while the labels said the upper
half was human. Both the labels AND the layout are now explicit, and human targets are
placed in the upper half because that is the claim the figure supports (lasso peptides do
reach human targets).
"""
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ROOT = r'D:\deepseek_harness\prp49'
FIG = os.path.join(ROOT, 'docs', 'figures')
MAIN, ACC, GREY = '#2b6cb0', '#dd6b20', '#718096'

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'figure.dpi': 200, 'savefig.bbox': 'tight'})

# human first => plotted at the top (y = n-1-k)
human = [
    ('RES-701-3', 'ETB receptor', 'A', 'cryo-EM 9KDF'),
    ('RES-701-1', 'ETB receptor', 'B', 'IC50 10 nM'),
    ('Anantin', 'NPR1', 'B', 'Kd 0.6 uM'),
    ('BI-32169', 'GCGR', 'B', 'cell cAMP assay'),
]
bact = [
    ('MccJ25', 'RNAP beta-prime', 'A', 'X-ray 6N60'),
    ('Capistruin', 'RNAP beta-prime', 'A', 'X-ray 6N61'),
    ('Lassomycin', 'ClpC1 (Mtb)', 'A', 'X-ray 8IBO'),
    ('Lariocidin', '30S ribosome', 'B', '2025 report'),
    ('Ubonodin', 'RNAP', 'B', 'in vitro',
     ),
    ('Klebsidin', 'RNAP', 'B', 'in vitro'),
]
rows = human + bact
n_h, n_b = len(human), len(bact)
gmap = {'A': ACC, 'B': MAIN}

fig, ax = plt.subplots(figsize=(8.2, 3.8))
for k, (pep, tg, g, note) in enumerate(rows):
    y = len(rows) - 1 - k
    ax.scatter(0, y, s=250, color=gmap[g], zorder=3)
    ax.text(0, y, g, ha='center', va='center', color='white', fontsize=8.5,
            fontweight='bold', zorder=4)
    ax.text(0.17, y, f'{pep}', va='center', fontsize=9, fontweight='bold')
    ax.text(0.62, y, f'\u2192  {tg}', va='center', fontsize=9)
    ax.text(1.32, y, note, va='center', fontsize=7.5, color=GREY, style='italic')

boundary = len(bact) - 0.5           # y between the last bacterial row and the first human row
ax.axhline(boundary, ls=':', color='#a0aec0', lw=1.2)
ax.text(-0.42, len(rows) - 1 - (n_h - 1) / 2, 'Human\ntargets', rotation=90,
        va='center', ha='center', fontsize=10, color=ACC, fontweight='bold')
ax.text(-0.42, (n_b - 1) / 2, 'Bacterial\ntargets', rotation=90,
        va='center', ha='center', fontsize=10, color=MAIN, fontweight='bold')

ax.set_xlim(-0.7, 2.15)
ax.set_ylim(-0.7, len(rows) - 0.3)
ax.set_xticks([]); ax.set_yticks([])
for sp in ax.spines.values():
    sp.set_visible(False)
ax.legend(handles=[Patch(color=ACC, label='grade A: co-crystal or cryo-EM structure'),
                   Patch(color=MAIN, label='grade B: binding or activity assay')],
          fontsize=8, loc='upper left', bbox_to_anchor=(0.0, -0.01),
          frameon=False, ncol=2)
ax.set_title('Characterised lasso peptide\u2013target pairs: 4 human targets, 6 bacterial targets',
             fontsize=10.5, pad=10)
fig.savefig(os.path.join(FIG, 'fig03_target_spectrum.png'))
plt.close(fig)
print('fig03 v2 written (human rows on top, labels match)')
print(f'  human rows y = {[len(rows)-1-k for k in range(n_h)]}')
print(f'  bacterial rows y = {[len(rows)-1-k for k in range(n_h, len(rows))]}')
print(f'  divider at y = {boundary}')
