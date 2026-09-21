"""Rebuild the three figures the user flagged.

fig01 - the old version drew a ring plus a detached line, so nothing said "one chain that
        closes on itself and threads its own tail through". Rebuilt as a SINGLE bead chain:
        beads 1-8 close the macrolactam ring, beads 9-16 are the C-terminal tail, drawn
        passing THROUGH the ring's hole (white halo puts the tail in front) and then
        bending away. The isopeptide bond is marked between bead 1 (N-terminus) and bead 6
        (the Asp/Glu side chain) - i.e. a side-chain-to-backbone bond, not head-to-tail.

fig03 - legend overlapped the last two rows; moved outside the axes.
fig07 - legend overlapped the first bar; moved above, and the right panel is reduced to the
        one comparison that is actually same-protocol (route A was only ever measured under
        family-grouped CV), with that stated on the panel.
"""
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

ROOT = r'D:\deepseek_harness\prp49'
FIG = os.path.join(ROOT, 'docs', 'figures')
MAIN, ACC, WARN, GREY, INK = '#2b6cb0', '#dd6b20', '#c53030', '#718096', '#1a202c'
FONT = 'Helvetica, Arial, sans-serif'

# =====================================================================
# fig01 - single bead chain forming a lasso
# =====================================================================
CX, CY, R = 700, 152, 52          # ring centre + radius
w, h = 900, 300
BEAD = 7.5

# ring beads 1..8 on a circle, starting at the top and going clockwise
ring = []
for i in range(8):
    th = -np.pi / 2 + i * (2 * np.pi / 8)
    ring.append((CX + R * np.cos(th), CY + R * np.sin(th)))

# tail beads 9..16: start inside the ring, thread up through the hole, then bend right.
# param a smooth path: from just below centre, up past the ring, then out to the right
t = np.linspace(0, 1, 8)
tail = []
for u in t:
    # quadratic-ish path
    x = CX - 14 + 96 * u ** 1.6
    y = CY + 34 - 118 * u + 22 * u ** 2
    tail.append((x, y))

s = f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}">
<rect width="{w}" height="{h}" fill="white"/>
<style>
  text {{ font-family: {FONT}; fill: {INK}; }}
  .t  {{ font-size: 15px; font-weight: bold; }}
  .xs {{ font-size: 10px; fill: #4a5568; }}
  .lbl {{ font-size: 12.5px; font-weight: bold; }}
  .num {{ font-size: 8.5px; fill: white; font-weight: bold; }}
</style>
<text x="{w/2}" y="24" text-anchor="middle" class="t">Lasso peptide: one chain closes a ring, then threads its own tail through it</text>
'''

# --- panel A: linear peptide (same bead style so the contrast is like-for-like)
s += '<text x="150" y="64" text-anchor="middle" class="lbl">Linear peptide</text>'
for i in range(9):
    x = 48 + i * 25.5
    s += f'<circle cx="{x:.1f}" cy="140" r="{BEAD}" fill="{MAIN}" opacity="0.9"/>'
    if i:
        s += f'<line x1="{x-25.5+BEAD:.1f}" y1="140" x2="{x-BEAD:.1f}" y2="140" stroke="{MAIN}" stroke-width="2.2"/>'
s += '<text x="150" y="182" text-anchor="middle" class="xs">flexible backbone</text>'
s += '<text x="150" y="198" text-anchor="middle" class="xs">proteases cleave anywhere</text>'

# --- panel B: cyclic peptide (beads closed head-to-tail)
s += '<text x="425" y="64" text-anchor="middle" class="lbl">Cyclic peptide</text>'
for i in range(10):
    th = -np.pi / 2 + i * (2 * np.pi / 10)
    x, y = 425 + 46 * np.cos(th), 152 + 46 * np.sin(th)
    s += f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{BEAD}" fill="{MAIN}" opacity="0.9"/>'
s += '<text x="425" y="222" text-anchor="middle" class="xs">closed head-to-tail</text>'
s += '<text x="425" y="238" text-anchor="middle" class="xs">rigid, but nothing threads it</text>'

# --- panel C: lasso (ring beads + tail beads, ONE chain)
s += f'<text x="{CX}" y="64" text-anchor="middle" class="lbl" fill="{ACC}">Lasso peptide</text>'

# white halo under the tail so it reads as passing IN FRONT of the ring
halo = 'M ' + ' L '.join(f'{x:.1f} {y:.1f}' for x, y in tail)
s += f'<path d="{halo}" fill="none" stroke="white" stroke-width="{BEAD*2.8:.0f}" stroke-linecap="round" stroke-linejoin="round"/>'
# ring beads drawn BEFORE the tail, so the tail's halo cuts them where it passes in front
for i, (x, y) in enumerate(ring, start=1):
    s += f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{BEAD}" fill="{ACC}" opacity="0.92"/>'
# connector lines inside the ring (the backbone)
for i in range(len(ring)):
    x1, y1 = ring[i]
    x2, y2 = ring[(i + 1) % len(ring)]
    s += (f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
          f'stroke="{ACC}" stroke-width="2.2" opacity="0.85"/>')
# tail beads on top
for i in range(len(tail) - 1):
    x1, y1 = tail[i]
    x2, y2 = tail[i + 1]
    s += f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{ACC}" stroke-width="2.2"/>'
for i, (x, y) in enumerate(tail, start=9):
    s += f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{BEAD}" fill="{ACC}"/>'
    s += f'<text x="{x:.1f}" y="{y+3:.1f}" text-anchor="middle" class="num">{i}</text>'

# number the ring beads too
for i, (x, y) in enumerate(ring, start=1):
    s += f'<text x="{x:.1f}" y="{y+3:.1f}" text-anchor="middle" class="num">{i}</text>'

# isopeptide bond: bead 1 (N-term alpha-amine) to bead 6 (Asp/Glu side chain)
x1, y1 = ring[0]
x6, y6 = ring[5]
s += (f'<path d="M {x1:.1f} {y1:.1f} Q {(x1+x6)/2-26:.1f} {(y1+y6)/2:.1f} {x6:.1f} {y6:.1f}" '
      f'fill="none" stroke="{WARN}" stroke-width="2.4" stroke-dasharray="5 3"/>')
s += (f'<text x="{(x1+x6)/2-58:.1f}" y="{(y1+y6)/2+4:.1f}" text-anchor="middle" class="xs" '
      f'fill="{WARN}">isopeptide bond</text>')
s += (f'<text x="{(x1+x6)/2-58:.1f}" y="{(y1+y6)/2+18:.1f}" text-anchor="middle" class="xs" '
      f'fill="{WARN}">(N-terminus &#8596; Asp/Glu side chain)</text>')

s += f'<text x="{CX+30}" y="228" text-anchor="middle" class="xs">beads 1-8: macrolactam ring &#183; beads 9-16: C-terminal tail</text>'
s += f'<text x="{CX+30}" y="246" text-anchor="middle" class="xs">the tail passes through the ring, then bends away</text>'
s += f'<text x="{w/2}" y="{h-10}" text-anchor="middle" class="xs">Trapped topology: proteases cannot reach the backbone, and the fold survives heat.</text>'
s += '</svg>'
open(os.path.join(FIG, 'fig01_lasso_topology.svg'), 'w', encoding='utf-8').write(s)
print('fig01 rebuilt as a single bead chain')

import cairosvg
cairosvg.svg2png(url=os.path.join(FIG, 'fig01_lasso_topology.svg'),
                 write_to=os.path.join(FIG, 'fig01_lasso_topology_preview.png'),
                 scale=2.0, background_color='white')

# =====================================================================
# fig03 - legend outside
# =====================================================================
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'figure.dpi': 200, 'savefig.bbox': 'tight'})
fig, ax = plt.subplots(figsize=(7.8, 3.4))
rows = [
    ('MccJ25', 'RNAP beta-prime', 'A', 'B'),
    ('Capistruin', 'RNAP beta-prime', 'A', 'B'),
    ('Lassomycin', 'ClpC1 (Mtb)', 'A', 'B'),
    ('Lariocidin', '30S ribosome', 'B', 'B'),
    ('Ubonodin', 'RNAP', 'B', 'B'),
    ('Klebsidin', 'RNAP', 'B', 'B'),
    ('RES-701-3', 'ETB receptor', 'A', 'H'),
    ('RES-701-1', 'ETB receptor', 'B', 'H'),
    ('Anantin', 'NPR1', 'B', 'H'),
    ('BI-32169', 'GCGR', 'B', 'H'),
]
gmap = {'A': ACC, 'B': MAIN}
for k, (pep, tg, g, grp) in enumerate(rows):
    y = len(rows) - 1 - k
    ax.scatter(0, y, s=230, color=gmap[g], zorder=3)
    ax.text(0, y, g, ha='center', va='center', color='white', fontsize=8,
            fontweight='bold', zorder=4)
    ax.text(0.16, y, f'{pep}  \u2192  {tg}', va='center', fontsize=8.5)
ax.axhline(5.5, ls=':', color='#a0aec0')
ax.text(-0.34, 7.5, 'Human\ntargets', rotation=90, va='center', ha='center',
        fontsize=9.5, color=ACC)
ax.text(-0.34, 2.5, 'Bacterial\ntargets', rotation=90, va='center', ha='center',
        fontsize=9.5, color=MAIN)
ax.set_xlim(-0.6, 1.75)
ax.set_ylim(-0.7, len(rows) - 0.3)
ax.set_xticks([]); ax.set_yticks([])
for sp in ax.spines.values():
    sp.set_visible(False)
ax.legend(handles=[Patch(color=ACC, label='grade A: co-crystal or cryo-EM'),
                   Patch(color=MAIN, label='grade B: binding assay')],
          fontsize=8, loc='upper left', bbox_to_anchor=(0.0, -0.02),
          frameon=False, ncol=2)
ax.set_title('Characterised lasso peptide\u2013target pairs: bacterial vs human targets',
             fontsize=10.5, pad=10)
fig.savefig(os.path.join(FIG, 'fig03_target_spectrum.png'))
plt.close(fig)
print('fig03 rebuilt (legend below the axes, no overlap)')

# =====================================================================
# fig07 - legend above, right panel reduced to the same-protocol comparison
# =====================================================================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.0, 3.4))
b = ax1.bar(['Reported\npositives', 'After family-homolog\naugmentation'], [352, 759],
            color=[WARN, MAIN], width=0.5)
for bb, v in zip(b, [352, 759]):
    ax1.text(bb.get_x() + bb.get_width() / 2, v + 14, str(v), ha='center', fontweight='bold')
ax1.set_ylabel('Training positives', fontsize=9.5)
ax1.set_ylim(0, 900)
ax1.set_title('Sample augmentation', fontsize=10)

x = np.arange(2)
ax2.bar(x, [0.8142, 0.7481], width=0.46, color=[MAIN, WARN])
for xi, v in zip(x, [0.8142, 0.7481]):
    ax2.text(xi, v + 0.012, f'{v:.4f}', ha='center', fontsize=9, fontweight='bold')
ax2.set_xticks(x)
ax2.set_xticklabels(['Family-grouped CV\nbaseline', 'Family-grouped CV\nwith propagation\n(route A)'],
                    fontsize=8.5)
ax2.set_ylim(0.6, 0.90)
ax2.set_ylabel('CV AUC', fontsize=9.5)
ax2.set_title('Same protocol, same folds: propagation hurts', fontsize=10)
ax2.annotate('', xy=(1, 0.7481), xytext=(0, 0.8142),
             arrowprops=dict(arrowstyle='->', color=GREY, lw=1.3, ls='--'))
ax2.text(0.5, 0.845, '\u22120.066', ha='center', fontsize=8.5, color=WARN, fontweight='bold')
fig.suptitle('Positive-sample scarcity and the family-homolog augmentation strategy',
             fontsize=10.5, y=1.04)
fig.savefig(os.path.join(FIG, 'fig07_family_augmentation.png'))
plt.close(fig)
print('fig07 rebuilt (no legend overlap; right panel is a same-protocol comparison)')
