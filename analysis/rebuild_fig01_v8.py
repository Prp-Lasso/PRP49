"""fig01 v8 - per-chain gradient from N-terminus to C-terminus.

v7 used one flat colour, which lost the sense of direction along the chain. Now every
panel is coloured by position along the chain:

  * bead colour   = colormap(t)          where t = 0 at the N-terminus, 1 at the C-terminus
  * segment colour= a linearGradient running from colormap(t_a) to colormap(t_b) along the
                    segment itself, so the gradient follows the chain rather than the page
  * the cyclic panel's closing bond runs from t=0.9 back to t=0, which is exactly the
    head-to-tail jump - the colour discontinuity marks the cyclisation
  * a colour bar under the figure states the reading direction (N -> C)

viridis is used because it is perceptually uniform and colour-blind safe.
"""
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.cm as cm
import numpy as np
from matplotlib.colors import to_hex

ROOT = r'D:\deepseek_harness\prp49'
FIG = os.path.join(ROOT, 'docs', 'figures')
WARN, GREY, INK = '#c53030', '#718096', '#1a202c'
FONT = 'Helvetica, Arial, sans-serif'
BEAD = 8.0


def cc(t):
    """colour at fractional chain position t (0 = N-terminus, 1 = C-terminus)"""
    return to_hex(cm.viridis(float(np.clip(t, 0.0, 1.0))))


def bead(x, y, t, label=None):
    col = cc(t)
    s = f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{BEAD}" fill="{col}"/>'
    if label is not None:
        s += (f'<text x="{x:.1f}" y="{y+3.2:.1f}" text-anchor="middle" font-size="8" '
              f'font-weight="bold" fill="white" stroke="#1a202c" stroke-width="2.2" '
              f'paint-order="stroke" stroke-linejoin="round">{label}</text>')
    return s


def segment(gid, p, q, ta, tb, width=3.0):
    """a bond drawn as a gradient along its own direction"""
    (x1, y1), (x2, y2) = p, q
    grad = (f'<linearGradient id="{gid}" gradientUnits="userSpaceOnUse" '
            f'x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}">'
            f'<stop offset="0" stop-color="{cc(ta)}"/>'
            f'<stop offset="1" stop-color="{cc(tb)}"/></linearGradient>')
    line = (f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
            f'stroke="url(#{gid})" stroke-width="{width}" stroke-linecap="round"/>')
    return grad, line


grads, body = [], []

# ---------------- A: linear peptide, 9 beads ----------------
ax_, ay_ = 42, 200
lin = [(ax_ + i * 24, ay_) for i in range(9)]
for i in range(8):
    g, l = segment(f'gl{i}', lin[i], lin[i + 1], i / 8, (i + 1) / 8)
    grads.append(g); body.append(l)
for i, (x, y) in enumerate(lin):
    body.append(bead(x, y, i / 8, 1 if i == 0 else (9 if i == 8 else None)))

# ---------------- B: cyclic peptide, 10 beads, closed ----------------
bx_, by_, br = 340, 200, 44
cyc = [(bx_ + br * np.cos(-np.pi / 2 + i * np.pi / 5),
        by_ + br * np.sin(-np.pi / 2 + i * np.pi / 5)) for i in range(10)]
for i in range(10):
    g, l = segment(f'gc{i}', cyc[i], cyc[(i + 1) % 10], i / 10, ((i + 1) % 10) / 10)
    grads.append(g); body.append(l)
for i, (x, y) in enumerate(cyc):
    body.append(bead(x, y, i / 10, 1 if i == 0 else (10 if i == 9 else None)))

# ---------------- C: lasso, 20 beads ----------------
CX, CY, R = 600, 200, 46
ring = [(CX + R * np.cos(-np.pi / 2 + i * np.pi / 4),
         CY + R * np.sin(-np.pi / 2 + i * np.pi / 4)) for i in range(8)]
tail = [(556, 132), (578, 108), (612, 100), (646, 106), (668, 128), (676, 160),
        (668, 190), (625, 202), (600, 210), (575, 218), (548, 246), (540, 278)]
NBEAD = 20
ts = [i / (NBEAD - 1) for i in range(NBEAD)]      # bead 1 -> t=0, bead 20 -> t=1

# assertions (topology unchanged from v7)
d_rt = min(np.hypot(tx - rx, ty - ry) for tx, ty in tail for rx, ry in ring)
inside = [i + 9 for i, (x, y) in enumerate(tail) if np.hypot(x - CX, y - CY) < R - 6]
gap89 = np.hypot(tail[0][0] - ring[7][0], tail[0][1] - ring[7][1])
cross = []
for (x1, y1), (x2, y2) in zip(tail[:-1], tail[1:]):
    f = lambda s: np.hypot(x1 + s * (x2 - x1) - CX, y1 + s * (y2 - y1) - CY) - R
    if f(0) * f(1) < 0:
        a, b = 0.0, 1.0
        for _ in range(40):
            m = (a + b) / 2
            a, b = (a, m) if f(a) * f(m) <= 0 else (m, b)
        cross.append((x1 + (a + b) / 2 * (x2 - x1), y1 + (a + b) / 2 * (y2 - y1)))
assert 20 < gap89 < 55 and d_rt > 14 and len(inside) >= 3 and len(cross) == 2
print(f'topology asserted: 8-9 bond {gap89:.0f}px | min sep {d_rt:.0f}px | '
      f'in-hole {inside} | crossings {len(cross)}')

# tail segments (lower z: ring will occlude at crossings)
for i in range(len(tail) - 1):
    g, l = segment(f'gt{i}', tail[i], tail[i + 1], ts[8 + i], ts[9 + i])
    grads.append(g); body.append(l)
# the 8->9 bond
g, l = segment('g89', ring[7], tail[0], ts[7], ts[8])
grads.append(g); body.append(l)
# ring segments on top
for i in range(8):
    g, l = segment(f'gr{i}', ring[i], ring[(i + 1) % 8], ts[i], ts[(i + 1) % 8])
    grads.append(g); body.append(l)
# beads: tail first then ring (ring on top = ring hides tail)
for i, (x, y) in enumerate(tail, start=9):
    body.append(bead(x, y, ts[i - 1], i if i in (9, 16, 18, 20) else None))
for i, (x, y) in enumerate(ring, start=1):
    body.append(bead(x, y, ts[i - 1], i if i in (1, 8) else None))

w, h = 900, 372
s = [f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}">
<rect width="{w}" height="{h}" fill="white"/>
<style>
  text {{ font-family: {FONT}; fill: {INK}; }}
  .t  {{ font-size: 15px; font-weight: bold; }}
  .xs {{ font-size: 10px; fill: #4a5568; }}
  .lbl {{ font-size: 12.5px; font-weight: bold; }}
</style>
<defs>
  <marker id="arw" markerWidth="10" markerHeight="10" refX="8" refY="5" orient="auto">
    <path d="M0,0 L10,5 L0,10 z" fill="{GREY}"/>
  </marker>
  {''.join(grads)}
</defs>
<text x="{w/2}" y="24" text-anchor="middle" class="t">Lasso peptide: colour runs along the chain from N-terminus to C-terminus</text>
''']
s += ['<text x="130" y="62" text-anchor="middle" class="lbl">Linear peptide</text>']
s += [g for g in body[:0]]
# assemble: linear block (first 8 grads/lines + 9 beads) is already in `body` in order
# easier: emit sections by re-walking the order we appended
idx = 0
# linear: 8 segments + 9 beads
lin_lines = body[0:8]; lin_beads = body[8:17]
s += lin_lines + lin_beads
s += ['<text x="130" y="248" text-anchor="middle" class="xs">flexible backbone</text>',
      '<text x="130" y="264" text-anchor="middle" class="xs">proteases cleave anywhere</text>']
# cyclic: 10 segments + 10 beads
cyc_lines = body[17:27]; cyc_beads = body[27:37]
s += ['<text x="340" y="62" text-anchor="middle" class="lbl">Cyclic peptide</text>']
s += cyc_lines + cyc_beads
s += ['<text x="340" y="282" text-anchor="middle" class="xs">closed head-to-tail (the closing bond jumps</text>',
      '<text x="340" y="296" text-anchor="middle" class="xs">from t=0.9 back to t=0 - visible as a colour step)</text>']
# lasso
s += [f'<text x="{CX+60}" y="62" text-anchor="middle" class="lbl">Lasso peptide</text>']
s += body[37:]
# isopeptide bond
x1, y1 = ring[0]; x8, y8 = ring[7]
s += [f'<path d="M {x1:.1f} {y1:.1f} A {R*0.62:.1f} {R*0.62} 0 0 1 {x8:.1f} {y8:.1f}" fill="none" '
      f'stroke="{WARN}" stroke-width="2.4" stroke-dasharray="5 3"/>',
      f'<text x="{CX-96}" y="{CY-46}" text-anchor="end" class="xs" fill="{WARN}">isopeptide bond</text>',
      f'<path d="M {CX-92} {CY-40} L {CX-18} {CY-40}" stroke="{WARN}" stroke-width="1" stroke-dasharray="2 2"/>']
# callouts
s += [f'<text x="716" y="118" text-anchor="start" class="xs">1. the chain leaves the ring at bead 9</text>',
      f'<path d="M 712 114 L {tail[1][0]+16} {tail[1][1]+4}" fill="none" stroke="{GREY}" '
      f'stroke-width="1.1" marker-end="url(#arw)"/>',
      f'<text x="716" y="150" text-anchor="start" class="xs">2. it turns back and comes down</text>',
      f'<path d="M 712 146 L {tail[4][0]+12} {tail[4][1]}" fill="none" stroke="{GREY}" '
      f'stroke-width="1.1" marker-end="url(#arw)"/>',
      f'<text x="716" y="182" text-anchor="start" class="xs">3. and threads through the ring</text>',
      f'<text x="716" y="196" text-anchor="start" class="xs">   (beads 16-18 sit in the hole;</text>',
      f'<text x="716" y="210" text-anchor="start" class="xs">   the ring hides the tail where</text>',
      f'<text x="716" y="224" text-anchor="start" class="xs">   they cross)</text>',
      f'<path d="M 712 190 L {tail[6][0]+14} {tail[6][1]}" fill="none" stroke="{GREY}" '
      f'stroke-width="1.1" marker-end="url(#arw)"/>']

# colour bar: N -> C
bar_x, bar_y, bar_w, bar_h = 250, 316, 400, 13
stops = ''.join(f'<stop offset="{i/20:.2f}" stop-color="{cc(i/20)}"/>' for i in range(21))
s += [f'<defs><linearGradient id="cbar" x1="{bar_x}" y1="0" x2="{bar_x+bar_w}" y2="0" '
      f'gradientUnits="userSpaceOnUse">{stops}</linearGradient></defs>',
      f'<rect x="{bar_x}" y="{bar_y}" width="{bar_w}" height="{bar_h}" rx="3" fill="url(#cbar)"/>',
      f'<text x="{bar_x-12}" y="{bar_y+11}" text-anchor="end" class="xs" font-weight="bold">N-terminus</text>',
      f'<text x="{bar_x+bar_w+12}" y="{bar_y+11}" text-anchor="start" class="xs" font-weight="bold">C-terminus</text>',
      f'<text x="{w/2}" y="{bar_y+34}" text-anchor="middle" class="xs">'
      f'Colour encodes position along the chain; each bond is drawn as its own gradient.</text>',
      f'<text x="{w/2}" y="{bar_y+48}" text-anchor="middle" class="xs">'
      f'The cyclic panel shows the closing bond jumping from t=0.9 back to t=0 - the head-to-tail step.</text>',
      '</svg>']
open(os.path.join(FIG, 'fig01_lasso_topology.svg'), 'w', encoding='utf-8').write('\n'.join(s))
print('fig01 v8 written')

import cairosvg
cairosvg.svg2png(url=os.path.join(FIG, 'fig01_lasso_topology.svg'),
                 write_to=os.path.join(FIG, 'fig01_lasso_topology_preview.png'),
                 scale=2.0, background_color='white')
print('rendered')
