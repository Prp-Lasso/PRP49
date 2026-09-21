"""fig01 v7 - same topology as v6 (which was correct), fixed layout:

  * the ring moves left (centre x 700 -> 600) so the callouts fit inside the canvas
    instead of running off the right edge
  * the tail's excursion is much smaller: it steps out above the ring and comes straight
    back down through the hole, instead of sweeping a second big circle that read as a
    "second ring"
  * the three in-hole beads sit well inside the hole radius so they no longer crowd the
    ring beads
  * everything is asserted again: 8-9 bonded, >=3 beads in the hole, exactly 2 crossings,
    minimum bead separation, apex clear of the title
"""
import os

import numpy as np

ROOT = r'D:\deepseek_harness\prp49'
FIG = os.path.join(ROOT, 'docs', 'figures')
CHAIN = '#dd6b20'
WARN, GREY, INK = '#c53030', '#718096', '#1a202c'
MAIN = '#2b6cb0'
FONT = 'Helvetica, Arial, sans-serif'
BEAD = 7.5

CX, CY, R = 600, 200, 46
ring = [(CX + R * np.cos(-np.pi / 2 + i * np.pi / 4),
         CY + R * np.sin(-np.pi / 2 + i * np.pi / 4)) for i in range(8)]

tail = [(556, 132),    # 9  - steps out above bead 8
        (578, 108),    # 10
        (612, 100),    # 11 - apex
        (646, 106),    # 12
        (668, 128),    # 13 - comes back down on the right
        (676, 160),    # 14 - outside the ring
        (668, 190),    # 15 - approaching the rim
        (625, 202),    # 16 - inside the hole
        (600, 210),    # 17 - inside the hole
        (575, 218),    # 18 - inside the hole
        (548, 246),    # 19 - exited
        (540, 278)]    # 20 - C-terminus

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

print(f'8->9 bond length          : {gap89:.1f}px')
print(f'min tail/ring separation  : {d_rt:.1f}px')
print(f'beads inside the hole     : {inside}  (hole radius {R}, safe inner zone {R-6})')
print(f'crossings                 : {[(round(x), round(y)) for x, y in cross]}')
print(f'tail apex y               : {min(y for _, y in tail):.0f}')

assert 20 < gap89 < 55, 'beads 8 and 9 must be bonded'
assert d_rt > 14, f'tail crowds ring beads ({d_rt:.1f}px)'
assert len(inside) >= 3, 'tail must pass inside the hole'
assert len(cross) == 2, 'tail must cross the rim exactly twice'
assert min(y for _, y in tail) > 80, 'tail hits the panel title'

w, h = 900, 356
s = f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}">
<rect width="{w}" height="{h}" fill="white"/>
<style>
  text {{ font-family: {FONT}; fill: {INK}; }}
  .t  {{ font-size: 15px; font-weight: bold; }}
  .xs {{ font-size: 10px; fill: #4a5568; }}
  .lbl {{ font-size: 12.5px; font-weight: bold; }}
  .num {{ font-size: 8.5px; fill: white; font-weight: bold; }}
</style>
<defs>
  <marker id="arw" markerWidth="10" markerHeight="10" refX="8" refY="5" orient="auto">
    <path d="M0,0 L10,5 L0,10 z" fill="{GREY}"/>
  </marker>
</defs>
<text x="{w/2}" y="24" text-anchor="middle" class="t">Lasso peptide: the tail leaves the ring, turns back, and threads through it</text>
'''

# A: linear
s += '<text x="130" y="62" text-anchor="middle" class="lbl">Linear peptide</text>'
for i in range(8):
    x = 42 + i * 24
    if i:
        s += f'<line x1="{x-24+BEAD:.1f}" y1="200" x2="{x-BEAD:.1f}" y2="200" stroke="{MAIN}" stroke-width="2.4"/>'
for i in range(8):
    x = 42 + i * 24
    s += f'<circle cx="{x:.1f}" cy="200" r="{BEAD}" fill="{MAIN}"/>'
s += '<text x="130" y="245" text-anchor="middle" class="xs">flexible backbone</text>'
s += '<text x="130" y="261" text-anchor="middle" class="xs">proteases cleave anywhere</text>'

# B: cyclic
s += '<text x="340" y="62" text-anchor="middle" class="lbl">Cyclic peptide</text>'
cr = 44
cyc = [(340 + cr * np.cos(-np.pi / 2 + i * np.pi / 5),
        200 + cr * np.sin(-np.pi / 2 + i * np.pi / 5)) for i in range(10)]
for i in range(10):
    x1, y1 = cyc[i]; x2, y2 = cyc[(i + 1) % 10]
    s += f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{MAIN}" stroke-width="2.4"/>'
for x, y in cyc:
    s += f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{BEAD}" fill="{MAIN}"/>'
s += '<text x="340" y="278" text-anchor="middle" class="xs">closed head-to-tail</text>'
s += '<text x="340" y="294" text-anchor="middle" class="xs">rigid, but nothing threads it</text>'

# C: lasso
s += f'<text x="{CX+60}" y="62" text-anchor="middle" class="lbl" fill="{CHAIN}">Lasso peptide</text>'
# tail first (ring will occlude it at the crossings)
for i in range(len(tail) - 1):
    x1, y1 = tail[i]; x2, y2 = tail[i + 1]
    s += (f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{CHAIN}" '
          f'stroke-width="3" stroke-linecap="round"/>')
# the 8->9 backbone bond
s += (f'<line x1="{ring[7][0]:.1f}" y1="{ring[7][1]:.1f}" x2="{tail[0][0]}" y2="{tail[0][1]}" '
      f'stroke="{CHAIN}" stroke-width="3" stroke-linecap="round"/>')
# ring on top
for i in range(8):
    x1, y1 = ring[i]; x2, y2 = ring[(i + 1) % 8]
    s += (f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
          f'stroke="{CHAIN}" stroke-width="3"/>')
# beads: tail then ring (ring on top = ring hides tail)
for i, (x, y) in enumerate(tail, start=9):
    s += f'<circle cx="{x}" cy="{y}" r="{BEAD}" fill="{CHAIN}"/>'
    s += f'<text x="{x}" y="{y+3}" text-anchor="middle" class="num">{i}</text>'
for i, (x, y) in enumerate(ring, start=1):
    s += f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{BEAD}" fill="{CHAIN}"/>'
    s += f'<text x="{x:.1f}" y="{y+3:.1f}" text-anchor="middle" class="num">{i}</text>'

# isopeptide bond: bead 1 (N-term) to bead 8 (Asp/Glu side chain)
x1, y1 = ring[0]; x8, y8 = ring[7]
s += (f'<path d="M {x1:.1f} {y1:.1f} A {R*0.62:.1f} {R*0.62} 0 0 1 {x8:.1f} {y8:.1f}" '
      f'fill="none" stroke="{WARN}" stroke-width="2.4" stroke-dasharray="5 3"/>')
s += (f'<text x="{CX-96}" y="{CY-46}" text-anchor="end" class="xs" fill="{WARN}">isopeptide bond</text>')
s += (f'<path d="M {CX-92} {CY-40} L {CX-18} {CY-40}" stroke="{WARN}" stroke-width="1" stroke-dasharray="2 2"/>')

# callouts on the right, all inside the canvas
s += f'<text x="716" y="118" text-anchor="start" class="xs">1. the chain leaves the ring at bead 9</text>'
s += (f'<path d="M 712 114 L {tail[1][0]+16} {tail[1][1]+4}" fill="none" stroke="{GREY}" '
      f'stroke-width="1.1" marker-end="url(#arw)"/>')
s += f'<text x="716" y="150" text-anchor="start" class="xs">2. it turns back and comes down</text>'
s += (f'<path d="M 712 146 L {tail[4][0]+12} {tail[4][1]}" fill="none" stroke="{GREY}" '
      f'stroke-width="1.1" marker-end="url(#arw)"/>')
s += f'<text x="716" y="182" text-anchor="start" class="xs">3. and threads through the ring</text>'
s += f'<text x="716" y="196" text-anchor="start" class="xs">   (beads 16-18 sit in the hole;</text>'
s += f'<text x="716" y="210" text-anchor="start" class="xs">   the ring hides the tail where</text>'
s += f'<text x="716" y="224" text-anchor="start" class="xs">   they cross)</text>'
s += (f'<path d="M 712 190 L {tail[6][0]+14} {tail[6][1]}" fill="none" stroke="{GREY}" '
      f'stroke-width="1.1" marker-end="url(#arw)"/>')

s += f'<text x="{w/2}" y="316" text-anchor="middle" class="xs">One chain, one colour: beads 1-8 close the macrolactam ring through an isopeptide bond between the N-terminus</text>'
s += f'<text x="{w/2}" y="330" text-anchor="middle" class="xs">and an Asp/Glu side chain. Beads 9-20 are that same chain continuing: pulled out of the ring, turned back,</text>'
s += f'<text x="{w/2}" y="344" text-anchor="middle" class="xs">and threaded through its own loop - the trapped topology that makes lasso peptides protease-resistant.</text>'
s += '</svg>'
open(os.path.join(FIG, 'fig01_lasso_topology.svg'), 'w', encoding='utf-8').write(s)
print('fig01 v7 written')

import cairosvg
cairosvg.svg2png(url=os.path.join(FIG, 'fig01_lasso_topology.svg'),
                 write_to=os.path.join(FIG, 'fig01_lasso_topology_preview.png'),
                 scale=2.0, background_color='white')
print('rendered')
