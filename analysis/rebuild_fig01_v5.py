"""fig01 v5 - the lasso geometry, done as a real lasso.

What was wrong in v4: the tail rose out of the ring's upper edge, which reads as "a ring
with a tail attached" - a balloon on a string. A lasso (and a lasso peptide) is different:
the rope's end is threaded BACK THROUGH the loop the rope itself forms. So the tail must
cross the ring's hole diagonally, with a visible segment on the entry side and another on
the exit side, and the two crossings must sit on opposite arcs of the ring.

Geometry used here (computed, then asserted):
  ring: centre (700,160), R=48, beads 1-8 every 45 deg starting at the top
  tail: a straight run through the ring centre along the 112.5deg -> -67.5deg chord,
        so it crosses the circumference halfway between beads 5-6 (lower left) and
        halfway between beads 1-2 (upper right); after exiting it bends to horizontal so
        it stays clear of the panel title.
  the lower crossing is where the tail ENTERS, the upper where it EXITS - that is the
  threading that makes it a lasso rather than a decoration.
"""
import os

import numpy as np

ROOT = r'D:\deepseek_harness\prp49'
FIG = os.path.join(ROOT, 'docs', 'figures')
MAIN, WARN, GREY, INK = '#2b6cb0', '#c53030', '#718096', '#1a202c'
RING_C, TAIL_C = '#dd6b20', '#9c4221'
FONT = 'Helvetica, Arial, sans-serif'
BEAD = 7.5

CX, CY, R = 700, 160, 48
ring = [(CX + R * np.cos(-np.pi / 2 + i * np.pi / 4),
         CY + R * np.sin(-np.pi / 2 + i * np.pi / 4)) for i in range(8)]

# tail: enters lower-left, passes through the centre, exits upper-right, then flattens
tail = [(660, 252),
        (672, 222),
        (686, 191),
        (700, 160),      # ring centre - the tail is inside the hole here
        (713, 130),
        (729, 108),      # just outside the ring on the upper right
        (764, 103),
        (800, 101)]

# assertions: geometry must actually be a threading
d_ring_tail = min(np.hypot(tx - rx, ty - ry) for tx, ty in tail for rx, ry in ring)
in_ring = [abs(np.hypot(x - CX, y - CY)) < R for x, y in tail]
print(f'min tail-ring bead distance: {d_ring_tail:.1f}px (need > {2*BEAD})')
print('tail beads inside the ring hole:', [i + 9 for i, v in enumerate(in_ring) if v])
print('tail apex y:', min(y for _, y in tail), '(panel title at y=62)')
assert d_ring_tail > 2 * BEAD + 1, 'tail beads collide with ring beads'
assert sum(in_ring) >= 2, 'the tail must be visibly inside the hole'
assert min(y for _, y in tail) > 80, 'tail would hit the panel title'

# where does the tail cross the circumference?
prev = None
crossings = []
for (x1, y1), (x2, y2) in zip(tail[:-1], tail[1:]):
    f = lambda s: np.hypot(x1 + s * (x2 - x1) - CX, y1 + s * (y2 - y1) - CY) - R
    a, b = 0.0, 1.0
    if f(a) * f(b) < 0:
        for _ in range(40):
            m = (a + b) / 2
            if f(a) * f(m) <= 0:
                b = m
            else:
                a = m
        s = (a + b) / 2
        crossings.append((x1 + s * (x2 - x1), y1 + s * (y2 - y1)))
print('circumference crossings at:', [(round(x), round(y)) for x, y in crossings])
assert len(crossings) == 2, 'the tail must cross the ring exactly twice'

w, h = 900, 322
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
<text x="{w/2}" y="24" text-anchor="middle" class="t">Lasso peptide: the chain's own tail is threaded back through its ring</text>
'''

# --- A linear
s += '<text x="150" y="62" text-anchor="middle" class="lbl">Linear peptide</text>'
for i in range(9):
    x = 48 + i * 25.5
    if i:
        s += f'<line x1="{x-25.5+BEAD:.1f}" y1="180" x2="{x-BEAD:.1f}" y2="180" stroke="{MAIN}" stroke-width="2.4"/>'
for i in range(9):
    x = 48 + i * 25.5
    s += f'<circle cx="{x:.1f}" cy="180" r="{BEAD}" fill="{MAIN}"/>'
s += '<text x="150" y="225" text-anchor="middle" class="xs">flexible backbone</text>'
s += '<text x="150" y="241" text-anchor="middle" class="xs">proteases cleave anywhere</text>'

# --- B cyclic
s += '<text x="405" y="62" text-anchor="middle" class="lbl">Cyclic peptide</text>'
cr = 46
cyc = [(405 + cr * np.cos(-np.pi / 2 + i * np.pi / 5),
        180 + cr * np.sin(-np.pi / 2 + i * np.pi / 5)) for i in range(10)]
for i in range(10):
    x1, y1 = cyc[i]; x2, y2 = cyc[(i + 1) % 10]
    s += f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{MAIN}" stroke-width="2.4"/>'
for x, y in cyc:
    s += f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{BEAD}" fill="{MAIN}"/>'
s += '<text x="405" y="256" text-anchor="middle" class="xs">closed head-to-tail</text>'
s += '<text x="405" y="272" text-anchor="middle" class="xs">rigid, but nothing threads it</text>'

# --- C lasso
s += f'<text x="{CX}" y="62" text-anchor="middle" class="lbl" fill="{TAIL_C}">Lasso peptide</text>'
# white halo first: the tail passes in front of the ring
halo = 'M ' + ' L '.join(f'{x} {y}' for x, y in tail)
s += (f'<path d="{halo}" fill="none" stroke="white" stroke-width="{BEAD*3:.0f}" '
      f'stroke-linecap="round" stroke-linejoin="round"/>')
# ring
for i in range(8):
    x1, y1 = ring[i]; x2, y2 = ring[(i + 1) % 8]
    s += f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{RING_C}" stroke-width="2.4"/>'
for x, y in ring:
    s += f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{BEAD}" fill="{RING_C}"/>'
# tail on top
for i in range(len(tail) - 1):
    x1, y1 = tail[i]; x2, y2 = tail[i + 1]
    s += f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{TAIL_C}" stroke-width="2.9"/>'
for i, (x, y) in enumerate(tail, start=9):
    s += f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{BEAD}" fill="{TAIL_C}"/>'
    s += f'<text x="{x:.1f}" y="{y+3:.1f}" text-anchor="middle" class="num">{i}</text>'
for i, (x, y) in enumerate(ring, start=1):
    s += f'<text x="{x:.1f}" y="{y+3:.1f}" text-anchor="middle" class="num">{i}</text>'

# isopeptide bond: bead 1 (N-term) to bead 6 (Asp/Glu side chain)
x1, y1 = ring[0]; x6, y6 = ring[5]
s += (f'<path d="M {x1:.1f} {y1:.1f} Q {(x1+x6)/2-30:.1f} {(y1+y6)/2:.1f} {x6:.1f} {y6:.1f}" '
      f'fill="none" stroke="{WARN}" stroke-width="2.4" stroke-dasharray="5 3"/>')
s += (f'<text x="{x1-64:.1f}" y="{y1+2:.1f}" text-anchor="middle" class="xs" fill="{WARN}">isopeptide</text>'
      f'<text x="{x1-64:.1f}" y="{y1+16:.1f}" text-anchor="middle" class="xs" fill="{WARN}">bond</text>')
s += (f'<text x="{x6-118:.1f}" y="{y6-6:.1f}" text-anchor="middle" class="xs" fill="{WARN}">N-terminus to</text>'
      f'<text x="{x6-118:.1f}" y="{y6+8:.1f}" text-anchor="middle" class="xs" fill="{WARN}">Asp/Glu side chain</text>')

# callouts on the two crossings: ENTER and EXIT
ex, ey = crossings[0]
xx, xy = crossings[1]
s += f'<circle cx="{ex:.1f}" cy="{ey:.1f}" r="13" fill="none" stroke="{GREY}" stroke-width="1.4" stroke-dasharray="4 3"/>'
s += (f'<path d="M {ex-15:.1f} {ey+14:.1f} L {ex-58:.1f} {ey+52:.1f}" stroke="{GREY}" '
      f'stroke-width="1.2" marker-end="url(#arw)"/>')
s += f'<text x="{ex-62:.1f}" y="{ey+52:.1f}" text-anchor="end" class="xs">tail enters here</text>'
s += f'<circle cx="{xx:.1f}" cy="{xy:.1f}" r="13" fill="none" stroke="{GREY}" stroke-width="1.4" stroke-dasharray="4 3"/>'
s += (f'<path d="M {xx+13:.1f} {xy-9:.1f} L {xx+58:.1f} {xy-34:.1f}" stroke="{GREY}" '
      f'stroke-width="1.2" marker-end="url(#arw)"/>')
s += f'<text x="{xx+62:.1f}" y="{xy-38:.1f}" class="xs">&#8230;and exits here</text>'

s += f'<text x="{CX}" y="268" text-anchor="middle" class="xs">beads 1-8: macrolactam ring &#183; beads 9-16: C-terminal tail (dark)</text>'
s += f'<text x="{CX}" y="286" text-anchor="middle" class="xs">like a lasso: the rope\'s end is threaded back through the loop it forms</text>'
s += f'<text x="{CX}" y="304" text-anchor="middle" class="xs">once pulled tight the topology is trapped - proteases cannot reach the backbone</text>'
s += '</svg>'
open(os.path.join(FIG, 'fig01_lasso_topology.svg'), 'w', encoding='utf-8').write(s)
print('fig01 v5 written')

import cairosvg
cairosvg.svg2png(url=os.path.join(FIG, 'fig01_lasso_topology.svg'),
                 write_to=os.path.join(FIG, 'fig01_lasso_topology_preview.png'),
                 scale=2.0, background_color='white')
print('rendered')
