"""fig01 v4: keep v3's (working) threading idea, fix the layout collisions.

v3 was right about the geometry - tail enters the ring, crosses between beads 1 and 2,
exits - but the tail climbed to y=52 and ran into the panel title, and the two bottom
captions overlapped because a single <text> mixed tspans.

v4 changes only the layout:
  * ring centre moved down (y=180) and the tail now flattens out at y=92, leaving 28px
    below the panel title
  * the crossing callout moves to the empty lower-right instead of sitting on the tail
  * bottom captions are two separate <text> lines with no tspan mixing
"""
import os

import numpy as np

ROOT = r'D:\deepseek_harness\prp49'
FIG = os.path.join(ROOT, 'docs', 'figures')
MAIN, WARN, GREY, INK = '#2b6cb0', '#c53030', '#718096', '#1a202c'
RING_C, TAIL_C = '#dd6b20', '#9c4221'
FONT = 'Helvetica, Arial, sans-serif'

CX, CY, R = 680, 180, 56
BEAD = 7.5
ring = [(CX + R * np.cos(-np.pi / 2 + i * np.pi / 4),
         CY + R * np.sin(-np.pi / 2 + i * np.pi / 4)) for i in range(8)]

tail = [(CX - 4, CY + 30),
        (CX + 2, CY + 10),
        (CX + 10, CY - 15),
        (CX + 21, CY - 52),      # crossing, between ring beads 1 and 2
        (CX + 42, CY - 72),
        (CX + 72, CY - 84),
        (CX + 108, CY - 88),
        (CX + 144, CY - 88)]

mind = min(np.hypot(tx - rx, ty - ry) for tx, ty in tail for rx, ry in ring)
top_y = min(y for _, y in tail)
print(f'min tail-ring distance {mind:.1f}px | tail apex y={top_y:.0f} (title at y=62)')
assert mind > 2 * BEAD + 2, 'tail beads collide with ring beads'
assert top_y > 80, 'tail would run into the panel title'

w, h = 900, 300
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

# A linear
s += '<text x="150" y="62" text-anchor="middle" class="lbl">Linear peptide</text>'
for i in range(9):
    x = 48 + i * 25.5
    if i:
        s += f'<line x1="{x-25.5+BEAD:.1f}" y1="170" x2="{x-BEAD:.1f}" y2="170" stroke="{MAIN}" stroke-width="2.4"/>'
for i in range(9):
    x = 48 + i * 25.5
    s += f'<circle cx="{x:.1f}" cy="170" r="{BEAD}" fill="{MAIN}"/>'
s += '<text x="150" y="215" text-anchor="middle" class="xs">flexible backbone</text>'
s += '<text x="150" y="231" text-anchor="middle" class="xs">proteases cleave anywhere</text>'

# B cyclic (with connecting lines)
s += '<text x="405" y="62" text-anchor="middle" class="lbl">Cyclic peptide</text>'
cr = 46
cyc = [(405 + cr * np.cos(-np.pi / 2 + i * np.pi / 5),
        180 + cr * np.sin(-np.pi / 2 + i * np.pi / 5)) for i in range(10)]
for i in range(10):
    x1, y1 = cyc[i]; x2, y2 = cyc[(i + 1) % 10]
    s += f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{MAIN}" stroke-width="2.4"/>'
for x, y in cyc:
    s += f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{BEAD}" fill="{MAIN}"/>'
s += '<text x="405" y="250" text-anchor="middle" class="xs">closed head-to-tail</text>'
s += '<text x="405" y="266" text-anchor="middle" class="xs">rigid, but nothing threads it</text>'

# C lasso
s += f'<text x="{CX+50}" y="62" text-anchor="middle" class="lbl" fill="{TAIL_C}">Lasso peptide</text>'
halo = 'M ' + ' L '.join(f'{x:.1f} {y:.1f}' for x, y in tail)
s += (f'<path d="{halo}" fill="none" stroke="white" stroke-width="{BEAD*3:.0f}" '
      f'stroke-linecap="round" stroke-linejoin="round"/>')
for i in range(8):
    x1, y1 = ring[i]; x2, y2 = ring[(i + 1) % 8]
    s += f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{RING_C}" stroke-width="2.4"/>'
for x, y in ring:
    s += f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{BEAD}" fill="{RING_C}"/>'
for i in range(len(tail) - 1):
    x1, y1 = tail[i]; x2, y2 = tail[i + 1]
    s += f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{TAIL_C}" stroke-width="2.8"/>'
for i, (x, y) in enumerate(tail, start=9):
    s += f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{BEAD}" fill="{TAIL_C}"/>'
    s += f'<text x="{x:.1f}" y="{y+3:.1f}" text-anchor="middle" class="num">{i}</text>'
for i, (x, y) in enumerate(ring, start=1):
    s += f'<text x="{x:.1f}" y="{y+3:.1f}" text-anchor="middle" class="num">{i}</text>'

# isopeptide bond between bead 1 and bead 6
x1, y1 = ring[0]; x6, y6 = ring[5]
s += (f'<path d="M {x1:.1f} {y1:.1f} Q {(x1+x6)/2-28:.1f} {(y1+y6)/2:.1f} {x6:.1f} {y6:.1f}" '
      f'fill="none" stroke="{WARN}" stroke-width="2.4" stroke-dasharray="5 3"/>')
s += (f'<text x="{x1-84:.1f}" y="{y1-4:.1f}" text-anchor="middle" class="xs" fill="{WARN}">isopeptide bond</text>')
s += (f'<text x="{x6-96:.1f}" y="{y6+20:.1f}" text-anchor="middle" class="xs" fill="{WARN}">N-term to Asp/Glu side chain</text>')

# crossing callout, placed in the clear area to the lower right of the ring
tx, ty = tail[3]
s += (f'<circle cx="{tx:.1f}" cy="{ty:.1f}" r="16" fill="none" stroke="{GREY}" '
      f'stroke-width="1.4" stroke-dasharray="4 3"/>')
s += f'<path d="M {tx+14:.1f} {ty+9:.1f} L {tx+40:.1f} {ty+48:.1f}" stroke="{GREY}" stroke-width="1.2"/>'
s += f'<text x="{tx+44:.1f}" y="{ty+58:.1f}" class="xs">tail passes through the ring here</text>'

# captions: separate text elements, no tspan mixing
s += f'<text x="{CX+50}" y="252" text-anchor="middle" class="xs">beads 1-8: macrolactam ring (orange)</text>'
s += f'<text x="{CX+50}" y="268" text-anchor="middle" class="xs">beads 9-16: C-terminal tail (dark) - it enters, crosses, and exits</text>'
s += f'<text x="{w/2}" y="{h-8}" text-anchor="middle" class="xs">Trapped topology: proteases cannot reach the backbone, and the fold survives heat.</text>'
s += '</svg>'
open(os.path.join(FIG, 'fig01_lasso_topology.svg'), 'w', encoding='utf-8').write(s)
print('fig01 v4 written')

import cairosvg
cairosvg.svg2png(url=os.path.join(FIG, 'fig01_lasso_topology.svg'),
                 write_to=os.path.join(FIG, 'fig01_lasso_topology_preview.png'),
                 scale=2.0, background_color='white')
print('rendered')
