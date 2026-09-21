"""fig01 v11 - make both crossings visible AND opposite.

Two problems found by inspecting v10:

  1. OCCLUSION WAS UNIFORM. v10 put the whole tail behind the ring, so the pair read as
     "tail lying against the ring". The user needs one crossing with the ring in front and
     the other with the tail in front - that is what makes it interlinked.

  2. THE CROSSINGS SAT ON RING BEADS. The crossing points were 4-5 px from the nearest ring
     bead, so the bead (a filled disc) covered the crossing and no stroke-on-stroke
     occlusion was visible at all. The crossings must land midway BETWEEN beads, where
     ring stroke meets tail stroke.

Fix: compute the tail geometrically. The crossings are placed at the angular midpoints
between ring beads (-22.5 deg and its antipode 157.5 deg), so the tail runs straight
through the hole from one midpoint to the other, and every crossing is >=17 px from any
bead centroid. Layering is then split: tail beads 9-15 (and the 15->16 link) go UNDER the
ring; beads 16-20 (and the 18->19 link) go OVER it.
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
    return to_hex(cm.viridis(float(np.clip(t, 0.0, 1.0)) * 0.88))


def bead(x, y, t, label=None):
    s = f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{BEAD}" fill="{cc(t)}"/>'
    if label is not None:
        s += (f'<text x="{x:.1f}" y="{y+3.2:.1f}" text-anchor="middle" font-size="8.5" '
              f'font-weight="bold" fill="#ffffff">{label}</text>')
    return s


def segment(gid, p, q, ta, tb, width=3.0):
    (x1, y1), (x2, y2) = p, q
    grad = (f'<linearGradient id="{gid}" gradientUnits="userSpaceOnUse" '
            f'x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}">'
            f'<stop offset="0" stop-color="{cc(ta)}"/>'
            f'<stop offset="1" stop-color="{cc(tb)}"/></linearGradient>')
    line = (f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
            f'stroke="url(#{gid})" stroke-width="{width}" stroke-linecap="round"/>')
    return grad, line


grads = []

# ---------------- A: linear ----------------
A = []
lin = [(42 + i * 24, 200) for i in range(9)]
for i in range(8):
    g, l = segment(f'gl{i}', lin[i], lin[i + 1], i / 8, (i + 1) / 8)
    grads.append(g); A.append(l)
for i, (x, y) in enumerate(lin):
    A.append(bead(x, y, i / 8, 1 if i == 0 else (9 if i == 8 else None)))

# ---------------- B: cyclic ----------------
B = []
cyc = [(340 + 44 * np.cos(-np.pi / 2 + i * np.pi / 5),
        200 + 44 * np.sin(-np.pi / 2 + i * np.pi / 5)) for i in range(10)]
for i in range(10):
    g, l = segment(f'gc{i}', cyc[i], cyc[(i + 1) % 10], i / 10, ((i + 1) % 10) / 10)
    grads.append(g); B.append(l)
for i, (x, y) in enumerate(cyc):
    B.append(bead(x, y, i / 10, 1 if i == 0 else (10 if i == 9 else None)))

# ---------------- C: lasso, tail computed ----------------
CX, CY, R = 600, 200, 46
ring = [(CX + R * np.cos(-np.pi / 2 + i * np.pi / 4),
         CY + R * np.sin(-np.pi / 2 + i * np.pi / 4)) for i in range(8)]

ANG1 = np.deg2rad(-22.5)          # midway between ring beads 2 and 3
ANG2 = ANG1 + np.pi               # antipode: midway between beads 6 and 7
X1 = np.array([CX + R * np.cos(ANG1), CY + R * np.sin(ANG1)])
X2 = np.array([CX + R * np.cos(ANG2), CY + R * np.sin(ANG2)])
u = (X2 - X1) / np.linalg.norm(X2 - X1)      # unit vector along the threading chord
chord = np.linalg.norm(X2 - X1)

b15 = X1 - u * 14                  # outside, approaching the first crossing
b16 = X1 + u * (chord * 0.30)      # inside the hole
b17 = X1 + u * (chord * 0.55)
b18 = X1 + u * (chord * 0.80)
b19 = X2 + u * 14                  # outside, past the second crossing
b20 = X2 + u * 40

tail = [(556, 132), (580, 106), (616, 98), (652, 106), (676, 128), (672, 158),
        tuple(b15), tuple(b16), tuple(b17), tuple(b18), tuple(b19), tuple(b20)]
ts = [i / 19 for i in range(20)]

# ---- assertions ----
d_rt = min(np.hypot(tx - rx, ty - ry) for tx, ty in tail for rx, ry in ring)
inside = [i + 9 for i, (x, y) in enumerate(tail) if np.hypot(x - CX, y - CY) < R - 6]
gap89 = np.hypot(tail[0][0] - ring[7][0], tail[0][1] - ring[7][1])
cross, own = [], []
for k, ((x1, y1), (x2, y2)) in enumerate(zip(tail[:-1], tail[1:])):
    f = lambda s: np.hypot(x1 + s * (x2 - x1) - CX, y1 + s * (y2 - y1) - CY) - R
    if f(0) * f(1) < 0:
        a, b = 0.0, 1.0
        for _ in range(40):
            m = (a + b) / 2
            a, b = (a, m) if f(a) * f(m) <= 0 else (m, b)
        cross.append((x1 + (a + b) / 2 * (x2 - x1), y1 + (a + b) / 2 * (y2 - y1)))
        own.append(k)                      # which tail link carries this crossing
cross_to_bead = [min(np.hypot(x - rx, y - ry) for rx, ry in ring) for x, y in cross]

print(f'8-9 bond {gap89:.1f}px | min tail/ring bead sep {d_rt:.1f}px | in-hole {inside}')
print(f'crossings: {[(round(x), round(y)) for x, y in cross]}')
print(f'  on tail links: {own}  (need one < 7 and one >= 7)')
print(f'  distance to nearest ring bead: {[round(d, 1) for d in cross_to_bead]}  (need > 16)')

assert 20 < gap89 < 55, 'beads 8 and 9 must be bonded'
assert d_rt > 14, 'tail crowds ring beads'
assert len(inside) >= 3, 'tail must pass through the hole'
assert len(cross) == 2, 'exactly two crossings'
assert min(cross_to_bead) > 16, 'crossings must sit between beads, not on them'
assert own[0] < own[1], 'crossings must be ordered along the chain'
# the layering split happens between the two crossings
SPLIT = own[0] + 1
print(f'  layering split after tail link {own[0]} (link {own[0]} behind, link {own[1]} in front)')

# ---- layer 1: tail behind the ring: bond 8-9, links 0..SPLIT-1, beads 9..(9+SPLIT-1)
C_behind = []
g, l = segment('g89', ring[7], tail[0], ts[7], ts[8]); grads.append(g); C_behind.append(l)
for i in range(SPLIT):
    g, l = segment(f'gtA{i}', tail[i], tail[i + 1], ts[8 + i], ts[9 + i])
    grads.append(g); C_behind.append(l)
for i in range(SPLIT):
    C_behind.append(bead(tail[i][0], tail[i][1], ts[8 + i], 9 if i == 0 else None))

# ---- layer 2: the ring
C_ring = []
for i in range(8):
    g, l = segment(f'gr{i}', ring[i], ring[(i + 1) % 8], ts[i], ts[(i + 1) % 8])
    grads.append(g); C_ring.append(l)
for i, (x, y) in enumerate(ring, start=1):
    C_ring.append(bead(x, y, ts[i - 1], i if i in (1, 8) else None))

# ---- layer 3: tail in front: links SPLIT..10, beads (9+SPLIT)..20
C_front = []
for i in range(SPLIT, 11):
    g, l = segment(f'gtB{i}', tail[i], tail[i + 1], ts[8 + i], ts[9 + i])
    grads.append(g); C_front.append(l)
for i in range(SPLIT, 12):
    n = 9 + i
    C_front.append(bead(tail[i][0], tail[i][1], ts[8 + i], n if n in (16, 18, 20) else None))

w, h = 900, 372
parts = [f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}">
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
parts += ['<text x="130" y="62" text-anchor="middle" class="lbl">Linear peptide</text>'] + A
parts += ['<text x="130" y="248" text-anchor="middle" class="xs">flexible backbone</text>',
          '<text x="130" y="264" text-anchor="middle" class="xs">proteases cleave anywhere</text>']
parts += ['<text x="340" y="62" text-anchor="middle" class="lbl">Cyclic peptide</text>'] + B
parts += ['<text x="340" y="282" text-anchor="middle" class="xs">closed head-to-tail (the closing bond jumps</text>',
          '<text x="340" y="296" text-anchor="middle" class="xs">from t=0.9 back to t=0 - visible as a colour step)</text>']
parts += [f'<text x="{CX+60}" y="62" text-anchor="middle" class="lbl">Lasso peptide</text>']
parts += C_behind + C_ring + C_front

x1, y1 = ring[0]; x8, y8 = ring[7]
parts += [f'<path d="M {x1:.1f} {y1:.1f} A {R*0.62:.1f} {R*0.62} 0 0 1 {x8:.1f} {y8:.1f}" fill="none" '
          f'stroke="{WARN}" stroke-width="2.4" stroke-dasharray="5 3"/>',
          f'<text x="{CX-96}" y="{CY-46}" text-anchor="end" class="xs" fill="{WARN}">isopeptide bond</text>',
          f'<path d="M {CX-92} {CY-40} L {CX-18} {CY-40}" stroke="{WARN}" stroke-width="1" stroke-dasharray="2 2"/>']
# callouts, now naming the two occlusions separately
parts += [f'<text x="716" y="106" text-anchor="start" class="xs">1. the chain leaves the ring at bead 9</text>',
          f'<path d="M 712 102 L {tail[1][0]+16} {tail[1][1]+4}" fill="none" stroke="{GREY}" '
          f'stroke-width="1.1" marker-end="url(#arw)"/>',
          f'<text x="716" y="138" text-anchor="start" class="xs">2. it turns back and comes down</text>',
          f'<path d="M 712 134 L {tail[4][0]+12} {tail[4][1]}" fill="none" stroke="{GREY}" '
          f'stroke-width="1.1" marker-end="url(#arw)"/>',
          f'<text x="716" y="170" text-anchor="start" class="xs">3. threading the ring - two crossings,</text>',
          f'<text x="716" y="184" text-anchor="start" class="xs">   opposite occlusions:</text>',
          f'<text x="716" y="202" text-anchor="start" class="xs">   (a) here the RING covers the TAIL</text>',
          f'<text x="716" y="220" text-anchor="start" class="xs">   (b) here the TAIL covers the RING</text>',
          f'<text x="716" y="238" text-anchor="start" class="xs">   - together they interlock the fold</text>',
          f'<path d="M 712 196 L {cross[0][0]+12} {cross[0][1]-6}" fill="none" stroke="{GREY}" '
          f'stroke-width="1.1" marker-end="url(#arw)"/>',
          f'<path d="M 700 236 L {cross[1][0]+2} {cross[1][1]+10}" fill="none" stroke="{GREY}" '
          f'stroke-width="1.1" marker-end="url(#arw)"/>']
bx, by, bw_, bh = 250, 316, 400, 13
stops = ''.join(f'<stop offset="{i/20:.2f}" stop-color="{cc(i/20)}"/>' for i in range(21))
parts += [f'<defs><linearGradient id="cbar" x1="{bx}" y1="0" x2="{bx+bw_}" y2="0" '
          f'gradientUnits="userSpaceOnUse">{stops}</linearGradient></defs>',
          f'<rect x="{bx}" y="{by}" width="{bw_}" height="{bh}" rx="3" fill="url(#cbar)"/>',
          f'<text x="{bx-12}" y="{by+11}" text-anchor="end" class="xs" font-weight="bold">N-terminus</text>',
          f'<text x="{bx+bw_+12}" y="{by+11}" text-anchor="start" class="xs" font-weight="bold">C-terminus</text>',
          f'<text x="{w/2}" y="{by+34}" text-anchor="middle" class="xs">'
          f'Colour encodes position along the chain; each bond is drawn as its own gradient.</text>',
          '</svg>']
open(os.path.join(FIG, 'fig01_lasso_topology.svg'), 'w', encoding='utf-8').write('\n'.join(parts))
print('fig01 v11 written')

import cairosvg
cairosvg.svg2png(url=os.path.join(FIG, 'fig01_lasso_topology.svg'),
                 write_to=os.path.join(FIG, 'fig01_lasso_topology_preview.png'),
                 scale=2.0, background_color='white')
print('rendered')
