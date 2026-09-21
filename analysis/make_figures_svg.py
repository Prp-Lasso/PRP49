"""Generate the 5 schematic figures as hand-authored SVG (crisp at any scale, editable).

English labels so they drop straight into the paper. Colours match the data figures.
"""
import os

ROOT = r'D:\deepseek_harness\prp49'
FIG = os.path.join(ROOT, 'docs', 'figures')
os.makedirs(FIG, exist_ok=True)

MAIN, ACC, WARN, GREY, INK = '#2b6cb0', '#dd6b20', '#c53030', '#718096', '#1a202c'
FONT = 'Helvetica, Arial, sans-serif'


def head(w, h, title):
    return f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}">
<rect width="{w}" height="{h}" fill="white"/>
<style>
  text {{ font-family: {FONT}; fill: {INK}; }}
  .t  {{ font-size: 15px; font-weight: bold; }}
  .s  {{ font-size: 11.5px; }}
  .xs {{ font-size: 10px; fill: #4a5568; }}
  .lbl {{ font-size: 12px; font-weight: bold; }}
</style>
<text x="{w/2}" y="26" text-anchor="middle" class="t">{title}</text>
'''


# ============ fig 01: lasso topology vs linear vs simple cyclic ============
w, h = 900, 340
s = head(w, h, 'Lasso peptide topology: an isopeptide ring with a threaded C-terminal tail')


def pep_chain(x0, y0, n, r=9, gap=26, color=MAIN):
    """draw a peptide as beads on a horizontal line; returns bead centres"""
    pts = [(x0 + i * gap, y0) for i in range(n)]
    out = []
    for i, (x, y) in enumerate(pts):
        out.append(f'<circle cx="{x}" cy="{y}" r="{r}" fill="{color}" opacity="0.9"/>')
        if i:
            px = pts[i - 1][0]
            out.append(f'<line x1="{px+r}" y1="{y0}" x2="{x-r}" y2="{y0}" stroke="{color}" stroke-width="2.4"/>')
    return '\n'.join(out), pts


# --- panel A: linear peptide
s += f'<text x="150" y="70" text-anchor="middle" class="lbl">Linear peptide</text>'
chain, pts = pep_chain(45, 115, 9)
s += chain
s += f'<text x="150" y="152" text-anchor="middle" class="xs">flexible backbone</text>'
s += f'<text x="150" y="168" text-anchor="middle" class="xs">proteases cleave anywhere</text>'

# --- panel B: simple cyclic peptide
s += f'<text x="450" y="70" text-anchor="middle" class="lbl">Cyclic peptide</text>'
cx, cy, R = 450, 125, 48
s += f'<circle cx="{cx}" cy="{cy}" r="{R}" fill="none" stroke="{MAIN}" stroke-width="9" opacity="0.9"/>'
s += f'<text x="{cx}" y="{cy+4}" text-anchor="middle" class="xs">head-to-tail</text>'
s += f'<text x="450" y="196" text-anchor="middle" class="xs">rigid, but no threaded topology</text>'

# --- panel C: lasso peptide (the real thing)
s += f'<text x="740" y="70" text-anchor="middle" class="lbl" fill="{ACC}">Lasso peptide</text>'
cx, cy, R = 740, 130, 52
# ring, drawn in two arcs so the tail can pass "through" it visually
s += f'<path d="M {cx-R} {cy} A {R} {R} 0 0 1 {cx+R} {cy}" fill="none" stroke="{ACC}" stroke-width="9" opacity="0.85"/>'
# tail: comes up from the ring, over the front, and out to the right (threading)
s += (f'<path d="M {cx} {cy+6} C {cx-24} {cy-30}, {cx+26} {cy-58}, {cx+4} {cy-78} '
      f'L {cx+92} {cy-96}" fill="none" stroke="{ACC}" stroke-width="7" stroke-linecap="round"/>')
# front arc of ring drawn after the tail => tail appears to pass behind then in front
s += f'<path d="M {cx-R} {cy} A {R} {R} 0 0 0 {cx+R} {cy}" fill="none" stroke="{ACC}" stroke-width="9" opacity="0.85"/>'
# isopeptide bond marker
s += (f'<circle cx="{cx-R*0.72}" cy="{cy-36}" r="6.5" fill="{WARN}"/>'
      f'<text x="{cx-R*0.72-14}" y="{cy-52}" text-anchor="middle" class="xs" fill="{WARN}">isopeptide</text>'
      f'<text x="{cx-R*0.72-14}" y="{cy-40}" text-anchor="middle" class="xs" fill="{WARN}">bond</text>')
s += f'<text x="740" y="206" text-anchor="middle" class="xs">tail threaded through the ring</text>'
s += f'<text x="740" y="222" text-anchor="middle" class="xs">&#8594; proteases blocked, thermally locked</text>'

s += f'<text x="{w/2}" y="{h-14}" text-anchor="middle" class="xs">The threaded [1]rotaxane fold is why generic structure predictors fail on this class.</text>'
s += '</svg>'
open(os.path.join(FIG, 'fig01_lasso_topology.svg'), 'w', encoding='utf-8').write(s)
print('fig01 done')


# ============ fig 02: biosynthesis ============
w, h = 940, 300
s = head(w, h, 'Lasso peptide biosynthesis: cyclisation then threading')
stages = [
    ('Precursor peptide', 'leader + core', GREY),
    ('Macrolactam ring', 'lasso cyclase\ncloses N-term to Asp/Glu', MAIN),
    ('Threading', 'C-terminal tail\npasses through the ring', ACC),
    ('Mature lasso peptide', 'topology locked', WARN),
]
bw, gap = 190, 30
x = 30
for i, (t, sub, col) in enumerate(stages):
    s += (f'<rect x="{x}" y="80" width="{bw}" height="110" rx="12" fill="{col}" opacity="0.10" '
          f'stroke="{col}" stroke-width="2"/>')
    s += f'<text x="{x+bw/2}" y="118" text-anchor="middle" class="lbl" fill="{col}">{t}</text>'
    for j, line in enumerate(sub.split('\n')):
        s += f'<text x="{x+bw/2}" y="{142+j*16}" text-anchor="middle" class="xs">{line}</text>'
    if i < len(stages) - 1:
        s += (f'<path d="M {x+bw+6} 135 L {x+bw+gap-8} 135" stroke="{GREY}" stroke-width="2.4" '
              f'marker-end="url(#ah)"/>')
    x += bw + gap
s += f'''<defs><marker id="ah" markerWidth="9" markerHeight="9" refX="7" refY="4.5" orient="auto">
<path d="M0,0 L9,4.5 L0,9 z" fill="{GREY}"/></marker></defs>'''
s += f'<text x="{w/2}" y="232" text-anchor="middle" class="xs">Only ~50 distinct lasso peptides have been structurally characterised in 30 years,</text>'
s += f'<text x="{w/2}" y="250" text-anchor="middle" class="xs">while genome mining proposes thousands of sequences.</text>'
s += '</svg>'
open(os.path.join(FIG, 'fig02_biosynthesis.svg'), 'w', encoding='utf-8').write(s)
print('fig02 done')


# ============ fig 04: generic predictors fail, LassoPred works ============
w, h = 880, 330
s = head(w, h, 'Why generic structure predictors fail, and what LassoPred adds')
s += f'<line x1="{w/2}" y1="52" x2="{w/2}" y2="{h-30}" stroke="#e2e8f0" stroke-width="1.5"/>'
s += f'<text x="{w/4}" y="76" text-anchor="middle" class="lbl" fill="{WARN}">Generic tools</text>'
tools = ['AlphaFold2', 'AlphaFold3', 'ESMFold']
for i, t in enumerate(tools):
    y = 112 + i * 46
    s += (f'<rect x="60" y="{y}" width="270" height="36" rx="8" fill="{WARN}" opacity="0.08" '
          f'stroke="{WARN}" stroke-width="1.6"/>')
    s += f'<text x="82" y="{y+24}" class="s">{t}</text>'
    s += f'<text x="316" y="{y+24}" text-anchor="end" class="s" fill="{WARN}">&#10007; knot + isopeptide</text>'
s += f'<text x="{w/4}" y="268" text-anchor="middle" class="xs">irregular lariat-knot fold</text>'
s += f'<text x="{w/4}" y="286" text-anchor="middle" class="xs">non-standard isopeptide bond</text>'

s += f'<text x="{3*w/4}" y="76" text-anchor="middle" class="lbl" fill="{MAIN}">LassoPred (Zhao / Vanderbilt / Princeton)</text>'
s += (f'<rect x="500" y="100" width="330" height="76" rx="10" fill="{MAIN}" opacity="0.08" '
      f'stroke="{MAIN}" stroke-width="1.8"/>')
s += f'<text x="520" y="126" class="s">1. classifier &#8594; ring / loop / tail</text>'
s += f'<text x="520" y="150" class="s">2. constructor &#8594; assemble 3D model</text>'
s += f'<text x="520" y="172" class="xs">(homology modelling + mutation + minimisation)</text>'
s += (f'<rect x="500" y="192" width="330" height="52" rx="10" fill="{ACC}" opacity="0.12" '
      f'stroke="{ACC}" stroke-width="1.8"/>')
s += f'<text x="520" y="216" class="lbl" fill="{ACC}">4,749 predicted structures</text>'
s += f'<text x="520" y="234" class="xs">vs &lt;50 experimentally solved in 30 years</text>'
s += f'<text x="{3*w/4}" y="268" text-anchor="middle" class="xs">near-experimental accuracy, minutes per sequence</text>'
s += '</svg>'
open(os.path.join(FIG, 'fig04_structure_prediction.svg'), 'w', encoding='utf-8').write(s)
print('fig04 done')


# ============ fig 05: project pipeline ============
w, h = 900, 520
s = head(w, h, 'Project pipeline: from structure prediction to drug-potential mining')
steps = [
    ('LassoPred output', '4,749 predicted lasso peptide structures + topology annotations', GREY),
    ('Preliminary docking', 'AutoDock Vina: lasso peptide x target feasibility, box handling for macrocycles', MAIN),
    ('Lasso vs linear classifier', 'sequence-level filter: is this really a lasso topology?', MAIN),
    ('Domain-specific embedder', 'LassoESM (650M) for peptides + ESM-2 for targets', ACC),
    ('Bilinear attention model', 'BAN fusion &#8594; binding classifier + affinity regressor', ACC),
    ('Scarcity handling', 'family-homolog augmentation of positives; family-grouped CV', WARN),
    ('Constrained regression', 'output-constraint variants to stop divergence', WARN),
    ('Screening + release', '13 peptides x 50 targets &#8594; 10 candidates; open-sourced', '#2f855a'),
]
y = 62
for i, (t, sub, col) in enumerate(steps):
    s += (f'<rect x="70" y="{y}" width="760" height="46" rx="10" fill="{col}" opacity="0.09" '
          f'stroke="{col}" stroke-width="1.8"/>')
    s += f'<circle cx="100" cy="{y+23}" r="15" fill="{col}" opacity="0.85"/>'
    s += f'<text x="100" y="{y+28}" text-anchor="middle" class="s" fill="white" font-weight="bold">{i+1}</text>'
    s += f'<text x="128" y="{y+21}" class="lbl" fill="{col}">{t}</text>'
    s += f'<text x="128" y="{y+37}" class="xs">{sub}</text>'
    if i < len(steps) - 1:
        s += f'<path d="M 450 {y+47} L 450 {y+58}" stroke="{GREY}" stroke-width="2.2" marker-end="url(#ah2)"/>'
    y += 56
s += f'''<defs><marker id="ah2" markerWidth="9" markerHeight="9" refX="7" refY="4.5" orient="auto">
<path d="M0,0 L9,4.5 L0,9 z" fill="{GREY}"/></marker></defs>'''
s += '</svg>'
open(os.path.join(FIG, 'fig05_pipeline.svg'), 'w', encoding='utf-8').write(s)
print('fig05 done')


# ============ fig 06: model architecture ============
w, h = 900, 430
s = head(w, h, 'Dual-encoder bilinear-attention architecture')
# inputs
s += f'<rect x="40" y="70" width="200" height="52" rx="10" fill="{MAIN}" opacity="0.10" stroke="{MAIN}" stroke-width="1.8"/>'
s += f'<text x="140" y="92" text-anchor="middle" class="lbl" fill="{MAIN}">Target sequence</text>'
s += f'<text x="140" y="110" text-anchor="middle" class="xs">up to 1024 aa</text>'
s += f'<rect x="40" y="230" width="200" height="52" rx="10" fill="{ACC}" opacity="0.10" stroke="{ACC}" stroke-width="1.8"/>'
s += f'<text x="140" y="252" text-anchor="middle" class="lbl" fill="{ACC}">Lasso peptide</text>'
s += f'<text x="140" y="270" text-anchor="middle" class="xs">15-40 aa, threaded topology</text>'
# encoders
s += f'<rect x="300" y="70" width="190" height="52" rx="10" fill="{MAIN}" opacity="0.85"/>'
s += f'<text x="395" y="94" text-anchor="middle" class="s" fill="white" font-weight="bold">ESM-2 (35M)</text>'
s += f'<text x="395" y="110" text-anchor="middle" class="xs" fill="white">protein encoder</text>'
s += f'<rect x="300" y="230" width="190" height="52" rx="10" fill="{ACC}" opacity="0.85"/>'
s += f'<text x="395" y="254" text-anchor="middle" class="s" fill="white" font-weight="bold">LassoESM (650M)</text>'
s += f'<text x="395" y="270" text-anchor="middle" class="xs" fill="white">lasso-specific embedder</text>'
# BAN
s += f'<rect x="560" y="140" width="200" height="76" rx="12" fill="#2f855a" opacity="0.12" stroke="#2f855a" stroke-width="2"/>'
s += f'<text x="660" y="168" text-anchor="middle" class="lbl" fill="#2f855a">Bilinear attention</text>'
s += f'<text x="660" y="186" text-anchor="middle" class="xs">h=256, heads=3, k=2</text>'
s += f'<text x="660" y="204" text-anchor="middle" class="xs">(BAN, low-rank factorised)</text>'
# heads
s += f'<rect x="560" y="300" width="200" height="46" rx="10" fill="{MAIN}" opacity="0.12" stroke="{MAIN}" stroke-width="1.8"/>'
s += f'<text x="660" y="320" text-anchor="middle" class="s" fill="{MAIN}">Binding classifier</text>'
s += f'<text x="660" y="337" text-anchor="middle" class="xs">BCE + alignment + ranking</text>'
s += f'<rect x="560" y="360" width="200" height="46" rx="10" fill="{WARN}" opacity="0.12" stroke="{WARN}" stroke-width="1.8"/>'
s += f'<text x="660" y="380" text-anchor="middle" class="s" fill="{WARN}">Affinity regressor</text>'
s += f'<text x="660" y="397" text-anchor="middle" class="xs">constrained output</text>'
# connectors
for y0, y1 in [(96, 140), (256, 216)]:
    s += f'<path d="M 240 {y0} L 300 {y0}" stroke="{GREY}" stroke-width="2.2" marker-end="url(#ah3)"/>'
s += f'<path d="M 490 96 C 530 96, 530 160, 560 168" fill="none" stroke="{GREY}" stroke-width="2.2" marker-end="url(#ah3)"/>'
s += f'<path d="M 490 256 C 530 256, 530 200, 560 192" fill="none" stroke="{GREY}" stroke-width="2.2" marker-end="url(#ah3)"/>'
s += f'<path d="M 660 216 L 660 300" stroke="{GREY}" stroke-width="2.2" marker-end="url(#ah3)"/>'
s += f'<path d="M 700 300 C 730 300, 730 380, 760 383" fill="none" stroke="{GREY}" stroke-width="2.2" marker-end="url(#ah3)"/>'
s += f'''<defs><marker id="ah3" markerWidth="9" markerHeight="9" refX="7" refY="4.5" orient="auto">
<path d="M0,0 L9,4.5 L0,9 z" fill="{GREY}"/></marker></defs>'''
s += '</svg>'
open(os.path.join(FIG, 'fig06_architecture.svg'), 'w', encoding='utf-8').write(s)
print('fig06 done')

print(f'\nall schematics in {FIG}')
