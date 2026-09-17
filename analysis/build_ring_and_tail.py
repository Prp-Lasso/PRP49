"""D3 steps 2+3: ring-size variants and C-terminal modification advice.

Ring/tail assignment comes from the literature topology of each characterised lasso
peptide (the macrolactam is formed between the alpha-amino group of the N-terminal
residue and the side-chain carboxyl of an Asp/Glu several residues in; the remainder
threads through). Where a peptide is not characterised, the default of an 8-residue
ring is used and that assumption is recorded per row.

Ring-size variants change only the ring segment, then the model scores them - the
stability side of the trade-off comes from the known topology rules (a 7-9 residue
ring is what permits threading), and is reported separately rather than invented by
the model, which has no notion of whether a macrocycle still closes.

C-terminal modifications: unlike an ordinary peptide, a lasso peptide's N-terminus is
already locked in the macrolactam, so N-terminal capping is not available; the exposed
C-terminal tail is the main protease handle, making amidation the highest-value edit.
"""
import os
import re

import pandas as pd

ROOT = r'D:\deepseek_harness\prp49'
AA_RE = re.compile(r'[^ACDEFGHIKLMNPQRSTVWYX]')

# literature topology: peptide -> (ring_end_index_1based, source)
TOPOLOGY = {
    'MccJ25':      (8,  'literature: lactam Gly1-E8, tail Y9-G21 threads; His5 H-bonds tail'),
    'Capistruin':  (9,  'literature: lactam Gly1-D9, tail R10-N19'),
    'Lassomycin':  (8,  'literature: lactam Gly1-D8, tail Q9-X16'),
    'RES-701-1':   (9,  'literature: 9-residue ring (RES series)'),
    'RES-701-3':   (9,  'literature: 9-residue ring (RES series)'),
    'Anantin':     (8,  'PREDICTED structure (LassoPred+tleap, 8 loop-length models)'),
    'Siamycin-I':  (8,  'assumed default (uncharacterised topology)'),
    'Sphingopyxin-I': (9, 'assumed default'),
    'Ubonodin':    (8,  'assumed default'),
    'Chaxapeptin': (8,  'assumed default'),
    'Lariocidin':  (8,  'assumed default'),
}
DEFAULT_RING = (8, 'assumed default')

cand = pd.read_csv(os.path.join(ROOT, 'results', 'candidates_final.csv'))
pep = pd.read_csv(os.path.join(ROOT, 'mvp_cpu', 'scan_peptides_ext.csv'))
seq_col = next(c for c in pep.columns if 'seq' in c.lower())
name_col = next(c for c in pep.columns if 'id' in c.lower() or 'name' in c.lower())
seq_of = {str(k): AA_RE.sub('', str(v).upper()) for k, v in zip(pep[name_col], pep[seq_col])}

top = cand.sort_values('final_rank').drop_duplicates('peptide_id').head(12)

# ---------------- (2) ring-size variants ----------------
variants, meta = [], []
for _, r in top.iterrows():
    s = seq_of.get(r.peptide_id, '')
    if len(s) < 10:
        continue
    ring_end, src = TOPOLOGY.get(r.peptide_id, DEFAULT_RING)
    ring_end = min(ring_end, len(s) - 3)
    ring = s[:ring_end]
    tail = s[ring_end:]
    for delta in (-1, +1, +2):
        if delta < 0:
            if len(ring) <= 5:
                continue
            new_ring = ring[:delta]                       # delete last ring residue
            kind = f'ring{len(ring)}to{len(new_ring)}'
        else:
            new_ring = ring + 'G' * delta                 # insert Gly to lengthen
            kind = f'ring{len(ring)}to{len(new_ring)}'
        m = new_ring + tail
        vid = f'{r.peptide_id}_{kind}'
        variants.append(dict(id=vid, seq=m))
        meta.append(dict(variant_id=vid, peptide_id=r.peptide_id, target_id=r.target_id,
                         wt_seq=s, variant_seq=m, ring_wt=len(ring), ring_new=len(new_ring),
                         delta=delta, topology_source=src,
                         stability_note=('7-9 residue ring: threading-competent'
                                         if 7 <= len(new_ring) <= 9 else
                                         'outside the 7-9 window: threading likely lost')))
rv = pd.DataFrame(variants)
rm = pd.DataFrame(meta)
print(f'ring-size variants: {len(rv)}')
print(rm[['variant_id', 'ring_wt', 'ring_new', 'stability_note']].to_string(index=False))
rv.to_csv(os.path.join(ROOT, 'mvp_cpu', 'ring_variants.csv'), index=False, lineterminator='\n')
rm.to_csv(os.path.join(ROOT, 'mvp_cpu', 'ring_meta.csv'), index=False, lineterminator='\n')

# ---------------- (3) C-terminal modification advice ----------------
GEN = [
    ('C_term_amidation', 'C-terminal carboxyl -> amide (CONH2)',
     'removes the carboxypeptidase substrate; the exposed tail is the main protease handle',
     'usually neutral: the tail points away from the binding interface in characterised complexes',
     'trivial (standard SPPS)', 'high'),
    ('C_term_methylation', 'C-terminal carboxyl -> methyl ester/amide',
     'same protection as amidation with a different synthesis route',
     'neutral to slightly negative (ester is hydrolytically labile)',
     'trivial', 'low'),
    ('N_term_capping', 'NOT APPLICABLE',
     "the N-terminal alpha-amino group is already consumed by the macrolactam",
     'n/a - this is an intrinsic advantage of the lasso fold',
     'n/a', 'n/a'),
    ('tail_truncation', 'delete 1-2 exposed tail residues',
     'shortens the protease-accessible segment',
     'RISK: the tail threads the ring; shortening can unthread the peptide and destroy the fold',
     'moderate', 'medium-low'),
    ('D-amino_acid_scan', 'replace tail residues with D-configured analogues',
     'proteases cannot cleave D-peptide bonds',
     'depends on the residue; tail is peripheral in known complexes',
     'moderate (chiral synthesis)', 'medium'),
    ('head-to-tail_cyclisation', 'covalent closure of the tail to the ring',
     'removes the free tail entirely - maximum protease resistance',
     'HIGH RISK: converts a lasso into a simple cycle, discarding the threaded topology',
     'hard', 'low'),
]
rows = []
for _, r in top.iterrows():
    s = seq_of.get(r.peptide_id, '')
    if len(s) < 8:
        continue
    ring_end, src = TOPOLOGY.get(r.peptide_id, DEFAULT_RING)
    ring_end = min(ring_end, len(s) - 3)
    cterm = s[-1]
    tail = s[ring_end:]
    for name, what, why, effect, cost, priority in GEN:
        rows.append(dict(peptide_id=r.peptide_id, final_rank=int(r.final_rank),
                         target_id=r.target_id, c_term_residue=cterm,
                         tail=f'{tail}({len(tail)} aa)', modification=name,
                         what_it_does=what, rationale=why,
                         expected_effect_on_binding=effect, synthesis_cost=cost,
                         priority=priority, topology_source=src))
adv = pd.DataFrame(rows)
adv.to_csv(os.path.join(ROOT, 'results', 'tail_modification_advice.csv'),
           index=False, lineterminator='\n')
print(f'\nwrote tail_modification_advice.csv: {len(adv)} rows')
print(adv.groupby('modification').agg(n=('peptide_id', 'size'),
                                      priority=('priority', 'first')).to_string())
print('\nC-terminal residues per candidate:')
print(adv.drop_duplicates('peptide_id')[['peptide_id', 'c_term_residue', 'tail', 'final_rank']]
      .to_string(index=False))
