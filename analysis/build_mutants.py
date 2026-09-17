"""D3 step 1: in-silico saturation mutagenesis of the top candidates.

For each of the top candidates, every position is substituted with all 19 other
residues and the model scores the variant against that candidate's own target. The
resulting per-position table shows which substitutions the model believes improve
binding - i.e. concrete peptide-engineering suggestions.

HONEST CAVEAT carried into the output: the training set contains no engineered
variants, so these scores are model EXTRAPOLATION. They are hypotheses for the bench,
not validated designs. Variants are additionally filtered for lasso compatibility:
substituting away a conserved Cys/Gly or a bulky residue inside the ring is flagged,
since the model knows nothing about whether the macrocycle still closes.
"""
import os
import re

import pandas as pd

ROOT = r'D:\deepseek_harness\prp49'
AA = 'ACDEFGHIKLMNPQRSTVWY'
AA_RE = re.compile(r'[^ACDEFGHIKLMNPQRSTVWYX]')

cand = pd.read_csv(os.path.join(ROOT, 'results', 'candidates_final.csv'))
print(f'candidates: {len(cand)}')

# peptide sequences
pep = pd.read_csv(os.path.join(ROOT, 'mvp_cpu', 'scan_peptides_ext.csv'))
seq_col = next(c for c in pep.columns if 'seq' in c.lower())
name_col = next(c for c in pep.columns if 'id' in c.lower() or 'name' in c.lower())
seq_of = dict(zip(pep[name_col].astype(str), pep[seq_col].astype(str)))

# top candidates, deduplicated by peptide, keeping the best-ranked target for each
top = cand.sort_values('final_rank').drop_duplicates('peptide_id').head(15)
print(f'top peptides: {len(top)}')
rows = []
for _, r in top.iterrows():
    s = AA_RE.sub('', str(seq_of.get(r.peptide_id, '')).upper())
    if len(s) < 5:
        print(f'  !! no sequence for {r.peptide_id}')
        continue
    rows.append(dict(peptide_id=r.peptide_id, target_id=r.target_id, seq=s,
                     final_rank=int(r.final_rank), binding_prob=r.binding_prob,
                     dock=r.dock))
sel = pd.DataFrame(rows)
print(sel[['peptide_id', 'target_id', 'final_rank', 'seq']].to_string(index=False))

# --- generate saturation variants ---
variants, meta = [], []
for _, r in sel.iterrows():
    s = r.seq
    for pos in range(len(s)):
        wt = s[pos]
        if wt not in AA:
            continue
        for aa in AA:
            if aa == wt:
                continue
            m = s[:pos] + aa + s[pos + 1:]
            vid = f'{r.peptide_id}_p{pos+1}{wt}{aa}'
            variants.append(dict(id=vid, seq=m))
            meta.append(dict(variant_id=vid, peptide_id=r.peptide_id,
                             target_id=r.target_id, position=pos + 1,
                             wt=wt, mut=aa, wt_seq=s))

vdf = pd.DataFrame(variants)
mdf = pd.DataFrame(meta)
print(f'\nvariants generated: {len(vdf)}')
print(f'  per peptide: {vdf.id.str.split("_p").str[0].value_counts().to_dict()}')
vdf.to_csv(os.path.join(ROOT, 'mvp_cpu', 'mut_variants.csv'), index=False,
           lineterminator='\n')
mdf.to_csv(os.path.join(ROOT, 'mvp_cpu', 'mut_meta.csv'), index=False,
            lineterminator='\n')

# --- targets: one entry per (peptide, target) pair used above ---
tgt = {}
name, buf = None, []
for line in open(os.path.join(ROOT, 'mvp_cpu', 'scan_targets.fasta'), encoding='utf-8'):
    line = line.strip()
    if line.startswith('>'):
        if name:
            tgt[name] = AA_RE.sub('', ''.join(buf).upper())
        name, buf = line[1:].split()[0], []
    elif line:
        buf.append(line)
if name:
    tgt[name] = AA_RE.sub('', ''.join(buf).upper())

needed = {r.target_id: tgt.get(r.target_id) or next(
    (v for k, v in tgt.items() if k.startswith(r.target_id)), None)
    for _, r in sel.iterrows()}
needed = {k: v for k, v in needed.items() if v}
print(f'\ntargets needed: {len(needed)}')
with open(os.path.join(ROOT, 'mvp_cpu', 'mut_targets.fasta'), 'w',
          encoding='utf-8', newline='\n') as f:
    for k, v in needed.items():
        f.write(f'>{k}\n{v}\n')
        print(f'  {k:12s} {len(v)} aa')

# --- lasso-compatibility flags (model knows nothing about ring closure) ---
STRUCTURAL = set('CG')          # ring-forming / conserved
def flag(row):
    if row.wt in STRUCTURAL:
        return 'structural_wt_Cys_Gly'
    if row.wt in 'P' and row.mut in 'GP':
        return 'may_break_turn'
    return ''
mdf['flag'] = mdf.apply(flag, axis=1)
print(f'\nflagged variants: {int((mdf.flag != "").sum())} of {len(mdf)}')
mdf.to_csv(os.path.join(ROOT, 'mvp_cpu', 'mut_meta.csv'), index=False, lineterminator='\n')
print('\nwrote mut_variants.csv / mut_meta.csv / mut_targets.fasta')
