"""Route D database, step 2: extract (lasso peptide, target protein) pairs.

Most of the 98 PDB entries are lasso peptides alone (NMR) or bound to their maturation
machinery / export pumps (isopeptidases, McjD). Route D needs the biologically relevant
complexes, so entries are filtered to those holding both a short peptide chain and a
larger protein chain, then the peptide is matched against the LassoPred library so the
pair is anchored to a known lasso sequence.

Output: lasso_target_db/lasso_target_pairs.csv
"""
import ast
import json
import os
import re

import pandas as pd

ROOT = r'D:\deepseek_harness\prp49'
DB = os.path.join(ROOT, 'lasso_target_db')
AA = re.compile(r'[^ACDEFGHIKLMNPQRSTVWY]')

ent = pd.read_csv(os.path.join(DB, 'pdb_lasso_entries.csv'))
print(f'entries: {len(ent)}')

# --- LassoPred library to anchor peptide identity ---
lib = pd.read_csv(os.path.join(ROOT, 'lassopred_database.csv'))
lib['core'] = lib.Core_Sequence.astype(str).str.upper().map(lambda s: AA.sub('', s))
lib = lib[lib.core.str.len() >= 10]
lib_seqs = set(lib.core)
lib_names = {}
for _, r in lib.iterrows():
    nm = str(r.get('Lasso_Peptide_Name', '') or '')
    if nm and nm.lower() not in ('nan', 'none'):
        lib_names.setdefault(r.core, nm)
print(f'LassoPred cores: {len(lib_seqs)}')

MACHINERY = ['isopeptidase', 'peptidase', 'transporter', 'atp-binding', 'abc transporter',
             'maturase', 'kinase', 'synthase', 'protease', 'isomerase', 'abc-type',
             'permease', 'atpase']
TARGET_HINT = ['polymerase', 'rnap', 'rna polymerase', 'clpc', 'clpp', 'clpb', 'ribosome',
               's50', 's30', 'lipid ii', 'pbp', 'dna gyrase', 'topoisomerase',
               'acetyltransferase', 'phosphatase', 'receptor', 'nmr structure of']


def peptide_matches_lib(seq):
    """Exact or near-exact containment of a known lasso core."""
    if seq in lib_seqs:
        return seq, 1.0
    for ls in lib_seqs:
        if len(ls) < 12:
            continue
        if ls in seq or seq in ls:
            shorter, longer = (ls, seq) if len(ls) < len(seq) else (seq, ls)
            return ls, len(shorter) / len(longer)
    return None, 0.0


rows = []
for _, e in ent.iterrows():
    try:
        chains = json.loads(e.chains)
    except Exception:
        continue
    if not isinstance(chains, list) or len(chains) < 2:
        continue
    short = [c for c in chains if c.get('len', 0) <= 80]
    long_ = [c for c in chains if c.get('len', 0) > 120]
    if not short or not long_:
        continue
    title = str(e.title).lower()
    is_machinery = any(m in title for m in MACHINERY)
    for s in short:
        core, cov = peptide_matches_lib(s['seq'])
        if core is None:
            continue
        for t in long_:
            desc = str(t.get('desc') or '').lower()
            hint = any(h in desc or h in title for h in TARGET_HINT)
            rows.append(dict(
                pdb_id=e.pdb_id, title=e.title, method=e.method,
                lasso_seq=core, lasso_len=len(core), match_cov=round(cov, 3),
                lasso_name=lib_names.get(core, ''),
                target_seq=t['seq'], target_len=t['len'],
                target_desc=t.get('desc'), target_organism=t.get('organism'),
                likely_machinery=is_machinery, target_keyword_hit=hint,
            ))

df = pd.DataFrame(rows).drop_duplicates(subset=['pdb_id', 'lasso_seq', 'target_seq'])
out = os.path.join(DB, 'lasso_target_pairs.csv')
df.to_csv(out, index=False, lineterminator='\n')
print(f'\nwrote {out}: {len(df)} candidate (lasso, target) pairs from {df.pdb_id.nunique() if len(df) else 0} entries')

if len(df):
    print(f'  flagged as maturation/transport machinery: {int(df.likely_machinery.sum())}')
    print(f'  with a target keyword in the description: {int(df.target_keyword_hit.sum())}')
    print('\nby peptide:')
    print(df.groupby('lasso_name').size().sort_values(ascending=False).head(12).to_string())
    bio = df[~df.likely_machinery]
    print(f'\n=== biologically plausible targets (not machinery): {len(bio)} ===')
    for _, r in bio.head(15).iterrows():
        print(f"  {r.pdb_id} {str(r.lasso_name)[:14]:14s} -> {str(r.target_desc)[:60]} "
              f"({r.target_len} aa, {r.target_organism})")
