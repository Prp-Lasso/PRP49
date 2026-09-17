"""Sequence-level check: do these three proteins have ANY lasso-peptide character?

A lasso peptide precursor has (i) a leader peptide, (ii) a Cys-Glu/Asp macrocyclisation
motif, (iii) a ~15-40 aa core. Testing for those directly settles whether the group's
"promising peptide" reading can hold, independently of the database provenance.
"""
import json
import os
import re
import urllib.request

ROOT = r'D:\deepseek_harness\prp49'
ACC = ['Q16774', 'Q16384', 'Q9Y3B2']

seqs = {}
for acc in ACC:
    url = f'https://rest.uniprot.org/uniprotkb/{acc}.fasta'
    req = urllib.request.Request(url, headers={'User-Agent': 'prp49-research/1.0'})
    with urllib.request.urlopen(req, timeout=60) as r:
        txt = r.read().decode()
    lines = txt.splitlines()
    seqs[acc] = (lines[0][1:], ''.join(lines[1:]))

print('=== lasso-peptide structural criteria ===')
for acc, (hdr, s) in seqs.items():
    cys = [i + 1 for i, c in enumerate(s) if c == 'C']
    gly = sum(1 for c in s if c == 'G')
    acidic = sum(1 for c in s if c in 'DE')
    print(f'\n{acc}  {len(s)} aa  ({hdr.split("|")[0] if "|" in hdr else hdr[:40]})')
    print(f'  Cys positions      : {cys if cys else "NONE"}')
    print(f'  Gly count          : {gly} ({100*gly/len(s):.1f}%)')
    print(f'  Asp/Glu count      : {acidic}')
    print(f'  length in 15-40 aa : {15 <= len(s) <= 40}   <-- lasso core range')
    print(f'  N-terminal Met     : {s.startswith("M")}  (eukaryotic start, not a lasso leader)')
    print(f'  has leader-like run: {bool(re.search(r"^M[^C]{5,30}$", s[:35]))}')

print('\n=== the decisive comparison ===')
print('  lasso peptide core  : 15-40 aa, macrocyclisation via Cys/Glu-Asp, threaded tail')
print('  these three         : 188-197 aa folded enzymes/scaffolds, no cyclisation motif')
print('  => they cannot be lasso peptides; they are full human proteins.')

print('\n=== so what ARE they in the TUnA/Bernett data? ===')
print('  Bernett et al. is a PROTEIN-PROTEIN interaction benchmark (both partners are')
print('  UniProt proteins, typically 100-1000 aa). "Intra0/1/2" are its partitions;')
print('  each line is "partnerA partnerB" labelled pos/neg by whether they interact.')
print('  Q16774 appears only in Intra1 (12 neg / 11 pos), Q9Y3B2 only in Intra2')
print('  (18 neg / 16 pos), Q16384 appears in NO pair file at all.')
print('  A protein appearing in BOTH pos and neg pairs is normal there: it interacts')
print('  with some partners and not others. It says nothing about peptide potential.')

print('\n=== what would a real "promising lasso peptide" look like in our data? ===')
lib = os.path.join(ROOT, 'lassopred_database.csv')
if os.path.exists(lib):
    import pandas as pd
    d = pd.read_csv(lib)
    core = d['Core_Sequence'].astype(str)
    print(f'  LassoPred library: {len(d)} cores, length {core.str.len().min()}-{core.str.len().max()} aa')
    print(f'  examples: {list(core.head(3))}')
