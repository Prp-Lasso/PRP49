"""Rebuild the route D datasets on the nucleic-acid-free entity table, then build
the verification scan inputs.

Clean peptide chains recovered:
  6N60  MccJ25 (21 aa)      x RNAP beta-prime -> human POLR2A
  6N61  capistruin (19 aa)  x RNAP beta-prime -> human POLR2A
  8IBO  lassomycin (16 aa)  x ClpC1            -> human CLPB
  8IBP  lassomycin (16 aa)  x ClpC1            -> human CLPB
  4CU4  MccJ25 (21 aa)      x FhuA             -> uptake receptor
  9KDF  RES-701-3 (16 aa)   x EDNRB            -> human EDNRB (grade A reference)
Both MccJ25 and capistruin are scan peptides, so this doubles as a held-out check.
"""
import os
import re

import pandas as pd

ROOT = r'D:\deepseek_harness\prp49'
DB = os.path.join(ROOT, 'lasso_target_db')
AA = re.compile(r'[^ACDEFGHIKLMNPQRSTVWYX]')

ent = pd.read_csv(os.path.join(DB, 'rnap_lasso_entities.csv'))
ent = ent[ent.entity_type.str.startswith('polypeptide')].copy()
ent['seq'] = ent.seq.astype(str).str.upper().map(lambda s: AA.sub('', s))

PRIMARY = {'6N60': "subunit beta'", '6N61': "subunit beta'", '6N62': "subunit beta'",
           '8IBO': 'ClpC/mecB', '8IBP': 'ClpC/mecB', '4CU4': 'FERRICHROME-IRON RECEPTOR',
           '9KDF': 'endothelin receptor type-B'}
HUMAN = {'6N60': 'POLR2A', '6N61': 'POLR2A', '6N62': 'POLR2A', '8IBO': 'CLPB',
         '8IBP': 'CLPB', '4CU4': 'FhuA(uptake)', '9KDF': 'EDNRB'}
WEIGHT = {'6N60': 1.0, '6N61': 1.0, '6N62': 1.0, '8IBO': 1.0, '8IBP': 1.0,
          '4CU4': 0.4, '9KDF': 1.0}

pos = []
for pdb, want in PRIMARY.items():
    peps = ent[(ent.pdb_id == pdb) & (ent.length <= 80)]
    tgt = ent[(ent.pdb_id == pdb) & (ent.length > 80) &
              (ent.description.astype(str).str.contains(want, case=False, regex=False))]
    if tgt.empty:
        tgt = ent[(ent.pdb_id == pdb) & (ent.length > 80)]
    if peps.empty or tgt.empty:
        continue
    tgt = tgt.sort_values('length', ascending=False).head(1).iloc[0]
    print(f'  {pdb}: {len(peps)} peptide(s) x {str(tgt.description)[:40]} ({tgt.length} aa)')
    for _, p in peps.iterrows():
        pos.append(dict(pep_seq=p.seq, prot_seq=tgt.seq, label=1,
                        prot_id=f'{pdb}_{want[:18]}'.replace(' ', '_').replace("'", ''),
                        pdb_id=pdb, pep_desc=p.description,
                        human_counterpart=HUMAN.get(pdb, ''),
                        weight=WEIGHT.get(pdb, 1.0), source='pdb_complex'))
posdf = pd.DataFrame(pos).drop_duplicates(subset=['pep_seq', 'prot_seq'])
print(f'\npositives: {len(posdf)} from {posdf.pdb_id.nunique()} structures')
print(posdf[['pdb_id', 'pep_desc', 'human_counterpart']].to_string(index=False))

# decoys: real lasso cores that are NOT in the true set, against the same targets
lib = pd.read_csv(os.path.join(ROOT, 'lassopred_database.csv'))
lib['core'] = lib.Core_Sequence.astype(str).str.upper().map(lambda s: AA.sub('', s))
true_peps = set(posdf.pep_seq)
pool = [s for s in lib.core.unique() if 12 <= len(s) <= 40 and s not in true_peps]
import random
random.seed(7)
neg = []
tgts = posdf[['prot_seq', 'prot_id', 'human_counterpart', 'pdb_id']].drop_duplicates('prot_seq')
for _, t in tgts.iterrows():
    for _ in range(4):
        neg.append(dict(pep_seq=random.choice(pool), prot_seq=t.prot_seq, label=0,
                        prot_id=t.prot_id, pdb_id=t.pdb_id, pep_desc='decoy',
                        human_counterpart=t.human_counterpart, weight=0.3,
                        source='decoy_lasso'))
aux = pd.concat([posdf, pd.DataFrame(neg)], ignore_index=True)
aux.to_csv(os.path.join(ROOT, 'mvp_cpu', 'db_aux_pairs.csv'), index=False, lineterminator='\n')
print(f'\nwrote db_aux_pairs.csv: {len(aux)} rows '
      f'({int((aux.label==1).sum())} pos / {int((aux.label==0).sum())} neg)')

# ---------------- verification scan inputs ----------------
peps = posdf[['pep_seq', 'pdb_id', 'pep_desc']].drop_duplicates('pep_seq').reset_index(drop=True)
pep_rows = []
for i, r in peps.iterrows():
    label = f'P{i+1}'
    if r.pep_desc and str(r.pep_desc).lower() not in ('nan', 'none'):
        label += '_' + re.sub(r'[^A-Za-z0-9]', '', str(r.pep_desc))[:14]
    pep_rows.append(dict(id=label, seq=r.pep_seq, pdb_id=r.pdb_id))
pep_df = pd.DataFrame(pep_rows)
pep_df.to_csv(os.path.join(ROOT, 'mvp_cpu', 'verify_peptides.csv'), index=False,
              lineterminator='\n')
print(f'\nverify peptides: {len(pep_df)}')
for _, r in pep_df.iterrows():
    print(f'  {r.id:22s} {len(r.seq):3d} aa  {r.seq}')

AA_RE = re.compile(r'[^ACDEFGHIKLMNPQRSTVWY]')
scan = {}
name, buf = None, []
for line in open(os.path.join(ROOT, 'mvp_cpu', 'scan_targets.fasta'), encoding='utf-8'):
    line = line.strip()
    if line.startswith('>'):
        if name:
            scan[name] = AA_RE.sub('', ''.join(buf).upper())
        name, buf = line[1:].split()[0], []
    elif line:
        buf.append(line)
if name:
    scan[name] = AA_RE.sub('', ''.join(buf).upper())

targets = {}
for _, r in posdf.drop_duplicates('prot_id').iterrows():
    targets[f'BACT_{r.prot_id}'] = r.prot_seq.replace('X', '')
for want in ['POLR2A', 'CLPB', 'EDNRB']:
    for k, v in scan.items():
        if k.startswith(want):
            targets[f'HUMAN_{want}'] = v
for want in ['MDM2', 'PPIA']:
    for k, v in scan.items():
        if k.startswith(want):
            targets[f'CTRL_{want}'] = v
with open(os.path.join(ROOT, 'mvp_cpu', 'verify_targets.fasta'), 'w',
          encoding='utf-8', newline='\n') as f:
    for k, v in targets.items():
        f.write(f'>{k}\n{v}\n')
print(f'\nverify targets: {len(targets)}')
for k, v in targets.items():
    print(f'  {k:28s} {len(v):5d} aa')

truth = []
for _, r in posdf.iterrows():
    pid = pep_df[pep_df.seq == r.pep_seq].id.iloc[0]
    truth.append(dict(peptide_id=pid, target_id=f'BACT_{r.prot_id}',
                      pdb_id=r.pdb_id, human_counterpart=r.human_counterpart))
pd.DataFrame(truth).drop_duplicates().to_csv(
    os.path.join(ROOT, 'mvp_cpu', 'verify_truth.csv'), index=False, lineterminator='\n')
print(f'\ntruth cells: {len(truth)}')
