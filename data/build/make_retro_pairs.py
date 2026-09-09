"""Build retro_pairs.csv: held-out known interaction pairs for retrospective validation."""
import os
import pandas as pd

root = r'D:\deepseek_harness\prp49'
exp = pd.read_csv(os.path.join(root, 'mvp_cpu', 'expanded_pairs.csv'))
# structural/literature anchor pairs (one row per peptide-target)
anchors = [
    ('RES701-3', '9KDF:R'),      # cryo-EM ETB complex
    ('PB1m7', '7VF3'),           # X-ray Plexin B1 complex
    ('MccJ25', '6N60'),          # cryo-EM RNAP complex
    ('Capistruin', '6N61'),      # cryo-EM RNAP complex
    ('Lassomycin', '8IBO:A'),    # ClpC1 NTD complex
    ('Anantin', 'NPR1'),         # biochemical antagonist (no structure)
]
rows = []
for pid, pdb in anchors:
    cand = exp[exp['pep_id'] == pid]
    if pdb.startswith('NPR'):
        cand = exp[exp['prot_id'].astype(str).str.startswith('P16066')]
    if len(cand):
        r = cand.iloc[0]
        rows.append({'pep_id': pid, 'pep_seq': r['pep_seq'],
                     'prot_seq': r['prot_seq'], 'prot_id': r['prot_id']})
    else:
        print('WARN: not found', pid, pdb)
out = pd.DataFrame(rows)
out.to_csv(os.path.join(root, 'mvp_cpu', 'retro_pairs.csv'), index=False)
print('retro_pairs.csv:', len(out), 'anchor pairs')
print(out[['pep_id', 'prot_id']].to_string(index=False))
