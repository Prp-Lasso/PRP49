"""Generate scan inputs: 10 peptides CSV + 10 human targets FASTA."""
import csv
import os

# --- scan_peptides.csv: 10 pre-selected lasso peptides (mature/core sequences) ---
peptides = [
    # id, seq, source note
    ("RES-701-3",   "GNWHGTSPDWFFNYYW",             "lassopred LP_434 core; 9KDF ETB complex verified"),
    ("MccJ25",      "GGAGHVPEYFVGIGTPISFYG",        "lassopred LP_471 core; 6N60/6N61 RNAP"),
    ("Capistruin",  "GTPGFQTPDARVISRFGFN",          "lassopred LP_105 core; 6N61 RNAP"),
    ("Anantin",     "GFIGWGKDIFGHYGG",              "lassopred LP_366 core; NPR-A antagonist"),
    ("Siamycin-I",  "CLGVGSCNDFAGCGYAIVCFW",        "lassopred LP_37 core; 1RPB"),
    ("Chaxapeptin", "GFGSKPLDSFGLNFF",              "lassopred LP_333 core; 2N5C"),
    ("Sphingopyxin-I", "GIEPLGPVDEDQGEHYLFAGGITADD", "lassopred LP_1300 core"),
    ("Lassomycin",  "GLRRLFANQLVGRRNI",             "8IBO chain (training-positive); lit GLRRLFADQLVGRRNI"),
    ("Ubonodin",    "GGDGSIAEYFNRPMHIHDWQIMDSGYYG", "BMRB 30625 mature 28 aa"),
    ("Lariocidin",  "SKKSKPGDGKFGRGVKRG",           "Nature 2025 640:1022, 18 aa; S1-D8 isopeptide"),
]
with open(r'D:\deepseek_harness\prp49\mvp_cpu\scan_peptides.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['id', 'seq', 'source'])
    w.writerows(peptides)
print('scan_peptides.csv written:', len(peptides), 'peptides')

# --- scan_targets.fasta: 10 human targets (11 chains, ITGAV+ITGB3 dimer) ---
targets = [
    ("EDNRB", "P24530"), ("PLXNB1", "O43157"), ("POLR2A", "P24928"),
    ("PPIA", "P62937"), ("FKBP1A", "P62942"), ("NPR1", "P16066"),
    ("C3", "P01024"), ("MDM2", "Q00987"),
    ("ITGAV", "P06756"), ("ITGB3", "P05106"), ("CLPB", "Q9H078"),
]
udir = r'D:\deepseek_harness\prp49\mvp_cpu\uniprot'
out_lines = []
for name, acc in targets:
    raw = open(os.path.join(udir, acc + '.fasta')).read().strip()
    seq = ''.join(line.strip() for line in raw.splitlines() if not line.startswith('>'))
    assert set(seq) <= set('ACDEFGHIKLMNPQRSTVWY'), (name, acc)
    out_lines.append(f'>{name}|{acc}|human')
    for i in range(0, len(seq), 60):
        out_lines.append(seq[i:i + 60])
with open(r'D:\deepseek_harness\prp49\mvp_cpu\scan_targets.fasta', 'w') as f:
    f.write('\n'.join(out_lines) + '\n')
print('scan_targets.fasta written:', len(targets), 'chains')
