"""Where do Q16774 / Q16384 / Q9Y3B2 come from in this project, and in what role?"""
import glob
import os

import pandas as pd

ROOT = r'D:\deepseek_harness\prp49'
TARGETS = ['Q16774', 'Q16384', 'Q9Y3B2']

# locate the files that mention them
files = {}
for pat in ['**/*.csv', '**/*.txt', '**/*.fasta', '**/*.json']:
    for p in glob.glob(os.path.join(ROOT, pat), recursive=True):
        if any(x in p for x in ('github_repo', '.git', 'node_modules')):
            continue
        try:
            head = open(p, encoding='utf-8', errors='ignore').read(400_000)
        except Exception:
            continue
        hit = [t for t in TARGETS if t in head]
        if hit:
            files[p] = (hit, os.path.getsize(p))
print('=== files mentioning the three accessions ===')
for p, (hit, sz) in sorted(files.items(), key=lambda x: -x[1][1]):
    print(f'  {os.path.relpath(p, ROOT):58s} {sz/1024:9.1f} KB  {hit}')

print('\n=== inspect the key files ===')
for name in ['bernett_pairs.csv']:
    cands = [p for p in files if os.path.basename(p) == name]
    for p in cands:
        d = pd.read_csv(p, nrows=5)
        print(f'\n  {os.path.relpath(p, ROOT)}')
        print(f'    columns: {list(d.columns)}')
        print(f'    rows   : {sum(1 for _ in open(p, encoding="utf-8", errors="ignore")) - 1}')
        full = pd.read_csv(p)
        for t in TARGETS:
            sub = full[full.astype(str).apply(lambda r: r.str.contains(t).any(), axis=1)]
            if len(sub):
                print(f'    {t}: {len(sub)} rows')
                print('      ', sub.head(2).to_dict('records'))
            else:
                print(f'    {t}: (not in this file)')

for name in ['Intra1_neg_rr.txt', 'Intra2_neg_rr.txt']:
    cands = [p for p in files if os.path.basename(p) == name]
    for p in cands:
        print(f'\n  {os.path.relpath(p, ROOT)}')
        lines = open(p, encoding='utf-8', errors='ignore').read().splitlines()
        print(f'    lines: {len(lines)} | first 3: {lines[:3]}')
        for t in TARGETS:
            n = sum(1 for l in lines if t in l)
            print(f'    {t}: appears in {n} lines')

# also check the fasta role
for name in ['human_swissprot_oneliner.fasta', 'data_human_proteome.fasta']:
    for p in [x for x in files if os.path.basename(x) == name]:
        print(f'\n  {os.path.relpath(p, ROOT)}  ({os.path.getsize(p)/1e6:.1f} MB)')
        txt = open(p, encoding='utf-8', errors='ignore').read()
        for t in TARGETS:
            i = txt.find(t)
            if i >= 0:
                seg = txt[max(0, i - 90):i + 150].replace('\n', ' | ')
                print(f'    {t}: ...{seg[:230]}...')
