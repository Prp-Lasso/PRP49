"""Resolve UniProt accessions for the 11 lasso-panel targets and fetch AFDB v6 models.

Why: the lasso panel currently docks into EXPERIMENTAL PDB structures while the
drug panel uses AlphaFold models - mixed provenance. This unifies the lasso panel
onto the same source (AFDB v6), matching what the drug panel already uses.

Note on "AlphaFold batch API": the AlphaFold Server has no official public API
(third-party wrappers exist but need an interactive Google login and daily quota).
For canonical UniProt proteins we do not need prediction at all - AFDB already
serves the models, which is also what keeps provenance identical across panels.
"""
import csv
import json
import os
import time
import urllib.parse
import urllib.request

ROOT = r'D:\deepseek_harness\prp49'
OUT = os.path.join(ROOT, 'docking_af')
os.makedirs(os.path.join(OUT, 'af'), exist_ok=True)

# lasso-panel targets -> gene symbol (ITGAVB3 is a heterodimer: two chains)
TARGETS = [('EDNRB', 'EDNRB'), ('PLXNB1', 'PLXNB1'), ('POLR2A', 'POLR2A'),
           ('PPIA', 'PPIA'), ('FKBP1A', 'FKBP1A'), ('NPR1', 'NPR1'),
           ('C3', 'C3'), ('MDM2', 'MDM2'), ('ITGAV', 'ITGAV'), ('ITGB3', 'ITGB3'),
           ('CLPB', 'CLPB')]

# known accessions (verified earlier in the project) - used as a fast path
KNOWN = {'EDNRB': 'P24530', 'NPR1': 'P16066', 'C3': 'P01024', 'MDM2': 'Q00987',
         'PPIA': 'P62937', 'FKBP1A': 'P62942', 'POLR2A': 'P24928'}


def uniprot_for(gene):
    if gene in KNOWN:
        return KNOWN[gene], 'known'
    q = (f'https://rest.uniprot.org/uniprotkb/search?query='
         f'(gene_exact:{gene})+AND+(organism_id:9606)+AND+(reviewed:true)'
         f'&fields=accession,protein_name,length&format=json&size=5')
    try:
        with urllib.request.urlopen(q, timeout=30) as r:
            data = json.loads(r.read())
        results = data.get('results', [])
        if results:
            acc = results[0]['primaryAccession']
            ln = results[0].get('sequence', {}).get('length')
            return acc, f'query({len(results)} hits, len {ln})'
    except Exception as exc:
        return None, f'query failed: {exc}'
    return None, 'no reviewed human hit'


def fetch_afdb(acc, gene):
    url = f'https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v6.pdb'
    dst = os.path.join(OUT, 'af', f'{gene}.pdb')
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'prp49/1.0'})
        with urllib.request.urlopen(req, timeout=60) as r:
            data = r.read()
        if len(data) < 1000:
            return None, f'suspiciously small ({len(data)} B)'
        open(dst, 'wb').write(data)
        n_atoms = data.count(b'\nATOM')
        return dst, f'{len(data)/1024:.0f} KB, {n_atoms} atoms'
    except Exception as exc:
        return None, f'download failed: {exc}'


rows = []
for gene, symbol in TARGETS:
    acc, how = uniprot_for(symbol)
    print(f'{symbol:9s} -> {acc or "??":10s} ({how})', end='  ')
    if not acc:
        rows.append(dict(gene=symbol, uniprot='', status='NO ACCESSION', detail=how))
        print()
        continue
    path, detail = fetch_afdb(acc, symbol)
    status = 'ok' if path else 'FAILED'
    print(f'{status}: {detail}')
    rows.append(dict(gene=symbol, uniprot=acc, status=status, detail=detail,
                     url=f'https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v6.pdb'))
    time.sleep(0.4)

with open(os.path.join(OUT, 'lasso_af_urls.txt'), 'w', newline='\n') as f:
    for r in rows:
        if r['uniprot']:
            f.write(f"{r['gene']} {r['uniprot']} {r['url']}\n")
with open(os.path.join(OUT, 'fetch_report.csv'), 'w', newline='\n') as f:
    w = csv.DictWriter(f, fieldnames=['gene', 'uniprot', 'status', 'detail', 'url'])
    w.writeheader(); w.writerows(rows)

ok = sum(1 for r in rows if r['status'] == 'ok')
print(f'\n=== summary: {ok}/{len(rows)} AFDB models fetched -> {OUT}/af ===')
for r in rows:
    print(f"  {r['gene']:9s} {r['uniprot']:10s} {r['status']:10s} {r['detail']}")
