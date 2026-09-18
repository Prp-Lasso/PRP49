"""(a) D v2 status for the scheduled check, and (b) where the three accessions sit in
this project - if anywhere.

User's question is about Q16774 / Q16384 / Q9Y3B2 as candidate PEPTIDES, but UniProt
says they are full human proteins (197/188/195 aa). Before concluding anything, check
whether they appear in the project's own data (as targets? as precursors?) and whether
any of them is annotated with a short peptide feature.
"""
import json
import os
import sys
sys.path.insert(0, r'D:\deepseek_harness\prp49\hpc_deploy\ops')
from opslib import connect, run

ROOT = r'D:\deepseek_harness\prp49'

# ---------- (a) D v2 ----------
cli = connect()
print('=== D v2 (62670602) ===')
st, out, err = run(cli, 'date "+%H:%M"; sacct -j 62670602 --format=JobID%12,State%12,Elapsed -X -n 2>/dev/null; '
                        'grep -E "^fold [0-9]+:|WARN" ~/LassoPep/prp49_homdb-62670602.out 2>/dev/null | tail -8; '
                        'if [ -f ~/LassoPep/results/homolog_db_results.json ]; then '
                        'python3 -c "import json;d=json.load(open(\'/dssg/home/acct-clswxl/clswxl-ccmbi1/LassoPep/results/homolog_db_results.json\'));'
                        'print(\'in_domain\',round(d[\'in_domain_auc_mean\'],4),\'cross_domain\',round(d[\'cross_domain_auc_mean\'],4))"; '
                        'else echo "(json not yet)"; fi')
print(out.strip() or err.strip()[:200])

# ---------- (b) the three accessions inside the project ----------
print('\n=== accession search across project data ===')
TERMS = ['Q16774', 'Q16384', 'Q9Y3B2', 'GUK1', 'SSX1', 'EXOSC1', 'Guanylate kinase',
         'CSL4', 'EXOSC1', 'SSX']
exts = ('.csv', '.fasta', '.fa', '.json', '.txt', '.md', '.tsv')
hits = {t: [] for t in TERMS}
for base, dirs, files in os.walk(ROOT):
    if '.git' in base or 'node_modules' in base:
        continue
    for f in files:
        if not f.endswith(exts):
            continue
        p = os.path.join(base, f)
        try:
            if os.path.getsize(p) > 40_000_000:
                continue
            txt = open(p, encoding='utf-8', errors='ignore').read()
        except Exception:
            continue
        for t in TERMS:
            if t in txt:
                hits[t].append(os.path.relpath(p, ROOT))
for t, v in hits.items():
    print(f'  {t:18s} {len(v):3d} file(s)' + (f'  e.g. {v[:2]}' if v else ''))

# ---------- (c) do these proteins contain short peptide features? ----------
print('\n=== peptide-like features in the three proteins ===')
import urllib.request
for acc in ['Q16774', 'Q16384', 'Q9Y3B2']:
    url = f'https://rest.uniprot.org/uniprotkb/{acc}.json'
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'prp49-research/1.0'})
        with urllib.request.urlopen(req, timeout=45) as r:
            d = json.loads(r.read())
    except Exception as exc:
        print(f'  {acc}: {exc}')
        continue
    feats = d.get('features', [])
    interesting = [f for f in feats if f.get('type') in
                   ('Signal', 'Transit peptide', 'Peptide', 'Propeptide')]
    print(f"  {acc} ({d.get('uniProtkbId')}): {len(interesting)} signal/peptide features")
    for f in interesting:
        loc = f.get('location', {})
        print(f"      {f.get('type')}: {loc.get('start',{}).get('value')}-"
              f"{loc.get('end',{}).get('value')} {f.get('description','')[:50]}")
    # any known interaction with lasso peptides / bacteria?
    kws = [k.get('name') for k in d.get('keywords', [])]
    print(f"      keywords: {kws[:8]}")
cli.close()
