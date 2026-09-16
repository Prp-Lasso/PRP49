"""Route D' step 1: Pfam domains for the 11 human targets, then fetch their HMMs.

Sequence identity between human POLR2A and the bacterial RNAP beta-prime subunit is
under 15%, so no BLAST-style search finds it - yet both carry the same Pfam domains
(RNA_pol_Rpb1_*), which is exactly the cross-species correspondence route D needs.
Domains are therefore the right level at which to link the human targets to the
bacterial receptors in Propedia.
"""
import json
import os
import time
import urllib.request

ROOT = r'D:\deepseek_harness\prp49'
OUT = os.path.join(ROOT, 'homolog_domains')
os.makedirs(OUT, exist_ok=True)

ACC = {'EDNRB': 'P24530', 'PLXNB1': 'O43157', 'POLR2A': 'P24928', 'PPIA': 'P62937',
       'FKBP1A': 'P62942', 'NPR1': 'P16066', 'C3': 'P01024', 'MDM2': 'Q00987',
       'ITGAV': 'P06756', 'ITGB3': 'P05106', 'CLPB': 'Q9H078'}


def get(url, timeout=60):
    req = urllib.request.Request(url, headers={'User-Agent': 'prp49-research/1.0'})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


domains = {}
for name, acc in ACC.items():
    url = (f'https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/{acc}/'
           f'?page_size=200')
    try:
        data = json.loads(get(url))
    except Exception as exc:
        print(f'{name} ({acc}): query failed - {exc}')
        continue
    found = []
    for item in data.get('results', []):
        md = item.get('metadata', {})
        pf = md.get('accession')
        nm = md.get('name')
        if pf and pf.startswith('PF'):
            found.append((pf, nm))
    found = sorted(set(found))
    domains[name] = found
    print(f'{name} ({acc}): {len(found)} Pfam domains')
    for pf, nm in found[:8]:
        print(f'    {pf}  {nm}')
    time.sleep(0.3)

json.dump(domains, open(os.path.join(OUT, 'target_domains.json'), 'w'), indent=2)

# --- download the HMM for every distinct domain ---
allpf = sorted({pf for v in domains.values() for pf, _ in v})
print(f'\ndistinct domains to fetch: {len(allpf)}')
ok, fail = 0, []
for pf in allpf:
    dst = os.path.join(OUT, f'{pf}.hmm')
    if os.path.exists(dst) and os.path.getsize(dst) > 500:
        ok += 1
        continue
    url = f'https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/{pf}/?annotation=hmm'
    try:
        blob = get(url, timeout=90)
        if blob.startswith(b'HMMER3') or b'HMMER3' in blob[:200]:
            open(dst, 'wb').write(blob)
            ok += 1
        else:
            # some responses are gzipped
            import gzip
            try:
                blob = gzip.decompress(blob)
                if b'HMMER3' in blob[:200]:
                    open(dst, 'wb').write(blob)
                    ok += 1
                else:
                    fail.append((pf, 'not HMMER3 after gunzip'))
            except Exception:
                fail.append((pf, 'not HMMER3'))
    except Exception as exc:
        fail.append((pf, str(exc)[:60]))
    time.sleep(0.2)

print(f'HMM downloaded: {ok}/{len(allpf)}')
if fail:
    print('failures:', fail[:10])
print('sizes (KB):', {f: round(os.path.getsize(os.path.join(OUT, f))/1024, 1)
                      for f in sorted(os.listdir(OUT)) if f.endswith('.hmm')})
