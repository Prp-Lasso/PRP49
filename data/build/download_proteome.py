"""Download UniProt human reference proteome (UP000005640) to data/human_proteome.fasta."""
import os
import time
import urllib.request

out = r'D:\deepseek_harness\prp49\data_human_proteome.fasta'
url = ('https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000005640'
       '+AND+reviewed:true&format=fasta&compressed=false')
# reviewed-only (Swiss-Prot, ~20k) keeps the file small; full proteome incl.
# TrEMBL (~80k) is available by dropping +AND+reviewed:true

for attempt in range(4):
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'prp49-prep/1.0'})
        with urllib.request.urlopen(req, timeout=300) as r, open(out, 'wb') as f:
            while True:
                chunk = r.read(1 << 20)
                if not chunk:
                    break
                f.write(chunk)
        size = os.path.getsize(out) / 1e6
        n = sum(1 for line in open(out, encoding='utf-8', errors='ignore') if line.startswith('>'))
        print(f'OK {out}: {size:.1f} MB, {n} entries')
        break
    except Exception as e:
        print(f'attempt {attempt} failed: {e}')
        time.sleep(10 * (attempt + 1))
else:
    print('DOWNLOAD FAILED')
