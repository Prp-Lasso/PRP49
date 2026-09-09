"""Build a larger positive-pair set from RCSB human peptide complexes.

Pipeline:
  1. Query polymer_entities: human + 5-40aa (paginated).
  2. For each hit entry: fetch FASTA (all chains), find peptide chain(s) 5-40aa
     and protein chains (>40aa).
  3. Filter out MHC/antibody complexes (description keywords) to keep
     physiologically-relevant peptide-protein interfaces.
  4. Sample N entries; save pairs.
"""
import sys
import json
import time
import random
import urllib.request
import urllib.parse
from collections import defaultdict

sys.stdout.reconfigure(encoding='utf-8')

AA1 = set('ACDEFGHIKLMNPQRSTVWY')
EXCLUDE_KW = ['MHC', 'HLA', 'class I histocompatibility', 'class II histocompatibility',
              'immunoglobulin', 'antibody', 'FAB', 'Fab', 'T-cell receptor',
              'T cell receptor', 'B2M', 'beta-2-microglobulin', 'MHC CLASS', 'antigen-binding']


def rcsb_query(start, rows=1000):
    q = {
        "query": {
            "type": "group", "logical_operator": "and", "nodes": [
                {"type": "terminal", "service": "text",
                 "parameters": {"attribute": "rcsb_entity_source_organism.ncbi_scientific_name",
                                "operator": "exact_match", "value": "Homo sapiens"}},
                {"type": "terminal", "service": "text",
                 "parameters": {"attribute": "entity_poly.rcsb_sample_sequence_length",
                                "operator": "range",
                                "value": {"from": 5, "to": 40,
                                          "include_lower": True, "include_upper": True}}},
            ]
        },
        "return_type": "polymer_entity",
        "request_options": {
            "paginate": {"start": start, "rows": rows},
            "results_content_type": ["experimental"]
        }
    }
    req = urllib.request.Request('https://search.rcsb.org/rcsbsearch/v2/query',
                                 data=json.dumps(q).encode(),
                                 headers={'Content-Type': 'application/json', 'User-Agent': 'Mozilla/5.0'})
    r = urllib.request.urlopen(req, timeout=120)
    return json.loads(r.read().decode())


def entry_fasta(pdb):
    url = f'https://www.rcsb.org/fasta/entry/{pdb}'
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    data = urllib.request.urlopen(req, timeout=60).read().decode('utf-8', errors='replace')
    seqs = {}
    cur = None
    for line in data.splitlines():
        if line.startswith('>'):
            cur = line[1:].strip()
            seqs[cur] = ''
        elif cur is not None:
            seqs[cur] += line.strip()
    return seqs


def entry_title(pdb):
    try:
        req = urllib.request.Request(f'https://data.rcsb.org/rest/v1/core/entry/{pdb}',
                                     headers={'User-Agent': 'Mozilla/5.0'})
        d = json.loads(urllib.request.urlopen(req, timeout=30).read().decode())
        return d.get('struct', {}).get('title', '') or ''
    except Exception:
        return ''


def main(max_entries=1500, max_hits=9000, max_per_protein=15):
    hits = []
    for start in range(0, max_hits, 1000):
        d = rcsb_query(start, 1000)
        total = d.get('total_count')
        batch = [e['identifier'] for e in d.get('result_set', [])]
        hits.extend(batch)
        print(f'fetched {len(hits)}/{total} hits', flush=True)
        if len(batch) < 1000:
            break
    print('total hits collected:', len(hits))

    # entry -> peptide chains
    entries = sorted(set(h.split('_')[0] for h in hits))
    random.seed(42)
    random.shuffle(entries)
    print('unique entries:', len(entries))

    pairs = []
    prot_count = defaultdict(int)
    titles_checked = 0
    for pdb in entries:
        if len(pairs) >= max_entries:
            break
        try:
            seqs = entry_fasta(pdb)
        except Exception:
            continue
        # parse header description for filtering
        descs = []
        short = []  # (chain, seq)
        longs = []  # (chain, seq)
        header_meta = []
        for hdr, seq in seqs.items():
            # FASTA header like: "4HHB_1|Chain A|HEMOGLOBIN|Homo sapiens"
            parts = hdr.split('|')
            desc = parts[2] if len(parts) > 2 else hdr
            descs.append(desc)
            header_meta.append((hdr, seq, desc))
        title = entry_title(pdb)
        titles_checked += 1
        blob = ' '.join(descs) + ' ' + title
        if any(kw.lower() in blob.lower() for kw in EXCLUDE_KW):
            continue
        # classify chains
        for hdr, seq, desc in header_meta:
            s = seq.replace('*', '').replace('-', '')
            if not s or not all(a in AA1 for a in s):
                continue
            if 5 <= len(s) <= 40:
                short.append((hdr.split('|')[1] if '|' in hdr else hdr, s, desc))
            elif len(s) > 40:
                longs.append((hdr.split('|')[1] if '|' in hdr else hdr, s, desc))
        if not short or not longs:
            continue
        # take the longest protein chain as receptor; all short chains as peptides
        rec = max(longs, key=lambda t: len(t[1]))
        if prot_count[rec[1]] >= max_per_protein:
            continue
        for ch, s, d in short:
            pairs.append(dict(pdb=pdb, pep_id=f'{pdb}:{ch}', pep_seq=s,
                              prot_id=rec[0], prot_seq=rec[1], title=title[:90]))
            prot_count[rec[1]] += 1
            if prot_count[rec[1]] >= max_per_protein:
                break
        print(f'{pdb}: {len(short)} peptide chains -> {len(pairs)} pairs (title: {title[:60]})', flush=True)

    with open('mvp_cpu/rcsb_pairs.json', 'w', encoding='utf-8') as f:
        json.dump(pairs, f, indent=1, ensure_ascii=False)
    print(f'\nsaved {len(pairs)} pairs from RCSB (entries checked: {titles_checked})')


if __name__ == '__main__':
    main()
