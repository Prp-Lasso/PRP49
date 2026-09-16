"""Route D' step 2 (v3): per-HMM search loop.

v2 called pyhmmer.hmmsearch(hmms, block) with all 53 HMMs at once and silently
returned nothing - yet a single HMM against a single sequence works perfectly
(PPIA x PF00160 -> E=4.2e-53). Rather than rely on the batched path, each HMM is
searched separately with cpus=1, and every HMM reports its own hit count so a
silent failure cannot hide again.

Usage: python search_domains_v3.py [--control]
  --control : search the 11 human targets themselves (must recover their own domains)
  default   : search all Propedia receptors
"""
import json
import os
import re
import sys

import pandas as pd
import pyhmmer

ROOT = r'D:\deepseek_harness\prp49'
DOM = os.path.join(ROOT, 'homolog_domains')
AA = re.compile(r'[^ACDEFGHIKLMNPQRSTVWY]')
CONTROL = '--control' in sys.argv


def _s(v):
    return v.decode() if isinstance(v, (bytes, bytearray)) else str(v)


def read_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path, encoding='utf-8'):
        line = line.strip()
        if line.startswith('>'):
            if name:
                seqs[name] = AA.sub('', ''.join(buf).upper())
            name, buf = line[1:].split()[0], []
        elif line:
            buf.append(line)
    if name:
        seqs[name] = AA.sub('', ''.join(buf).upper())
    return seqs


tdom = json.load(open(os.path.join(DOM, 'target_domains.json')))
target_domains = {k: {pf for pf, _ in v} for k, v in tdom.items()}

hmms = []
for f in sorted(os.listdir(DOM)):
    if f.endswith('.hmm'):
        with pyhmmer.plan7.HMMFile(os.path.join(DOM, f)) as hf:
            hmms.append(hf.read())
print(f'loaded {len(hmms)} HMMs')

if CONTROL:
    tgt = read_fasta(os.path.join(ROOT, 'mvp_cpu', 'scan_targets.fasta'))
    items = list(tgt.items())
    label = 'human targets (CONTROL)'
else:
    prop = pd.read_csv(os.path.join(ROOT, 'docking_propedia', 'all_scores.csv'))
    recs = prop[['rec_seq']].drop_duplicates().reset_index(drop=True)
    recs['rec_seq'] = recs.rec_seq.astype(str).map(lambda s: AA.sub('', s.upper()))
    recs = recs[recs.rec_seq.str.len() >= 30].reset_index(drop=True)
    items = [(str(i), s) for i, s in enumerate(recs.rec_seq)]
    label = f'Propedia receptors'
print(f'scanning {len(items)} sequences: {label}')

alphabet = pyhmmer.easel.Alphabet.amino()
idx2name, dig = {}, []
for i, (k, v) in enumerate(items):
    idx2name[str(i)] = k
    dig.append(pyhmmer.easel.TextSequence(name=str(i).encode(), sequence=v).digitize(alphabet))
block = pyhmmer.easel.DigitalSequenceBlock(alphabet, dig)

hits = {}
for n, hmm in enumerate(hmms, start=1):
    acc = hmm.accession or hmm.name
    pf = _s(acc).split('.')[0]
    n_hit = 0
    try:
        res = list(pyhmmer.hmmsearch([hmm], block, cpus=1, incE=1e-3, incdomE=1e-3))
    except Exception as exc:
        print(f'  [{n:2d}/{len(hmms)}] {pf}: search failed {type(exc).__name__}: {exc}')
        continue
    for hh in res:
        for hit in hh:
            i = _s(hit.name)
            d = hits.setdefault(i, {})
            ev, bits = float(hit.evalue), float(hit.score)
            prev = d.get(pf)
            if prev is None or ev < prev[0]:
                d[pf] = (ev, bits, len(hit.domains))
            n_hit += 1
    if n_hit or n % 10 == 0:
        print(f'  [{n:2d}/{len(hmms)}] {pf:10s}: {n_hit} hits '
              f'({len(hits)} sequences hit so far)', flush=True)

print(f'\ntotal sequences with >=1 domain hit: {len(hits)} / {len(items)}')

if CONTROL:
    ok = 0
    for k, doms in tdom.items():
        expected = {pf for pf, _ in doms}
        i = str([j for j, (nm, _) in enumerate(items) if nm == k][0])
        got = set(hits.get(i, {}))
        m = len(expected & got)
        ok += bool(m)
        print(f'  {k:9s} expected {len(expected):2d} | matched {m:2d} '
              f'{"OK" if m else "MISS"}')
    print(f'\ncontrol: {ok}/{len(tdom)} targets recovered -> '
          f'{"pipeline verified" if ok >= len(tdom)*0.7 else "PIPELINE STILL BROKEN"}')
    sys.exit(0 if ok >= len(tdom) * 0.7 else 1)

rows = []
for name, dset in target_domains.items():
    if not dset:
        continue
    for i, hd in hits.items():
        shared = dset & set(hd)
        if not shared:
            continue
        cov = len(shared) / len(dset)
        best_bits = max(hd[pf][1] for pf in shared)
        w = float(min(1.0, cov * (0.5 + 0.5 * min(1.0, best_bits / 200.0))))
        rows.append(dict(human_target=name, aux_rec_seq=items[int(i)][1],
                         n_shared=len(shared), n_target_domains=len(dset),
                         coverage=round(cov, 3), domains=','.join(sorted(shared)),
                         best_bits=round(best_bits, 1), weight=round(max(w, 0.05), 3)))

h = pd.DataFrame(rows)
if len(h):
    h = h.sort_values(['human_target', 'weight'], ascending=[True, False])
out = os.path.join(ROOT, 'mvp_cpu', 'homolog_domains.csv')
h.to_csv(out, index=False, lineterminator='\n')
print(f'\nwrote {out}: {len(h)} rows')
if len(h):
    print('\nper-target receptors:')
    print(h.groupby('human_target').agg(receptors=('aux_rec_seq', 'nunique'),
                                        max_w=('weight', 'max'),
                                        med_cov=('coverage', 'median')).round(3).to_string())
    print('\ntop correspondences:')
    print(h.head(12)[['human_target', 'n_shared', 'n_target_domains', 'coverage',
                      'best_bits', 'weight', 'domains']].to_string(index=False))
