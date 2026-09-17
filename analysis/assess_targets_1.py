"""Target assessment, part 1: hard evidence for GUK1 / SSX1 / EXOSC1 as drug targets.

Checks, per protein:
  * what ligands are already crystallised (a bound drug-like molecule is the strongest
    druggability evidence there is)
  * subcellular localisation (a peptide drug must physically reach the target)
  * disease annotation and genetic evidence
  * whether the protein is already inside this project's docking panels
"""
import json
import os
import urllib.request

ROOT = r'D:\deepseek_harness\prp49'
ACC = {'Q16774': 'GUK1', 'Q16384': 'SSX1', 'Q9Y3B2': 'EXOSC1'}


def gql(query):
    req = urllib.request.Request('https://data.rcsb.org/graphql',
                                 data=json.dumps({'query': query}).encode(),
                                 headers={'Content-Type': 'application/json'})
    with urllib.request.urlopen(req, timeout=90) as r:
        return json.loads(r.read()).get('data', {})


print('=== 1. PDB entries and their bound ligands (druggability evidence) ===')
for acc, gene in ACC.items():
    url = (f'https://rest.uniprot.org/uniprotkb/{acc}.json?fields=accession,'
           f'xref_pdb,cc_subcellular_location,cc_disease,cc_function,keyword,'
           f'cc_interaction,ft_binding,ft_act_site')
    req = urllib.request.Request(url, headers={'User-Agent': 'prp49/1.0'})
    with urllib.request.urlopen(req, timeout=60) as r:
        d = json.loads(r.read())
    pdbs = [x['id'] for x in d.get('uniProtKBCrossReferences', []) if x['database'] == 'PDB']
    print(f'\n{acc} ({gene}): {len(pdbs)} PDB entries')
    if pdbs:
        q = '{ entries(entry_ids: [%s]) { rcsb_id struct{title} rcsb_entry_info{resolution_combined} ' \
            'nonpolymer_entities { nonpolymer_comp { chem_comp { id name formula } } } } }' % \
            ','.join(f'"{p}"' for p in pdbs[:12])
        try:
            ents = gql(q).get('entries', [])
            for e in ents:
                if not e:
                    continue
                ligs = []
                for ne in (e.get('nonpolymer_entities') or []):
                    cc = ((ne or {}).get('nonpolymer_comp') or {}).get('chem_comp') or {}
                    if cc.get('id'):
                        ligs.append(f"{cc['id']}({cc.get('name','')[:28]})")
                res = (e.get('rcsb_entry_info', {}).get('resolution_combined') or [None])[0]
                print(f"   {e['rcsb_id']}  {str(res):>6} A  lig: {ligs[:4] if ligs else 'none'}")
                print(f"        {str((e.get('struct') or {}).get('title'))[:96]}")
        except Exception as exc:
            print('   ligand query failed:', type(exc).__name__, exc)

    for c in d.get('comments', []):
        if c.get('commentType') == 'SUBCELLULAR LOCATION':
            locs = [l.get('location', {}).get('value') for l in c.get('subcellularLocations', [])]
            print(f'   LOCALISATION: {[x for x in locs if x]}')
        if c.get('commentType') == 'DISEASE':
            dis = c.get('disease', {})
            print(f"   DISEASE: {dis.get('diseaseId')} - {str(dis.get('diseaseName'))[:60]}")
        if c.get('commentType') == 'INTERACTION':
            n = len(c.get('interactions', []))
            print(f'   INTERACTIONS: {n} annotated partners')

print('\n=== 2. are they already in this project\'s panels? ===')
for f, label in [('mvp_cpu/scan_targets.fasta', '11 scan targets'),
                 ('docking_drugpanel/tasks.csv', '39 druggable-target docking panel'),
                 ('mvp_cpu/affinity_pairs.csv', '564 affinity targets')]:
    p = os.path.join(ROOT, f)
    if not os.path.exists(p):
        print(f'  {label:38s} ({f}) NOT FOUND')
        continue
    txt = open(p, encoding='utf-8', errors='ignore').read()
    hits = [f'{acc}({g})' for acc, g in ACC.items() if acc in txt]
    hits2 = [g for g in ACC.values() if g in txt]
    print(f'  {label:38s} accession hits: {hits or "none"} | gene-name hits: {hits2 or "none"}')

# the drug panel target list, if present
p = os.path.join(ROOT, 'docking_drugpanel', 'tasks.csv')
if os.path.exists(p):
    import pandas as pd
    t = pd.read_csv(p)
    print(f'\n  drug panel columns: {list(t.columns)}')
    if 'target' in t.columns:
        print(f"  panel targets ({t.target.nunique()}): {sorted(t.target.unique())[:20]}")
