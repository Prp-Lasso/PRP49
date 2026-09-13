"""Locate druggable pockets for the 41-target panel using RCSB co-crystal ligands.

Strategy (better than blind docking on AlphaFold models):
  for each UniProt accession, find a human structure solved with a small-molecule
  ligand (resolution <= 3.0 A), take the largest drug-like ligand as the pocket
  marker, and derive the docking box from its coordinates.
  Targets without a usable holo structure fall back to the AlphaFold model with a
  whole-protein box (flagged in the output).

Outputs (docking_drugpanel/):
  poc/*.pdb          receptor chains with the ligand stripped
  poc_meta.csv       target, source(holo|af), pdb_id, ligand, center, size
  tasks_drugpanel.csv docking tasks ready for job_dock_drugpanel.sh
"""
import csv
import json
import os
import time
import urllib.request

ROOT = r'D:\deepseek_harness\prp49'
OUT = os.path.join(ROOT, 'docking_drugpanel')
os.makedirs(os.path.join(OUT, 'poc'), exist_ok=True)

LIG_IGNORE = {'HOH', 'WAT', 'DOD', 'NA', 'K', 'CL', 'MG', 'CA', 'ZN', 'MN', 'FE', 'CU',
              'SO4', 'PO4', 'GOL', 'EDO', 'PEG', 'ACT', 'DMS', 'TRS', 'MES', 'IOD', 'BR',
              'FMT', 'CIT', 'EPE', 'MPD', 'BME', 'NAG', 'BMA', 'MAN', 'FUC', 'GAL', 'SIA',
              # lipids / detergents / cryo additives: NOT drug-pocket markers
              'CLR', 'LBN', 'OLA', 'PLM', 'STE', 'MYR', 'PC1', 'LPP', 'DGA', 'TCH', 'CHL',
              'DDM', 'LMT', 'LMG', 'C8E', 'BOG', 'PE5', 'PGE', 'PG4', '1PE', 'P6G', 'XTA',
              'TYS', 'B3P', 'FT8'}


def fetch(url, timeout=30):
    req = urllib.request.Request(url, headers={'User-Agent': 'prp49-pocket-finder/1.0'})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def search_holo(uniprot, max_res=3.0):
    """Return candidate PDB ids with a drug-like ligand, newest/highest-res first."""
    q = {
        "query": {"type": "group", "logical_operator": "and", "nodes": [
            {"type": "terminal", "service": "text",
             "parameters": {"attribute": "rcsb_polymer_entity_container_identifiers."
                                        "reference_sequence_identifiers.database_accession",
                            "operator": "exact_match", "value": uniprot}},
            {"type": "terminal", "service": "text",
             "parameters": {"attribute": "rcsb_entry_info.resolution_combined",
                            "operator": "less_or_equal", "value": max_res}},
            {"type": "terminal", "service": "text",
             "parameters": {"attribute": "rcsb_entry_info.nonpolymer_entity_count",
                            "operator": "greater_or_equal", "value": 1}},
        ]},
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": 8},
                            "sort": [{"sort_by": "rcsb_entry_info.resolution_combined",
                                      "direction": "asc"}]},
    }
    try:
        raw = fetch('https://search.rcsb.org/rcsbsearch/v2/query?json=' +
                    urllib.parse.quote(json.dumps(q)))
        data = json.loads(raw)
        return [hit['identifier'] for hit in data.get('result_set', [])]
    except Exception as exc:
        print(f'    search failed: {exc}')
        return []


import urllib.parse  # noqa: E402  (used inside search_holo)


def get_pdb_text(pdb_id):
    for url in (f'https://files.rcsb.org/download/{pdb_id}.pdb',
                f'https://files.rcsb.org/download/{pdb_id}.cif'):
        try:
            txt = fetch(url).decode('utf-8', 'replace')
            if url.endswith('.cif'):
                return None      # keep it simple: PDB format only
            return txt
        except Exception:
            continue
    return None


def biggest_ligand(pdb_text):
    """Return (resname, chain, center, n_atoms) for the largest drug-like HETATM group."""
    groups = {}
    for line in pdb_text.splitlines():
        if not line.startswith('HETATM'):
            continue
        resn = line[17:20].strip()
        if resn in LIG_IGNORE:
            continue
        key = (resn, line[21], line[22:27])
        try:
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue
        groups.setdefault(key, []).append(xyz)
    if not groups:
        return None
    key, atoms = max(groups.items(), key=lambda kv: len(kv[1]))
    if len(atoms) < 8:          # too small to be a drug-like pocket marker
        return None
    cx = sum(a[0] for a in atoms) / len(atoms)
    cy = sum(a[1] for a in atoms) / len(atoms)
    cz = sum(a[2] for a in atoms) / len(atoms)
    return key[0], key[1], (cx, cy, cz), len(atoms)


def strip_and_save(pdb_text, path, keep_ligand=False, lig_resn=None):
    """Keep protein ATOM records; optionally drop waters.

    NOTE: PDB record names occupy columns 1-6 padded with spaces ("ATOM  "),
    so the slice must be stripped before comparison - a raw `line[:6] == 'ATOM'`
    test silently discards every atom (this bug produced 2-byte receptor files).
    """
    out = []
    for line in pdb_text.splitlines():
        rec = line[:6].strip()
        if rec == 'ATOM':
            out.append(line)
        elif rec == 'HETATM' and keep_ligand:
            out.append(line)
    with open(path, 'w', newline='\n') as f:
        f.write('\n'.join(out) + '\n')
    return len(out)


def main():
    panel = list(csv.DictReader(open(os.path.join(OUT, 'target_panel.csv'))))
    print(f'panel: {len(panel)} targets')
    meta, tasks = [], []
    tid = 5000
    for row in panel:
        gene, acc = row['gene'], row['uniprot']
        print(f'  {gene} ({acc}) ...', flush=True)
        try:
            pdb_id = lig = None
            center = size = None
            source = 'af'
            for cand in search_holo(acc):
                txt = get_pdb_text(cand)
                if not txt:
                    continue
                hit = biggest_ligand(txt)
                if not hit:
                    continue
                pdb_id, lig, center, n_atoms = cand, hit[0], hit[2], hit[3]
                # box: ligand extent + 12 A padding, clamped to a sane range
                xs, ys, zs = [], [], []
                for line in txt.splitlines():
                    if line.startswith('HETATM') and line[17:20].strip() == lig:
                        xs.append(float(line[30:38])); ys.append(float(line[38:46])); zs.append(float(line[46:54]))
                ext = (max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs)) if xs else (0, 0, 0)
                # min edge 38 A: large lasso peptides need room (smaller boxes gave
                # spurious positive scores in the earlier 110-pair matrix)
                size = tuple(min(max(e + 12, 38), 60) for e in ext)
                strip_and_save(txt, os.path.join(OUT, 'poc', f'{gene}.pdb'))
                source = 'holo'
                print(f'    -> holo {cand} ligand {lig} ({n_atoms} atoms) box {size}')
                break
        except Exception as exc:
            print(f'    !! error: {type(exc).__name__}: {exc}')
            source, pdb_id, lig, center, size = 'af', None, None, None, None
        if source == 'af':
            af = os.path.join(OUT, 'af', f'{gene}.pdb')
            if os.path.exists(af):
                import shutil
                shutil.copy2(af, os.path.join(OUT, 'poc', f'{gene}.pdb'))
                print('    -> AF fallback (whole-protein box)')
            else:
                print('    -> skipped (no structure)')
                continue
        meta.append(dict(target=gene, uniprot=acc, source=source, pdb_id=pdb_id or '',
                         ligand=lig or '', cx=round(center[0], 1) if center else '',
                         cy=round(center[1], 1) if center else '',
                         cz=round(center[2], 1) if center else '',
                         sx=round(size[0]) if size else 90, sy=round(size[1]) if size else 90,
                         sz=round(size[2]) if size else 90))

    with open(os.path.join(OUT, 'poc_meta.csv'), 'w', newline='\n') as f:
        w = csv.DictWriter(f, fieldnames=list(meta[0].keys()))
        w.writeheader(); w.writerows(meta)

    # tasks: 13 lasso peptides x every panel target
    peps = list(csv.DictReader(open(os.path.join(ROOT, 'mvp_cpu', 'scan_peptides_ext.csv'))))
    rows = []
    for m in meta:
        for p in peps:
            rows.append(dict(id=tid, pep=f'dp_{p["id"]}', rec=f'dp_{m["target"]}',
                             cx=m['cx'], cy=m['cy'], cz=m['cz'],
                             sx=m['sx'], sy=m['sy'], sz=m['sz']))
            tid += 1
    with open(os.path.join(OUT, 'tasks_drugpanel.csv'), 'w', newline='\n') as f:
        w = csv.writer(f)
        w.writerow(['id', 'pep', 'rec', 'cx', 'cy', 'cz', 'sx', 'sy', 'sz'])
        for r in rows:
            w.writerow([r['id'], r['pep'], r['rec'], r['cx'], r['cy'], r['cz'], r['sx'], r['sy'], r['sz']])

    holo = sum(1 for m in meta if m['source'] == 'holo')
    print(f'\nmeta written: {len(meta)} targets ({holo} with co-crystal pocket, {len(meta)-holo} AF fallback)')
    print(f'tasks written: {len(rows)} pairs -> tasks_drugpanel.csv')


if __name__ == '__main__':
    main()
