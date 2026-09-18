"""Merge a teammate's evidence table into the screening shortlist.

Teammates fill `templates/teammate_evidence_template.csv`; this script joins it onto the
project's candidate table so external evidence appears as first-class columns rather than
being pasted into a chat.

Usage:
  python scratch/merge_teammate_evidence.py <filled.csv> [--candidates results/candidates_screen50.csv]
                                             [--out results/candidates_merged.csv]

Rules enforced:
  * peptide/target names must match this project's spelling (case-insensitive match is
    attempted, but any row that matches nothing is REPORTED, not silently dropped)
  * confidence is restricted to high/medium/low
  * evidence_type is restricted to experiment/literature/docking/prediction/assay
  * the merge is left-join on the candidate table: a row that names an unknown pair is
    printed as UNMATCHED so it can be fixed or added manually
"""
import argparse
import os
import sys

import pandas as pd

CONF = {'high', 'medium', 'low'}
KIND = {'experiment', 'literature', 'docking', 'prediction', 'assay'}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('evidence_csv')
    ap.add_argument('--candidates', default='results/candidates_screen50.csv')
    ap.add_argument('--out', default='results/candidates_merged.csv')
    args = ap.parse_args()

    ev = pd.read_csv(args.evidence_csv, comment='#')
    ev.columns = [c.strip().lower() for c in ev.columns]
    need = {'peptide', 'target', 'evidence_type', 'source', 'confidence', 'note'}
    missing = need - set(ev.columns)
    if missing:
        print(f'ERROR: evidence file is missing columns: {sorted(missing)}')
        sys.exit(1)
    ev = ev.dropna(subset=['peptide', 'target'])
    print(f'evidence rows: {len(ev)}')

    bad_kind = ev[~ev.evidence_type.str.lower().isin(KIND)]
    bad_conf = ev[~ev.confidence.astype(str).str.lower().isin(CONF)]
    for tag, bad, allowed in [('evidence_type', bad_kind, KIND), ('confidence', bad_conf, CONF)]:
        if len(bad):
            print(f'WARNING: {len(bad)} rows have an invalid {tag}; allowed = {sorted(allowed)}')
            print(bad.head(5).to_string(index=False))

    cand = pd.read_csv(args.candidates)
    pcol = next((c for c in ['peptide_id', 'peptide'] if c in cand.columns), None)
    tcol = next((c for c in ['target_id', 'target'] if c in cand.columns), None)
    if pcol is None or tcol is None:
        print(f'ERROR: cannot find peptide/target columns in {args.candidates}')
        sys.exit(1)
    print(f'candidate table: {len(cand)} rows, keyed on {pcol}/{tcol}')

    ev['_p'] = ev.peptide.astype(str).str.strip().str.lower()
    ev['_t'] = ev.target.astype(str).str.strip().str.lower()
    cand['_p'] = cand[pcol].astype(str).str.strip().str.lower()
    cand['_t'] = cand[tcol].astype(str).str.strip().str.lower()

    agg = (ev.groupby(['_p', '_t'])
             .agg(ext_evidence=('evidence_type', lambda s: ';'.join(sorted(set(s)))),
                  ext_source=('source', lambda s: ';'.join(sorted(set(map(str, s))))),
                  ext_confidence=('confidence', lambda s: ';'.join(sorted(set(map(str, s))))),
                  ext_note=('note', lambda s: ' | '.join(map(str, s.dropna()))))
             .reset_index())

    merged = cand.merge(agg, on=['_p', '_t'], how='left')
    matched = int(merged.ext_evidence.notna().sum())
    print(f'\nmatched onto candidates: {matched} / {len(merged)}')

    keys = set(zip(cand._p, cand._t))
    unmatched = [(p, t) for p, t in zip(ev._p, ev._t) if (p, t) not in keys]
    if unmatched:
        print(f'UNMATCHED evidence rows ({len(unmatched)}) - fix the spelling or add the pair:')
        for p, t in sorted(set(unmatched))[:20]:
            print(f'   {p} x {t}')

    merged = merged.drop(columns=['_p', '_t'])
    merged.to_csv(args.out, index=False, lineterminator='\n')
    print(f'\nwrote {args.out} ({len(merged)} rows, {merged.shape[1]} columns)')


if __name__ == '__main__':
    main()
