#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import argparse
import csv
import re
import subprocess
from collections import defaultdict

EFF_RE = re.compile(r'^([0-9]+(?:\.[0-9]+)?)\s+amplicon', re.IGNORECASE)

def read_fasta(path):
    seq = []
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('>'):
                continue
            seq.append(line)
    s = ''.join(seq).upper()
    if not s:
        raise RuntimeError("Empty FASTA: %s" % path)
    return s

def norm_strand(s):
    s = ('' if s is None else str(s)).strip()
    if s in ('f', 'F', '+'):
        return 'F'
    if s in ('r', 'R', '-'):
        return 'R'
    return None

def run_py2_eff(py2, testing_py, template_seq, fwd, rev):
    """
    Python2-compatible subprocess call.
    Runs: py2 testing_py template_seq fwd rev
    Returns efficiency float or None if parse failed.
    """
    cmd = [py2, testing_py, template_seq, fwd, rev]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()

    txt = (out or '') + '\n' + (err or '')
    for line in txt.splitlines():
        m = EFF_RE.match(line.strip())
        if m:
            try:
                return float(m.group(1))
            except Exception:
                return None
    return None  # parse failed

def sort_key_pid(x):
    # mimic: key=lambda x: int(x) if x.isdigit() else x
    try:
        return (0, int(x))
    except Exception:
        return (1, x)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--null_csv", required=True)
    ap.add_argument("--template_fasta", required=True)
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--python2", default="python")
    ap.add_argument("--path_col", default="null_path_id")
    args = ap.parse_args()

    testing_py = 'testing.py'  # assumed to be in same dir as this script

    template = read_fasta(args.template_fasta)

    # paths[path_id][pair_id]['F'/'R'] = primer
    paths = defaultdict(lambda: defaultdict(dict))

    with open(args.null_csv, 'r') as f:
        r = csv.DictReader(f)
        for row in r:
            pid = str(row.get(args.path_col, '')).strip()
            pair = str(row.get('pair_id', '')).strip()
            strand = norm_strand(row.get('strand', ''))
            primer = str(row.get('primer_seq_5to3', '')).strip().upper()
            if (not pid) or (not pair) or (not strand) or (not primer):
                continue
            # keep first
            if strand not in paths[pid][pair]:
                paths[pid][pair][strand] = primer

    out_rows = []

    for pid in sorted(paths.keys(), key=sort_key_pid):
        effs = []
        missing = 0
        failed = 0

        # NOTE: dict iteration order is not guaranteed in py2.
        # If you want deterministic order, sort pair_ids:
        pair_ids = sorted(paths[pid].keys(), key=sort_key_pid)

        for pair_id in pair_ids:
            d = paths[pid][pair_id]
            fwd = d.get('F')
            rev = d.get('R')
            if (not fwd) or (not rev):
                missing += 1
                continue

            eff = run_py2_eff(args.python2, testing_py, template, fwd, rev)
            if eff is None:
                failed += 1
            else:
                effs.append(eff)

        avg = (sum(effs) / float(len(effs))) if effs else float('nan')

        out_rows.append({
            'path_id': pid,
            'n_pairs': len(paths[pid]),
            'n_success': len(effs),
            'n_missing_F_or_R': missing,
            'n_failed_parse': failed,
            'avg_efficiency': avg,
        })

        print("path %s: avg_eff=%s (ok=%d, missing=%d, failed=%d)" %
              (pid, ("%.6f" % avg) if effs else "NA", len(effs), missing, failed))

    if not out_rows:
        raise RuntimeError("No paths parsed from %s (check columns / content)." % args.null_csv)

    with open(args.out_csv, 'w') as f:
        fieldnames = list(out_rows[0].keys())
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(out_rows)

if __name__ == "__main__":
    main()