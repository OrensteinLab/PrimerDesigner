#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import argparse
import csv
import re
import subprocess
from collections import defaultdict

EFF_RE = re.compile(r'^([0-9]+(?:\.[0-9]+)?)\s+amplicon', re.IGNORECASE)

def read_fasta_one_sequence(path):
    """Read first FASTA record; return sequence as uppercase string (no whitespace)."""
    seq_parts = []
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if seq_parts:  # already read first record
                    break
                continue
            seq_parts.append(line)
    seq = ''.join(seq_parts).upper()
    if not seq:
        raise RuntimeError("No sequence found in FASTA: %s" % path)
    return seq

def norm_strand(s):
    s = (s or '').strip()
    if s in ('f', 'F', '+'):
        return 'F'
    if s in ('r', 'R', '-'):
        return 'R'
    return None

def run_testing_py(python2_exe, testing_py, template_seq, fwd, rev):
    """
    Run: python2 testing.py <template_seq> <fwd> <rev>
    Parse efficiency from output lines like: "<num> amplicon ..."
    """
    cmd = [python2_exe, testing_py, template_seq, fwd, rev]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()

    txt = (out or '') + "\n" + (err or '')
    for line in txt.splitlines():
        m = EFF_RE.match(line.strip())
        if m:
            try:
                return float(m.group(1))
            except Exception:
                return None
    return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--primers_csv", required=True,
                    help="CSV with columns: method,pair_id,strand,primer_seq_5to3")
    ap.add_argument("--template_fasta", required=True,
                    help="FASTA containing the full template sequence")
    ap.add_argument("--testing_py", default="testing.py",
                    help="Path to testing.py (python2 script)")
    ap.add_argument("--python2", default="python",
                    help="Python2 executable to run testing.py (default: python)")
    ap.add_argument("--verbose", action="store_true",
                    help="Print per-pair results")
    args = ap.parse_args()

    template_seq = read_fasta_one_sequence(args.template_fasta)

    # pairs[(method, pair_id)] = {'F': seq, 'R': seq}
    pairs = defaultdict(dict)

    with open(args.primers_csv, 'r') as f:
        r = csv.DictReader(f)
        required = set(['method', 'pair_id', 'strand', 'primer_seq_5to3'])
        fieldnames = set(r.fieldnames or [])
        missing = required - fieldnames
        if missing:
            raise RuntimeError("Missing columns in CSV: %s" % sorted(list(missing)))

        for row in r:
            method = (row.get('method') or '').strip()
            pair_id = (row.get('pair_id') or '').strip()
            strand = norm_strand(row.get('strand'))
            primer = (row.get('primer_seq_5to3') or '').strip().replace(' ', '').upper()
            if not method or not pair_id or not strand or not primer:
                continue
            # keep first occurrence if duplicates
            if strand not in pairs[(method, pair_id)]:
                pairs[(method, pair_id)][strand] = primer

    methods = sorted(set([m for (m, _) in pairs.keys()]))
    if not methods:
        raise RuntimeError("No (method, pair_id) entries parsed. Check your CSV content.")

    print("Template length: %d" % len(template_seq))
    print("Methods found: %s" % ", ".join(methods))

    for m in methods:
        keys = sorted([k for k in pairs.keys() if k[0] == m], key=lambda x: x[1])

        effs = []
        missing_fr = 0
        failed_parse = 0

        for (_, pid) in keys:
            d = pairs[(m, pid)]
            fwd = d.get('F')
            rev = d.get('R')
            if not fwd or not rev:
                missing_fr += 1
                if args.verbose:
                    print("[%s] pair_id=%s SKIP (missing F or R)" % (m, pid))
                continue

            eff = run_testing_py(args.python2, args.testing_py, template_seq, fwd, rev)
            if eff is None:
                failed_parse += 1
                if args.verbose:
                    print("[%s] pair_id=%s FAILED (parse)" % (m, pid))
            else:
                effs.append(eff)
                print("[%s] pair_id=%s eff=%.4f" % (m, pid, eff))

        ok = len(effs)
        total = len(keys)
        avg = (sum(effs) / float(ok)) if ok else float('nan')

        print("\n====================")
        print("Method: %s" % m)
        print("====================")
        print("Total pair_ids:     %d" % total)
        print("Successful runs:    %d" % ok)
        print("Missing F/R:        %d" % missing_fr)
        print("Failed parse:       %d" % failed_parse)
        if ok:
            print("Average efficiency: %.6f" % avg)
        else:
            print("Average efficiency: NA")

if __name__ == "__main__":
    main()