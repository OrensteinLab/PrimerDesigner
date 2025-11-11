#!/usr/bin/env python3
from __future__ import annotations
import numpy as np
import multiprocessing as mp
from dataclasses import dataclass
from General.utils import * 


# ---- calc_tm wrapper (heterodimer Tm) ----
def calc_tm(seq1, seq2):
    return PCR.calc_heterodimer(
        seq1, seq2,).tm

# ---------------------------
# Worker globals (per-process)
# ---------------------------
_G = {}

def _init_worker_stage1(seq1, seq2, lmin, lmax, max_tm, allowed_overlap):
    _G["seq1"] = seq1
    _G["seq2"] = seq2
    _G["lmin"] = lmin
    _G["lmax"] = lmax
    _G["max_tm"] = max_tm
    _G["allowed_overlap"] = allowed_overlap

def _init_worker_stage2(seq1, seq2, lmin, lmax, max_tm, upstream_len):
    _G["seq1"] = seq1
    _G["seq2"] = seq2
    _G["lmin"] = lmin
    _G["lmax"] = lmax
    _G["max_tm"] = max_tm
    _G["up"] = upstream_len

# ---------------------------
# Stage 1: candidates at Lmax
# ---------------------------
def _stage1_chunk_inter(i_range):
    """Find Lmax×Lmax candidate (p1,p2) across two sequences."""
    seq1 = _G["seq1"]; seq2 = _G["seq2"]
    lmax = _G["lmax"]; max_tm = _G["max_tm"]
    n1, n2 = len(seq1), len(seq2)
    i_start, i_end = i_range

    buf = []
    for p1 in range(i_start, i_end):
        if p1 + lmax > n1:
            break
        sub1 = seq1[p1:p1 + lmax]
        for p2 in range(0, n2 - lmax + 1):
            sub2 = seq2[p2:p2 + lmax]
            # NOTE: you must define calc_tm() globally
            if calc_tm(sub1, sub2) >= max_tm:
                buf.extend((p1, p1 + lmax, p2, p2 + lmax))
    if not buf:
        return np.empty((0, 4), dtype=np.int32)
    return np.fromiter(buf, dtype=np.int32).reshape(-1, 4)

def _stage1_chunk_intra(i_range):
    """Find Lmax×Lmax candidate (p1,p2) within the same sequence (seq1==seq2)."""
    seq = _G["seq1"]              # same object as _G["seq2"]
    lmax = _G["lmax"]
    max_tm = _G["max_tm"]
    sigma = _G["allowed_overlap"]
    L = len(seq)
    i_start, i_end = i_range

    buf = []
    for p1_start in range(i_start, i_end):
        p1_end = p1_start + lmax
        if p1_end > L:
            break
        sub1 = seq[p1_start:p1_end]

        # Lower bound = max(i+1, i+lmax−sigma), so:
        #  - only scan to the right (avoid duplicates)
        #  - window overlap ≤ sigma
        j_lo = max(p1_start + 1, p1_end - sigma)
        j_hi = L - lmax + 1

        for p2_start in range(j_lo, j_hi):
            sub2 = seq[p2_start:p2_start + lmax]
            if calc_tm(sub1, sub2)  >= max_tm:
                buf.extend((p1_start, p1_end, p2_start, p2_start + lmax))

    if not buf:
        return np.empty((0, 4), dtype=np.int32)
    return np.fromiter(buf, dtype=np.int32).reshape(-1, 4)


# ---------------------------
# Stage 2: expand candidates
# ---------------------------
def _stage2_expand(cand_row):
    """Expand one Lmax×Lmax candidate into all sub-length cross-hyb hits."""
    seq1 = _G["seq1"]; seq2 = _G["seq2"]
    lmin = _G["lmin"]; lmax = _G["lmax"]; max_tm = _G["max_tm"]; up = _G["up"]
    p1_start, p1_end, p2_start, p2_end = map(int, cand_row)
    primer1 = seq1[p1_start:p1_end]
    primer2 = seq2[p2_start:p2_end]

    out_buf = []
    for len1 in range(lmin, lmax + 1):
        max_s1 = lmax - len1
        for len2 in range(lmin, lmax + 1):
            max_s2 = lmax - len2
            for s1 in range(max_s1 + 1):
                sub1 = primer1[s1:s1 + len1]
                for s2 in range(max_s2 + 1):
                    sub2 = primer2[s2:s2 + len2]
                    if calc_tm(sub1, sub2)  >= max_tm:
                        a = p1_start + s1 - up
                        b = a + len1
                        c = p2_start + s2 - up
                        d = c + len2
                        out_buf.extend((a, b, c, d))
    if not out_buf:
        return np.empty((0, 4), dtype=np.int32)
    return np.fromiter(out_buf, dtype=np.int32).reshape(-1, 4)

# ---------------------------
# Generic two-stage driver
# ---------------------------
def _two_stage_driver(
    seq1: str,
    seq2: str,
    args,
    stage1_func,
) -> set[tuple[tuple[int,int], tuple[int,int]]]:
    """Two-stage pipeline with parallel Stage 1 (candidates) and Stage 2 (expansion)."""

    lmin, lmax, max_tm = args.primer_lmin, args.primer_lmax, MAX_TM
    n1 = len(seq1)
    last = n1 - lmax + 1  # p1 starts that fit Lmax
    if last <= 0:
        return set()

    # find number of available processors
    procs = max(1, mp.cpu_count() - 1)

    # stage-1 chunk (how many p1 positions per task)
    chunk = last // max(1, procs)

    ranges = [(i, min(i + chunk, last)) for i in range(0, last, chunk)]
    if not ranges:
        return set()

    ctx = mp.get_context("spawn")

    # ----- Stage 1: Lmax candidates -----
    with ctx.Pool(
        processes=procs,
        initializer=_init_worker_stage1,
        initargs=(seq1, seq2, lmin, lmax, max_tm, args.allowed_overlap),
    ) as pool:
        cand_iter = pool.imap_unordered(stage1_func, ranges, chunksize=1)
        cand_list = [arr for arr in cand_iter if arr.size]

    if not cand_list:
        return set()

    candidates = np.vstack(cand_list)

    # ----- Stage 2: expand survivors -----
    with ctx.Pool(
        processes=procs,
        initializer=_init_worker_stage2,
        initargs=(seq1, seq2, lmin, lmax, max_tm, len(upstream_nt)),
    ) as pool:
        parts_iter = pool.imap_unordered(_stage2_expand, candidates, chunksize=1)
        out: set[tuple[tuple[int,int], tuple[int,int]]] = set()
        for arr in parts_iter:
            if arr.size:
                out.update(((int(a), int(b)), (int(c), int(d))) for a, b, c, d in arr)
    return out

# ---------------------------
# Public API (two functions)
# ---------------------------
def find_forbidden_pairs_inter(seq1: str, seq2: str, args):
    """
    Find all cross-hybridizing primer pairs BETWEEN two sequences (inter-sequence).

    """
    return _two_stage_driver(seq1, seq2, args, stage1_func=_stage1_chunk_inter)

def find_forbidden_pairs_intra(seq: str, args):
    """
    Find all cross-hybridizing primer pairs WITHIN the same sequence (intra-sequence).
  
    """
    return _two_stage_driver(seq, seq, args, stage1_func=_stage1_chunk_intra)

