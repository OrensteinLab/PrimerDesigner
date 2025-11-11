#!/usr/bin/env python3
"""
Assemble P2/P3 tiles using primer-anchored joins.

For each pair (Tile i -> Tile i+1):
- Find where RC(Rev_Tile_i) occurs in the **LAST** oligo sequence of Tile i+1.
- Append from that position in the **FIRST** oligo sequence of Tile i+1.
- Keep primers (no trimming before the join pos). Raise if anchor not found.

Inputs:
- Folder "Geller_Study" with Tile_1_oligos.fa ... Tile_18_oligos.fa
- PRIMERS dict below.

Output:
- assembled_primer_anchored_first_seq_last_pos.fasta
"""

import os, re, glob
from pathlib import Path

TILES_DIR = "Geller_Study"
OUTPUT_FASTA = "assembled_primer_anchored_first_seq_last_pos.fasta"

PRIMERS = {
    1:  ("taccactactaggcaaagcatcact",           "ctcttggacctctactagacctgg"),
    2:  ("gcactacccaatttcgtttgaagga",           "ccgaatgcatttccaagctgttc"),
    3:  ("gaacagggagtgaaggactatgtg",            "cagccatagggattccgtaatattg"),
    4:  ("gtggctcaaacagaaggtgtca",              "gttcctggtcactttgggatgg"),
    5:  ("cacaatcgagcagagcgcg",                 "tttctcagcaagcgaccttcc"),
    6:  ("caagtcggtggcaacaaacttaatt",           "cggtgaggtgaacagaatgcc"),
    7:  ("catggctgccctagaagagaaa",              "aattgaccgggcaacactcatc"),
    8:  ("cccatgtcagtcaagacttgtgac",            "gtgcaacgctaattttgatctctct"),
    9:  ("caccgagatgtttagggagtacaat",           "ccactgacacaaatgtggtcaa"),
    10: ("ggctttcatttgcttacaggca",              "gcccacctgtcatagatgcc"),
    11: ("gaatatggcgagtttaccatgctg",            "gttgggaaacttgctggtgttaat"),
    12: ("ggttaatgaggcagtgctagca",              "caccttgctcatcattgaagtagtg"),
    13: ("cttctcagcagcactcctcaaa",              "gcgtagtggtccactgcttc"),
    14: ("acacgtggatgagtacatgctg",              "cttctctatggacctgagctcat"),
    15: ("cctgaacctaccaatggtgacttat",           "ccagacagggcttaagctagc"),
    16: ("agcatttgattactctgggtacgat",           "caccatatgcgatcatcctgaattg"),
    17: ("gtgtacaaagggattgacttggac",            "ttgattcgtgtatgtctttcatggg"),
    18: ("cttcctggtgcatcctgttatg",              "agcacagtagggttaagccaa"),
}

def rc(s: str) -> str:
    return s.translate(str.maketrans("ACGTacgt","TGCAtgca"))[::-1]

def extract_tile_index(path: str) -> int:
    m = re.search(r"Tile[_-]?(\d+)", os.path.basename(path))
    return int(m.group(1)) if m else 10**9

def read_all_fasta(path: str):
    recs, h, buf = [], None, []
    with open(path) as f:
        for line in f:
            line=line.strip()
            if not line: continue
            if line.startswith(">"):
                if h is not None:
                    recs.append((h, "".join(buf).upper()))
                h=line[1:].strip(); buf=[]
            else:
                if h is not None: buf.append(line)
        if h is not None:
            recs.append((h, "".join(buf).upper()))
    if not recs:
        raise ValueError(f"No FASTA records in {path}")
    return recs

# load tiles, keep both FIRST and LAST sequences
tile_paths = sorted(glob.glob(str(Path(TILES_DIR)/"Tile_*_oligos.fa")), key=extract_tile_index)
if not tile_paths:
    raise SystemExit(f"No files matched {TILES_DIR}/Tile_*_oligos.fa")

tiles = []  # (idx, file_name, first_header, first_seq, last_header, last_seq)
for p in tile_paths:
    idx = extract_tile_index(p)
    if idx not in PRIMERS:
        raise SystemExit(f"No primers for Tile_{idx}")
    recs = read_all_fasta(p)
    first_h, first_s = recs[0]
    last_h,  last_s  = recs[-1]
    tiles.append((idx, os.path.basename(p), first_h, first_s, last_h, last_s))

# assemble:
assembled = tiles[0][3]  # start from FIRST oligo of Tile 1
joins = []

for i in range(len(tiles)-1):
    idx_i, name_i, _, first_i, _, last_i = tiles[i]
    idx_j, name_j, first_h_j, first_s_j, last_h_j, last_s_j = tiles[i+1]

    anchor = rc(PRIMERS[idx_i][1]).upper()  # RC(Rev_Tile_i)
    pos_in_last = last_s_j.find(anchor)
    if pos_in_last == -1:
        raise RuntimeError(
            "Primer-anchored join failed:\n"
            f"  prev: {name_i} (Tile {idx_i})\n"
            f"  next: {name_j} (Tile {idx_j})\n"
            f"  anchor: {anchor}\n"
            f"  not found in LAST oligo header: {last_h_j}\n"
            f"  next LAST prefix: {last_s_j[:60]}..."
        )

    # append from the SAME position in the FIRST oligo’s sequence
    if pos_in_last > len(first_s_j):
        raise RuntimeError(
            f"Join position beyond FIRST sequence length for {name_j}: "
            f"pos={pos_in_last}, first_len={len(first_s_j)}"
        )

    joins.append((name_i, name_j, pos_in_last))
    assembled += first_s_j[pos_in_last:]

with open(OUTPUT_FASTA, "w") as out:
    out.write(">assembled_primer_anchored_first_seq_last_pos\n")
    for i in range(0, len(assembled), 80):
        out.write(assembled[i:i+80] + "\n")

print(f"Wrote {OUTPUT_FASTA}")
print("Joins (file_i, file_j, position_in_last_of_j):")
for j in joins:
    print("  ", j)
