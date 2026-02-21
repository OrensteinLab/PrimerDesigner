import re

def parse_tile_id(name: str) -> int:
    # expects something like: 0b388d08_5_LEFT_1
    m = re.search(r"_(\d+)_", name)
    if not m:
        raise ValueError(f"Can't parse tile id from name: {name}")
    return int(m.group(1))

def load_tiles_from_bed(bed_path: str):
    tiles = {}  # tile_id -> dict(left=(start,end), right=(start,end))
    with open(bed_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 6:
                continue  # skip weird lines

            chrom, start, end, name, score, strand = parts[:6]
            start, end = int(start), int(end)
            tile = parse_tile_id(name)

            tiles.setdefault(tile, {"left": None, "right": None})

            if strand == "+":
                tiles[tile]["left"] = (start, end)
            elif strand == "-":
                tiles[tile]["right"] = (start, end)

    return tiles

def main():
    bed_path = "Comparisons/PrimalScheme_primers.bed"  # change path if needed
    tiles = load_tiles_from_bed(bed_path)

    # build oligo spans
    oligos = []  # (tile, start, end)
    for tile in sorted(tiles.keys()):
        left = tiles[tile]["left"]
        right = tiles[tile]["right"]
        if left is None or right is None:
            print(f"[WARN] tile {tile}: missing left or right primer (left={left}, right={right})")
            continue

        oligo_start = left[0]
        oligo_end = right[1]   # span until right primer end
        if oligo_end <= oligo_start:
            print(f"[WARN] tile {tile}: weird span [{oligo_start},{oligo_end})")
            continue

        oligos.append((tile, oligo_start, oligo_end))

    oligos.sort(key=lambda x: x[1])

    # check overlap >= 1 between consecutive oligos
    ok = True
    for (t1, s1, e1), (t2, s2, e2) in zip(oligos, oligos[1:]):
        overlap = min(e1, e2) - max(s1, s2)
        if overlap < 1:
            ok = False
            gap = max(s1, s2) - min(e1, e2)
            print(f"[FAIL] tile {t1} [{s1},{e1}) vs tile {t2} [{s2},{e2}): "
                  f"overlap={overlap} (gap={gap})")

    if ok:
        print("[OK] Every consecutive oligo overlaps the next by at least 1 bp.")
    else:
        print("[DONE] Found at least one pair that does not overlap by >= 1 bp.")

if __name__ == "__main__":
    main()