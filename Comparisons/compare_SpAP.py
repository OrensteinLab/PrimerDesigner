import re
import time
import tracemalloc
from pathlib import Path
import pandas as pd
from General.primer_graphs import create_primer_df, create_graph
import General.utils as GU
from General.args import get_args
import csv

# Match Olivar-style name: SPAP_<pair_num>_LEFT_1 or SPAP_<pair_num>_RIGHT_1
_BED_PAIR_NAME = re.compile(r"^\w+_(\d+)_(?:LEFT|RIGHT)_\d+$")


def read_bed_primers(path, cfg):
    """
    Returns (nodes, pair_ids).
    nodes: list of (start, end, strand) with strand 'f' or 'r'.
    pair_ids: list of same length; int (0-based pair id) when parseable from BED name column, else None.
    """
    primers = []
    pair_ids = []
    with open(path) as f:
        for lineno, line in enumerate(f, 1):
            line = line.strip()

            # skip blanks + header/comment lines
            if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue

            fields = line.split()  # split on ANY whitespace (tabs/spaces)

            if len(fields) < 6:
                raise ValueError(f"Malformed BED line {lineno}: {line}")

            start = int(fields[1]) - len(cfg.upstream)
            end   = int(fields[2]) - len(cfg.upstream)
            strand = fields[5]

            direction = "f" if strand == "+" else "r"
            primers.append((start, end, direction))

            # Parse pair id from name (column 4), 1-based in BED -> 0-based in CSV
            pid = None
            if len(fields) >= 4:
                m = _BED_PAIR_NAME.match(fields[3].strip())
                if m:
                    pid = int(m.group(1)) - 1
            pair_ids.append(pid)

    return primers, pair_ids

def export_null_paths_primers_py(
    paths,
    template_5to3,
    out_csv_path,
    protein_name="SPAP",
    method="NullWeighted",
    cfg = None
):
    """
    paths: list of paths, each path is ['s', (start,end,strand), ..., 'd']
    Writes a CSV of all primers from all null paths (Python 3 compatible).
    """

    with open(out_csv_path, "w", newline="") as f:
        writer = csv.writer(f)

        writer.writerow([
            "protein_name", "method", "null_path_id",
            "primer_order", "pair_id",
            "start", "end", "strand",
            "primer_seq_5to3", "length"
        ])

        for path_id, path in enumerate(paths):

            # strip source / sink
            core = path
            if core and core[0] == "s":
                core = core[1:]
            if core and core[-1] == "d":
                core = core[:-1]

            pair_id = -1

            for primer_order, node in enumerate(core):
                start, end, strand = node

                seq = primer_seq_from_template(
                    template_5to3, start, end, strand, cfg
                ).upper()

                if strand.lower() == "f":
                    pair_id += 1

                writer.writerow([
                    protein_name,
                    method,
                    path_id,
                    primer_order,
                    pair_id,
                    start,
                    end,
                    strand,
                    seq,
                    len(seq)
                ])

def primer_seq_from_template(template_5to3: str, start: int, end: int, strand: str, cfg) -> str:
    """
    start/end are 0-based, python slicing semantics [start:end)
    strand: 'f' or 'r'
    returns primer sequence in 5'->3' orientation as ordered.
    """
    subseq = template_5to3[start +len(cfg.upstream):end + len(cfg.upstream)]
    if strand.lower() == "f":
        return subseq
    elif strand.lower() == "r":
        return GU.revcomp(subseq)
    else:
        raise ValueError(f"Unknown strand: {strand}")



# -----------------------------
# Export function (same as Geller/PrimalScheme)
# -----------------------------
def export_primer_set(primer_df, nodes, protein_name, method, pair_ids=None):
    """
    pair_ids: optional list of int (0-based), same length as nodes.
    When provided and entry is not None, use it; otherwise increment pair_id on each forward primer (f).
    """
    rows = []
    pair_id = -1

    for order, node in enumerate(nodes):
        r = primer_df.loc[node]

        start, end, strand = node
        seq = r["seq"].upper()
        length = len(seq)

        if pair_ids is not None and order < len(pair_ids) and pair_ids[order] is not None:
            pair_id = pair_ids[order]
        elif strand == "f":
            pair_id += 1

        rows.append({
            "protein_name": protein_name,
            "method": method,
            "primer_order": order,
            "pair_id": pair_id,
            "start": start,
            "end": end,
            "strand": strand,
            "primer_seq_5to3": seq,
            "length": length
        })

    return pd.DataFrame(rows)


def main():
    
    # Load config
    cfg = GU.load_config("configs/SPAP_experiment.json")

    # ---- Input sequence ----
    mutreg_nt = GU.read_fasta("data/SPAP_reference.fa")

    # add upstream/downstream padding to the sequence for primer design
    sequence_nt = cfg.upstream + mutreg_nt + cfg.downstream

    protein_name = "SPAP"

    # get user arguments
    args = get_args()

    t0 = time.time()

    # ---- Build primer table and graph ----
    primer_df = create_primer_df(sequence_nt, args, cfg)

    t_graph0 = time.time()
    tracemalloc.start()
    graph = create_graph(primer_df, len(mutreg_nt), args)
    graph_time = time.time() - t_graph0
    graph_peak_mb = tracemalloc.get_traced_memory()[1] / 1e6
    tracemalloc.stop()

    longest_path_t0 = time.time()

    # find longest path in the graph (which corresponds to the primer set)
    full_path = GU.longest_path_dag(graph, "s", "d")
    primer_path_nodes = full_path[1:-1]

    # select the primers from the DataFrame that correspond to the longest path nodes
    primer_designer_set = primer_df.loc[primer_path_nodes].copy().reset_index()
    primer_designer_efficiency = float(primer_designer_set["efficiency"].mean())

    longest_path_time = time.time() - longest_path_t0
    total_time = time.time() - t0

    # ---- Load QuickChange primers (from paper) ----
    quick_primers, _ = read_bed_primers("Comparisons/Primers/QuickChange_primers.bed", cfg)

    # select only those primers from the dataframe that match the QuickChange primers
    quick_set = primer_df.loc[quick_primers].copy().reset_index()
    quick_efficiency = float(quick_set["efficiency"].mean())


    # ---- Load PrimalScheme primers (bed file) ----
    primal_scheme_primers, _ = read_bed_primers("Comparisons/Primers/PrimalScheme_primers.bed", cfg)

    # select primers from primer_df that match the PrimalScheme's primers
    primal_scheme_set = primer_df.loc[primal_scheme_primers].copy().reset_index()
    primal_scheme_efficiency = float(primal_scheme_set["efficiency"].mean())

    # extend primer design to include shorter primers
    args.primer_lmin = 15
    args.primer_lmax = 30
    extended_primer_set = create_primer_df(sequence_nt, args, cfg)

    olivar_primers, olivar_pair_ids = read_bed_primers("Comparisons/Primers/olivar-design.primer.bed", cfg)
    # select only those primers from the dataframe that match the QuickChange primers
    olivar_set = extended_primer_set.loc[olivar_primers].copy().reset_index()
    olivar_efficiency = float(olivar_set["efficiency"].mean())

    # ---- Create results directory ----
    results_dir = Path("Comparisons/Results")
    results_dir.mkdir(parents=True, exist_ok=True)

    # ---- Summary CSV ----
    row = {
        "protein_name": protein_name,
        "graph_nodes": len(graph.nodes),
        "graph_edges": len(graph.edges),
        "graph_time_sec": round(graph_time, 3),
        "graph_peak_mem_MB": round(graph_peak_mb, 1),
        "PD_single_avg_efficiency": primer_designer_efficiency,
        "QuickChange_avg_efficiency": quick_efficiency,
        "PrimalScheme_avg_efficiency": primal_scheme_efficiency,
        "PD_single_primers": len(primer_path_nodes),
        "QuickChange_primers": len(quick_primers),
        "PrimalScheme_primers": len(primal_scheme_primers),
        "olivar_primers": len(olivar_primers),
        "olivar_avg_efficiency": olivar_efficiency,
        "shortest_path_time": longest_path_time,
        "total_time_sec": round(total_time, 3),
    }

    pd.DataFrame([row]).to_csv(results_dir / "SpAP_comparison.csv", index=False)

    # ---- Primer lists CSV ----
    df_pd = export_primer_set(
        primer_df=primer_df,
        nodes=primer_path_nodes,
        protein_name=protein_name,
        method="PD_single",
    )

    df_olivar = export_primer_set(
        primer_df=extended_primer_set,
        nodes=olivar_primers,
        protein_name=protein_name,
        method="Olivar",
        pair_ids=olivar_pair_ids,
    )

    df_quick = export_primer_set(
        primer_df=primer_df,
        nodes=quick_primers,
        protein_name=protein_name,
        method="QuickChange",
    )

    df_primal = export_primer_set(
        primer_df=primer_df,
        nodes=primal_scheme_primers,
        protein_name=protein_name,
        method="PrimalScheme",
    )

    # concatenate all primer sets into one DataFrame and save to CSV
    pd.concat([df_pd, df_olivar, df_quick, df_primal], ignore_index=True).to_csv(
        results_dir / "SPAP_primers.csv",
        index=False,
    )

    # sample 1000 random paths from the graph for null distribution 
    sampled_paths = GU.sample_paths_dag_uniform(graph,  's', 'd', k=1000, max_tries=1000, seed=42)

    export_null_paths_primers_py(
        paths=sampled_paths,
        template_5to3=sequence_nt,
        out_csv_path=results_dir / "null_paths_primers_SPAP.csv",
        protein_name=protein_name,
        method="NullWeighted", cfg=cfg
    )
    print("✔ Null primer paths written")

if __name__ == "__main__":
    main()