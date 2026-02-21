import time
import tracemalloc
from pathlib import Path
import pandas as pd
import networkx as nx
import General.utils as GU
from General.primer_graphs import create_primer_df, create_graph
from General.args import get_args
import csv


def read_primers_csv(path):
    primers = {}
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            tile = int(row["tile"])
            primers[tile] = (row["fwd"].upper(), row["rev"].upper())
    return primers


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


def export_null_paths_primers_py(
    paths,
    template_5to3,
    out_csv_path,
    protein_name="CVB3",
    method="NullWeighted",
    cfg=None,
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
                    template_5to3, start, end, strand, cfg,
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


def export_primer_set(primer_df, nodes, protein_name, method):
    rows = []
    pair_id = -1

    for order, node in enumerate(nodes):
        r = primer_df.loc[node]

        start, end, strand = node
        seq = r["seq"].upper()
        length = len(seq)

        # Increment pair_id on forward primer
        if strand == "f":
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


def rc(seq: str) -> str:
    tbl = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(tbl)[::-1]


def find_unique(haystack: str, needle: str, prev_pos):
    """Return 0-based index of first occurrence at/after prev_pos; raise if not found."""
    p = haystack.find(needle, prev_pos)
    if p == -1:
        raise SystemExit(f"Primer {needle} not found!")
    return p



def main():

    CVB3_PRIMERS = read_primers_csv("Comparisons/CVB3_primers.csv")

    # load upstream/downstream context and other parameters from config
    cfg = GU.load_config("configs/CVB3_experiment.json")

    mutreg_nt = GU.read_fasta("data/CVB3_reference.fa")

    # append upstream/downstream padding to the sequence for primer design
    sequence_nt = cfg.upstream + mutreg_nt + cfg.downstream

    sequence_nt = sequence_nt.upper()

    # get user arguments (e.g. for primer design parameters like Tm, GC content, etc.)
    args = get_args()

    # adjust oligo lengths to match CVB3 primers of length ~250bp
    args.oligo_lmin = 240
    args.oligo_lmax = 260

    t0 = time.time()

    # ---- Build primer table ----
    primer_df = create_primer_df(sequence_nt, args, cfg)

    competing_primers = []
    search_pos = 0  # pointer for searching forward primers

    for tile in sorted(CVB3_PRIMERS.keys()):
        fwd = CVB3_PRIMERS[tile][0].upper()
        rev = CVB3_PRIMERS[tile][1].upper()
        rc_rev = rc(rev)

        # find forward primer starting at search_pos
        f_hit = find_unique(sequence_nt, fwd, search_pos)

        # find reverse primer AFTER the forward primer
        r_hit = find_unique(sequence_nt, rc_rev, f_hit + len(fwd))

        # append correct positions (convert back to coords relative to mutreg start)
        competing_primers.append(
            (f_hit - len(cfg.upstream), f_hit + len(fwd) - len(cfg.upstream), "f")
        )
        competing_primers.append(
            (r_hit - len(cfg.upstream), r_hit + len(rev) - len(cfg.upstream), "r")
        )

        # advance search position for next forward primer to be after this one
        search_pos = f_hit + len(fwd)

    print("All primers were found in sequence!")

    # choose only those primers from the graph that match the CVB3 primers
    CVB3_set = primer_df.loc[competing_primers].copy().reset_index()
    CVB3_efficiency = float(CVB3_set["efficiency"].mean())

    # ---- Build graph (time + peak mem) ----
    t_graph0 = time.time()
    tracemalloc.start()
    graph = create_graph(primer_df, len(mutreg_nt), args)
    graph_time = time.time() - t_graph0
    graph_peak_mb = tracemalloc.get_traced_memory()[1] / 1e6
    tracemalloc.stop()

    # ---- find longest path ----
    longest_path_t0 = time.time()

    full_path = GU.longest_path_dag(graph, "s", "d")

    primer_path_nodes = full_path[1:-1]  # strip s/d vertices

    # check that all nodes in path are in primer_df
    primer_set = primer_df.loc[primer_path_nodes].copy().reset_index()
    primer_efficiency = float(primer_set["efficiency"].mean())

    longest_path_time = time.time() - longest_path_t0
    total_time = time.time() - t0

    # ---- Create results directory ----
    results_dir = Path("Comparisons/Results")
    results_dir.mkdir(parents=True, exist_ok=True)

    # ---- Summary CSV ----
    row = {
        "graph_nodes": len(graph.nodes),
        "graph_edges": len(graph.edges),
        "graph_time_sec": round(graph_time, 3),
        "graph_peak_mem_MB": round(graph_peak_mb, 1),
        "PD_avg_efficiency": primer_efficiency,
        "CVB3_avg_efficiency": CVB3_efficiency,
        "PD_single_primers": len(primer_path_nodes),
        "CVB3_primers": len(competing_primers),
        "longest_path_time": longest_path_time,
        "total_time_sec": round(total_time, 3),
    }

    df = pd.DataFrame([row])
    summary_path = results_dir / "CVB3_comparison.csv"
    df.to_csv(summary_path, index=False)

    # ---- Primer lists CSV ----
    pd_nodes = primer_path_nodes
    CVB3_nodes = competing_primers
    protein_name = "CVB3"

    df_pd = export_primer_set(
        primer_df=primer_df,
        nodes=pd_nodes,
        protein_name=protein_name,
        method="PD_single",
    )

    df_CVB3 = export_primer_set(
        primer_df=primer_df,
        nodes=CVB3_nodes,
        protein_name=protein_name,
        method="CVB3",
    )

    df_all = pd.concat([df_pd, df_CVB3], ignore_index=True)

    primer_path = results_dir / "CVB3_primers.csv"
    df_all.to_csv(primer_path, index=False)

    # sample 1000 random paths for null distribution 
    sampled_paths = GU.sample_paths_dag_uniform(graph,  's', 'd', k=1000, max_tries=1000, seed=42)

    # save null paths
    export_null_paths_primers_py(
        paths=sampled_paths,
        template_5to3=sequence_nt,
        out_csv_path=results_dir / "null_paths_primers_CVB3.csv",
        protein_name=protein_name,
        method="NullWeighted", cfg=cfg
    )
    print("✔ Null primer paths written")


if __name__ == "__main__":
    main()