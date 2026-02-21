import time
from pathlib import Path
import pandas as pd
import networkx as nx
from General.primer_graphs import create_primer_df, create_graph
from General.utils import *

def run_longest_path(sequence_nt, mutreg_nt, protein_name, args, cfg):
    """Single-protein shortest-path run:
    - Build graph
    - Compute shortest path (s->d)
    - Save only metadata to CSV
    - Save the actual path to JSON
    """
    t0 = time.time()

    # ---- Build primer table and graph ----
    primer_df = create_primer_df(sequence_nt, args, cfg)

    t_graph0 = time.time()
    graph = create_graph(primer_df, len(mutreg_nt), args)
    graph_time = time.time() - t_graph0

    # ---- Longest path (by edge weight) ----
    try:
        full_path = longest_path_dag(graph, 's', 'd')
        # strip s/d for primer nodes only)
        primer_path_nodes = full_path[1:-1]
        # select primers in path and compute total efficiency
        primer_set = primer_df.loc[primer_path_nodes].copy().reset_index()
        primer_cost = float(primer_set['efficiency'].sum())

    except nx.NetworkXNoPath:
        print("[WARN] No valid path found for {protein_name}. Skipping.")
        primer_path_nodes = []
        primer_cost = float('nan')

    total_time = time.time() - t0

    # ---- CSV data ----
    row = {
        "protein_name": protein_name,
        "graph_nodes": len(graph.nodes),
        "graph_edges": len(graph.edges),
        "graph_time_sec": round(graph_time, 3),
        "longest_path_efficiency": primer_cost,
        "total_time_sec": round(total_time, 3),
        "num_primers": len(primer_path_nodes),
    }
    df = pd.DataFrame([row])

    # Ensure output directory exists
    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ---- CSV summary ----
    csv_path = out_dir / "PD_single_LPath_results.csv"
    df.to_csv(csv_path, index=False)

    # ---- Paths CSV ----
    paths_csv_path = out_dir / "PD_single_LPath_selected_primers.csv"

    # save selected primers in csv
    if len(primer_path_nodes) > 0:
        primer_set = primer_set.rename(columns={"fr": "orientation"})
        primer_set.insert(0, "protein_name", protein_name)
        primer_set.insert(1, "primer_index", range(len(primer_set)))
        primer_set = primer_set[['protein_name','primer_index','start','stop','orientation','seq','efficiency']]
        primer_set.to_csv(paths_csv_path, index=False)
    else:
        pd.DataFrame().to_csv(paths_csv_path, index=False)

    print(f"Saved selected primers to: {paths_csv_path}")

    return df, primer_path_nodes

