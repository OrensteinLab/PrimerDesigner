import time
import json
import tracemalloc
from pathlib import Path
import pandas as pd
import networkx as nx
from General.primer_graphs import create_primer_df, create_graph
from General.utils import *

def run_shortest_path(sequence_nt, mutreg_nt, protein_name, args):
    """Single-protein shortest-path run:
    - Build graph
    - Compute shortest path (s->d)
    - Save only metadata to CSV
    - Save the actual path to JSON
    """
    t0 = time.time()

    # ---- Build primer table and graph ----
    primer_df = create_primer_df(sequence_nt, args)

    t_graph0 = time.time()
    tracemalloc.start()
    graph = create_graph(primer_df, len(mutreg_nt), args)
    graph_time = time.time() - t_graph0
    graph_peak_mb = tracemalloc.get_traced_memory()[1] / 1e6
    tracemalloc.stop()

    # ---- Shortest path (by edge weight) ----
    # NetworkX returns list of nodes: ['s', <primer_nodes...>, 'd']
    try:
        full_path = longest_path_dag(graph, 's', 'd')
        # strip s/d for primer nodes only (they should match primer_df index)
        primer_path_nodes = full_path[1:-1]
        # cost = sum of node weights as you stored them on edges entering those nodes
        # but you already stored the per-primer cost in primer_df, so:
        primer_set = primer_df.loc[primer_path_nodes].copy().reset_index()
        primer_cost = float(primer_set['efficiency'].sum())
        status = "OK"
    except nx.NetworkXNoPath:
        primer_path_nodes = []
        primer_cost = float('nan')
        status = "NO_PATH"

    total_time = time.time() - t0

    # ---- CSV (no path in CSV) ----
    row = {
        "protein_name": protein_name,
        "graph_nodes": len(graph.nodes),
        "graph_edges": len(graph.edges),
        "graph_time_sec": round(graph_time, 3),
        "graph_peak_mem_MB": round(graph_peak_mb, 1),
        "shortest_path_cost": primer_cost,
        "status": status,
        "total_time_sec": round(total_time, 3),
    }
    df = pd.DataFrame([row])

    # Ensure output directory exists
    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ---- CSV summary ----
    csv_path = out_dir / "PD_single_results.csv"
    df.to_csv(csv_path, index=False)

    # ---- Paths JSON ----
    paths_json_path = out_dir / "primer_paths.json"
    with open(paths_json_path, "w") as f:
        json.dump({
            "protein_name": protein_name,
            "shortest_path_nodes": primer_path_nodes
        }, f, indent=2)

    print(f"✅ Saved CSV summary to: {csv_path}")
    print(f"✅ Saved path details to: {paths_json_path}")

