import os
import time
import json
import tracemalloc
import pandas as pd
from pathlib import Path
from General.primer_graphs import *
from PD_var_ILP.greedy import *
from PD_var_ILP.ilp_model import *
from General.args import *
import sys


sys.argv = [
    sys.argv[0],
    "--file_path", "data/10_protein_coding_sequences.txt",
    "--output", "Experiment_results/var_ILP_different_proteins",
]

args = get_args()

# Create output directory if not exists
output_dir = Path(args.output)
output_dir.mkdir(parents=True, exist_ok=True)
print(f"[INFO] Output directory: {output_dir.resolve()}")

# ============================================================
# LOAD SEQUENCES
# ============================================================
print(f"[INFO] Reading protein coding sequences from: {args.file_path}")
all_mutreg_regions, all_full_sequences, all_protein_names = read_sequences(args.file_path)
print(f"[INFO] Total proteins loaded: {len(all_protein_names)}")

# ============================================================
# MAIN LOOP
# ============================================================
overall_start = time.time()

for i,protein in enumerate(all_protein_names):

    print(f"\n[INFO] Running for protein: {protein}. Protein length: {len(all_mutreg_regions[i])}")

    mutreg_nt = all_mutreg_regions[i]
    sequence_nt = UPSTREAM_NT + mutreg_nt + DOWNSTREAM_NT

    # ---- Step 1: Build primers and graph ----
    print("[STEP 1] Creating primer DataFrame...")
    primer_df = create_primer_df(sequence_nt, args)

    print(primer_df.head(5))

    print("[STEP 2] Creating primer graph...")
    t0 = time.time()
    tracemalloc.start()
    graph = create_graph(primer_df, len(mutreg_nt), args)
    graph_time = time.time() - t0
    graph_peak_mb = tracemalloc.get_traced_memory()[1] / 1e6
    tracemalloc.stop()

    # ---- Step 2: Greedy solution ----
    print("[STEP 3] Running greedy algorithm...")
    t1 = time.time()
    greedy_solution, greedy_obj = run_greedy(graph, primer_df, args)
    greedy_time = time.time() - t1

    # ---- Step 3: ILP model ----
    print("[STEP 4] Solving ILP model...")
    ilp_res: ILPResult = ilp_model(graph, sequence_nt, mutreg_nt, args)

    # ---- Step 4: Save CSV summary ----
    results = {
    "protein_name": protein,
    "seq_length": len(sequence_nt),
    "graph_nodes": len(graph.nodes),
    "graph_edges": len(graph.edges),
    "graph_time_sec": round(graph_time, 3),
    "graph_peak_mem_MB": round(graph_peak_mb, 1),
    

    # ---- ILP (updated field names) ----
    "ilp_num_vars": ilp_res.num_vars,
    "ilp_num_constraints": ilp_res.num_constraints,
    "ilp_setup_time_sec": round(ilp_res.setup_time, 3),
    "ilp_setup_peak_mem_MB": round(ilp_res.setup_peak_mem_mb, 1),
    "ilp_optimize_time_sec": round(ilp_res.optimize_time, 3),
    "ilp_optimize_peak_mem_MB": round(ilp_res.optimize_peak_mem_mb, 1),
    "ilp_objective": ilp_res.objective,
    "ilp_status": ilp_res.status,
    "ilp_path_length": sum(len(path) for path in ilp_res.paths),
    

    # ---- Greedy ----
    "greedy_objective": greedy_obj,
    "greedy_time_sec": round(greedy_time, 3),
    "greedy_path_length": sum(len(path) for path in greedy_solution),
    }

    csv_path = output_dir / f"results_protein_{protein}.csv"
    pd.DataFrame([results]).to_csv(csv_path, index=False)
    print(f"Saved CSV summary: {csv_path}")

    # ---- Step 5: Save JSON paths ----
    paths_out = {
        "protein_name": protein,
        "seq_length": len(sequence_nt),
        "ilp_path": ilp_res.paths,
        "greedy_path": greedy_solution,
    }

    json_path = output_dir / f"paths_protein_{protein}.json"
    with open(json_path, "w") as f:
        json.dump(paths_out, f, indent=2)
    print(f"Saved JSON paths: {json_path}")


overall_end = time.time()
total_minutes = (overall_end - overall_start) / 60
print(f"\n[ALL DONE] All protein sequences processed successfully.")
print(f"[TOTAL RUNTIME] {total_minutes:.2f} minutes ({overall_end - overall_start:.1f} seconds total)")