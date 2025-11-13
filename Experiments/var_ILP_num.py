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

# ============================================================
# CONFIGURATION
# ============================================================
sys.argv = [
    sys.argv[0],
    "--file_path", "None",
    "--output", "Experiment_results/var_ILP_num_variants",
]

args = get_args()

# Create output directory if not exists
output_dir = Path(args.output)
output_dir.mkdir(parents=True, exist_ok=True)
print(f"[INFO] Output directory: {output_dir.resolve()}")

# ============================================================
# CONSTANT SEQUENCE AND PARAMETERS
# ============================================================
mutreg_nt_full = (
    "ATGCAAAGCCCAGCACCTGCCGCAGCGCCTGCCCCTGCGGCACGTTCCATCGCAGCTACGCCTCCTAAACTGATCGTGGCAATTAGCGTGGACCAGTTTAGTGCAGACTTGTTCTCGGAGTATCGTCAATATTACACCGGAGGTTTAAAGCGTCTTACATCCGAAGGAGCTGTGTTCCCACGTGGTTATCAGAGTCATGCGGCAACAGAAACGTGTCCTGGTCACTCAACGATCCTGACAGGATCACGTCCGTCACGTACGGGTATTATCGCTAATAACTGGTTCGACTTGGACGCAAAGCGTGAGGATAAAAATCTGTACTGTGCTGAGGATGAATCCCAACCCGGTAGTTCGTCTGACAAGTACGAAGCTTCGCCACTGCACTTAAAGGTACCCACCCTGGGGGGACGCATGAAAGCCGCCAATCCTGCGACTCGTGTCGTCTCTGTTGCCGGCAAGGATCGCGCGGCCATTATGATGGGTGGCGCCACAGCGGATCAGGTCTGGTGGTTAGGGGGGCCTCAGGGGTATGTTTCGTATAAGGGTGTAGCGCCAACTCCCCTTGTAACACAGGTCAATCAGGCCTTTGCACAGCGCTTAGCTCAGCCGAACCCGGGATTTGAGTTGCCTGCTCAGTGCGTCAGCAAGGACTTTCCTGTTCAAGCGGGAAATCGCACAGTGGGTACCGGCCGCTTCGCCCGTGATGCTGGTGACTACAAAGGTTTTCGCATTTCCCCGGAGCAGGATGCTATGACGCTTGCATTCGCTGCCGCGGCCATTGAAAATATGCAATTAGGGAAGCAGGCCCAGACCGATATTATTAGCATTGGACTGAGCGCTACGGATTACGTGGGACACACCTTCGGCACGGAGGGTACGGAGAGTTGCATCCAAGTGGATCGTTTAGACACGGAGCTTGGTGCATTCTTTGATAAACTGGATAAGGATGGGATTGACTACGTAGTAGTGCTGACTGCAGATCATGGAGGACACGATCTGCCCGAACGTCATCGTATGAATGCCATGCCGATGGAACAGCGCGTAGACATGGCCCTGACACCTAAAGCTCTGAATGCTACCATCGCTGAGAAAGCTGGCCTTCCGGGCAAAAAGGTTATTTGGTCAGATGGACCTTCTGGCGATATTTACTATGATAAGGGCCTTACAGCCGCTCAACGTGCCCGTGTTGAAACCGAGGCGTTAAAATACTTGCGCGCGCATCCCCAAGTACAGACTGTATTCACTAAGGCGGAAATCGCGGCTACCCCTTCTCCGTCGGGACCACCTGAGAGCTGGAGTTTGATCCAGGAAGCTCGCGCGTCATTTTACCCGTCGCGCTCCGGGGACCTGTTACTTTTATTGAAACCTCGTGTGATGAGCATTCCTGAGCAAGCAGTCATGGGCTCGGTTGCAACCCATGGATCTCCATGGGATACGGATCGCCGTGTGCCTATCCTGTTTTGGCGCAAAGGTATGCAGCATTTCGAACAACCCTTAGGAGTAGAGACTGTTGATATTTTGCCCTCCTTGGCTGCACTTATTAAGCTTCCTGTTCCTAAGGATCAGATCGACGGCCGCTGTCTGGACTTGGTCGCCGGCAAGGATGATTCCTGTGCTGGACAGGGA"
)

# ============================================================
# MAIN LOOP
# ============================================================
overall_start = time.time()

for num in range(2,9):

    print(f"\n[INFO] Running for num variants: {num}")

    args.num_proteins = num

    mutreg_nt = mutreg_nt_full
    sequence_nt = UPSTREAM_NT + mutreg_nt + DOWNSTREAM_NT

    # ---- Step 1: Build primers and graph ----
    print("[STEP 1] Creating primer DataFrame...")
    primer_df = create_primer_df(sequence_nt, args)

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
    "protein_name": "SPAP",
    "num_variants": num,
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
    "ilp_path_length": sum(len(path) for path in ilp_res.paths),
    "ilp_status": ilp_res.status,

    # ---- Greedy ----
    "greedy_objective": greedy_obj,
    "greedy_time_sec": round(greedy_time, 3),
    "greedy_path_length": sum(len(path) for path in greedy_solution)
    }

    csv_path = output_dir / f"results_{num}_variants.csv"
    pd.DataFrame([results]).to_csv(csv_path, index=False)
    print(f"Saved CSV summary: {csv_path}")

    # ---- Step 5: Save JSON paths ----
    paths_out = {
        "protein_name": "SPAP",
        "num_variants": num,
        "ilp_path": ilp_res.paths,
        "greedy_path": greedy_solution,
    }

    json_path = output_dir / f"paths_{num}_variants.json"
    with open(json_path, "w") as f:
        json.dump(paths_out, f, indent=2)
    print(f"Saved JSON paths: {json_path}")


overall_end = time.time()
total_minutes = (overall_end - overall_start) / 60
print(f"\n[ALL DONE] All sequence lengths processed successfully.")
print(f"[TOTAL RUNTIME] {total_minutes:.2f} minutes ({overall_end - overall_start:.1f} seconds total)")