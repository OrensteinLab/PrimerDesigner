import os
import time
import json
import pandas as pd
from pathlib import Path
from PD_mul_ILP.create_graphs import *
from PD_mul_ILP.ilp_model import *
from General.args import *

import sys

sys.argv = [
        sys.argv[0],
        "--file_path", "data/10_protein_coding_sequences_example.txt",
        "--output", "Experiment_results/mul_ilp",
        ]

args = get_args()

def main():

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
    overall_start = time.time()
    # ============================================================
    # MAIN LOOP: number of proteins 2–10
    # ============================================================
    summary_rows = []  # <-- collect rows for a final aggregated CSV
    for i in range(2,len(all_protein_names)):
        print(f"\n[INFO] Processing {i} protein(s)...")

        # Subset of first i proteins
        mutreg_regions = all_mutreg_regions[:i]
        sequences_nt = all_full_sequences[:i]
        protein_names = all_protein_names[:i]

        # --------------------------------------------------------
        # GRAPH CREATION
        # --------------------------------------------------------
        print(f"[STEP] Creating graphs for {i} proteins...")
        graphs, graph_time, graph_memory, primer_dfs = create_graphs(
            mutreg_regions, sequences_nt, protein_names, args
        )
        print(f"[DONE] Graphs created in {graph_time:.2f} sec (peak {graph_memory:.1f} MB).")

        
        print("Finding forbidden pairs across proteins...")
        t_forbid = time.time()
        single_forbidden, multiple_forbidden, single_pair_cnt, multi_pairs_cnt = find_forbidden_pairs(
            protein_names, sequences_nt, args
        )
        forbidden_time = time.time() - t_forbid
        print(f"[DONE] Forbidden pairs in {forbidden_time:.2f} sec "
              f"(single={single_pair_cnt}, inter={multi_pairs_cnt}).")


        # --------------------------------------------------------
        # ILP OPTIMIZATION
        # --------------------------------------------------------
        print(f"[STEP] Running ILP model...")
        ilp_res = run_ilp(single_forbidden, multiple_forbidden, protein_names, graphs)
        print(
            f"[DONE] ILP finished with status={ilp_res.status}, "
            f"objective={ilp_res.objective:.2f}, time={ilp_res.ilp_time:.2f} sec."
        )

        del single_forbidden, multiple_forbidden

        # ============================================================
        # SAVE RESULTS (CSV + JSON per iteration)
        # ============================================================

        # ---- CSV ----
        csv_path = output_dir / f"results_{i:02d}_proteins.csv"
        results = {
            "num_proteins": len(protein_names),
            "graph_time_sec": graph_time,
            "graph_peak_mem_MB": graph_memory,
            "ilp_num_vars": ilp_res.num_vars,
            "ilp_num_constraints": ilp_res.num_constraints,
            "ilp_single_forbidden_cnt": single_pair_cnt,
            "ilp_inter_forbidden_cnt": multi_pairs_cnt,
            "forbidden_time_sec": forbidden_time,    
            "ilp_setup_time_sec": ilp_res.setup_time,
            "ilp_setup_peak_mem_MB": ilp_res.setup_memory,
            "ilp_optimize_time_sec": ilp_res.ilp_time,
            "ilp_optimize_peak_mem_MB": ilp_res.ilp_memory,
            "ilp_objective": ilp_res.objective,
            "ilp_path_length": sum(len(path) for path in ilp_res.protein_paths.values()),
            "ilp_status": ilp_res.status,
        }

        pd.DataFrame([results]).to_csv(csv_path, index=False)
        print(f"[SAVE] CSV saved → {csv_path}")

        summary_rows.append(results.copy())

        # ---- JSON ----
        json_path = output_dir / f"paths_{i:02d}_proteins.json"

        def to_jsonable_paths(paths_dict):
            out = {}
            for name, path in (paths_dict or {}).items():
                out[name] = [
                    node if isinstance(node, str) else list(node)
                    for node in path
                ]
            return out

        combined_paths = {
            "ILP_paths": to_jsonable_paths(ilp_res.protein_paths),
        }

        with open(json_path, "w") as f:
            json.dump(combined_paths, f, indent=2)
        print(f"[SAVE] JSON saved → {json_path}")

        print(f"[DONE] Iteration {i} complete.\n")

    # ============================================================
    # FINAL SUMMARY CSV (all i)
    # ============================================================
    final_df = pd.DataFrame(summary_rows)
    final_csv = output_dir / "results_all_proteins.csv"
    final_df.to_csv(final_csv, index=False)
    print(f"[SAVE] Aggregated CSV saved → {final_csv}")

    overall_end = time.time()
    total_minutes = (overall_end - overall_start) / 60
    print("[INFO] All iterations complete.")
    print(f"[TOTAL RUNTIME] {total_minutes:.2f} minutes ({overall_end - overall_start:.1f} seconds total)")

if __name__ == '__main__':
    mp.set_start_method("spawn", force=True)
    main()