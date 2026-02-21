import time
import json
import pandas as pd
from pathlib import Path
from PD_mul_ILP.create_graphs import *
from PD_mul_ILP.ilp_model import *
from General.args import *


def main():

    args = get_args()

    args.output = "Results"

    args.file_path = "data/10_protein_coding_sequences.txt"

    cfg = GU.load_config()

    # Create output directory if not exists
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"[INFO] Output directory: {output_dir.resolve()}")


    # load sequences
    print(f"[INFO] Reading protein coding sequences from: {args.file_path}")
    all_mutreg_regions, all_full_sequences, all_protein_names = read_sequences(args.file_path, cfg)
    print(f"[INFO] Total proteins loaded: {len(all_protein_names)}")
    overall_start = time.time()

    # main loop
    summary_rows = []  
    for i in range(1,len(all_protein_names)+1):
        print(f"\n[INFO] Processing {i} protein(s)...")

        # Subset of first i proteins
        mutreg_regions = all_mutreg_regions[:i]
        sequences_nt = all_full_sequences[:i]
        protein_names = all_protein_names[:i]


        # graph creation
        print(f"[STEP] Creating graphs for {i} proteins...")
        graphs, graph_time, graph_memory = create_graphs(
            mutreg_regions, sequences_nt, protein_names, args, cfg
        )
        print(f"[DONE] Graphs created in {graph_time:.2f} sec (peak {graph_memory:.1f} MB).")

        
        print("Finding forbidden pairs across proteins...")
        t_forbid = time.time()
        single_forbidden, multiple_forbidden, single_pair_cnt, multi_pairs_cnt = find_forbidden_pairs(
            protein_names, sequences_nt, args,cfg
        )
        forbidden_time = time.time() - t_forbid
        print(f"[DONE] Forbidden pairs in {forbidden_time:.2f} sec "
              f"(single={single_pair_cnt}, inter={multi_pairs_cnt}).")


        # ILP optimization
        print(f"[STEP] Running ILP model...")
        ilp_res = run_ilp(single_forbidden, multiple_forbidden, protein_names, graphs)
        print(
            f"[DONE] ILP finished with status={ilp_res.status}, "
            f"objective={ilp_res.objective:.2f}, time={ilp_res.ilp_time:.2f} sec."
        )

        del single_forbidden, multiple_forbidden

        # save results for this iteration
        # ---- CSV ----
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

        summary_rows.append(results.copy())

    # FINAL SUMMARY CSV (all i)
    final_df = pd.DataFrame(summary_rows)
    final_csv = output_dir / "PD-mul-ILP.csv"
    final_df.to_csv(final_csv, index=False)
    print(f"[SAVE] Aggregated CSV saved → {final_csv}")

    overall_end = time.time()
    total_minutes = (overall_end - overall_start) / 60
    print("[INFO] All iterations complete.")
    print(f"[TOTAL RUNTIME] {total_minutes:.2f} minutes ({overall_end - overall_start:.1f} seconds total)")

if __name__ == '__main__':
    mp.set_start_method("spawn", force=True)
    main()