from General.primer_graphs import *
from PD_var_ILP.greedy import *
from PD_var_ILP.ilp_model import *
import time
import tracemalloc
import pandas as pd
import json
from pathlib import Path


def run_var_ilp(sequence_nt, mutreg_nt, protein_name, args):
    """Run primer design with ILP and greedy methods, log performance, and save outputs."""

    # ---- Step 1: Build primers and graph ----
    primer_df = create_primer_df(sequence_nt, args)

    t0 = time.time()
    graph = create_graph(primer_df, len(mutreg_nt), args)
    graph_time = time.time() - t0

    # ---- Step 2: Greedy solution ----
    t1 = time.time()
    greedy_solution, greedy_obj = run_greedy(graph, primer_df, args)
    greedy_time = time.time() - t1

    # ---- Step 3: ILP model ----
    ilp_res: ILPResult = ilp_model(graph, sequence_nt, mutreg_nt, args)

    # ---- Step 4: CSV summary (no large paths) ----
    results = {
    "protein_name": protein_name,
    "num_variants": getattr(args, "num_proteins", None),
    "seq_length": len(mutreg_nt),
    "graph_nodes": len(graph.nodes),
    "graph_edges": len(graph.edges),
    "graph_time_sec": round(graph_time, 3),
    "greedy_path_length": sum(len(path) for path in greedy_solution),

    # ---- ILP (updated field names) ----
    "ilp_num_vars": ilp_res.num_vars,
    "ilp_path_length": sum(len(path) for path in ilp_res.paths),
    "ilp_num_constraints": ilp_res.num_constraints,
    "ilp_setup_time_sec": round(ilp_res.setup_time, 3),
    "ilp_optimize_time_sec": round(ilp_res.optimize_time, 3),
    "ILP_primer_efficiency": ilp_res.objective,
    "ILP_solution_feasibility": "FEASIBLE" if ilp_res.objective != None else "INFEASIBLE",

    # ---- Greedy ----
    "greedy_solution_feasibility": "FEASIBLE" if greedy_obj!= None else "INFEASIBLE",
    "greedy_primer_efficiency": greedy_obj,
    "greedy_time_sec": round(greedy_time, 3),
    }


    # Ensure output directory exists
    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ---- Step 4: Save CSV ----
    csv_path = out_dir / "var_ILP_results.csv"

    results_df = pd.DataFrame([results])
    results_df.to_csv(csv_path, index=False)

    # ---- Step 5: Save JSON paths ----
    paths_out = {
        "protein_name": protein_name,
        "num_variants": getattr(args, "num_proteins", None),
        "ilp_path": ilp_res.paths,
        "greedy_path": greedy_solution,
    }

    json_path = out_dir / "primer_paths.json"
    with open(json_path, "w") as f:
        json.dump(paths_out, f, indent=2)

    # ---- Step 6: Print confirmation ----
    print(f"Saved summary to: {csv_path}")
    print(f"Saved paths to: {json_path}")

    return results_df, ilp_res.paths, greedy_solution



