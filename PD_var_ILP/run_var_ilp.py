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
    tracemalloc.start()
    graph = create_graph(primer_df, len(mutreg_nt), args)
    graph_time = time.time() - t0
    graph_peak_mb = tracemalloc.get_traced_memory()[1] / 1e6
    tracemalloc.stop()

    # ---- Step 2: Greedy solution ----
    t1 = time.time()
    greedy_solution, greedy_obj = run_greedy(graph, primer_df, args)
    greedy_time = time.time() - t1

    # ---- Step 3: ILP model ----
    ilp_res: ILPResult = ilp_model(graph, sequence_nt, mutreg_nt, args)

    # ---- Step 4: CSV summary (no large paths) ----
    results = {
    "protein_name": "SPAP",
    "seq_length": len(mutreg_nt),
    "graph_nodes": len(graph.nodes),
    "graph_edges": len(graph.edges),
    "graph_time_sec": round(graph_time, 3),
    "graph_peak_mem_MB": round(graph_peak_mb, 1),
    "greedy_path_length": sum(len(path) for path in greedy_solution),

    # ---- ILP (updated field names) ----
    "ilp_num_vars": ilp_res.num_vars,
    "ilp_path_length": sum(len(path) for path in ilp_res.paths),
    "ilp_num_constraints": ilp_res.num_constraints,
    "ilp_setup_time_sec": round(ilp_res.setup_time, 3),
    "ilp_setup_peak_mem_MB": round(ilp_res.setup_peak_mem_mb, 1),
    "ilp_optimize_time_sec": round(ilp_res.optimize_time, 3),
    "ilp_optimize_peak_mem_MB": round(ilp_res.optimize_peak_mem_mb, 1),
    "ilp_objective": ilp_res.objective,
    "ilp_status": ilp_res.status,

    # ---- Greedy ----
    "greedy_objective": greedy_obj,
    "greedy_time_sec": round(greedy_time, 3),
    }


    out_base = Path(args.output)
    pd.DataFrame([results]).to_csv(f"{args.output}.csv", index=False)

    # ---- Step 5: JSON paths only ----
    paths_out = {
        "protein_name": protein_name,
        "num_variants": getattr(args, "num_proteins", None),
        "ilp_path": ilp_res.paths,
        "greedy_path": greedy_solution,
    }

    with open(out_base.with_suffix('').as_posix() + "_paths.json", "w") as f:
        json.dump(paths_out, f, indent=2)

    print(f"✅ Saved summary to: {out_base.with_suffix('.csv')}")
    print(f"✅ Saved paths to: {out_base.with_suffix('.json')}")

