from General.primer_graphs import *
from PD_var_ILP.greedy import *
from PD_var_ILP.ilp_model import *
import time
import pandas as pd
from pathlib import Path
import ast

def run_var_ilp(sequence_nt, mutreg_nt, protein_name, args,cfg):
    """Run primer design with ILP and greedy methods, log performance, and save outputs."""

    # ---- Step 1: Build primers and graph ----
    primer_df = create_primer_df(sequence_nt, args, cfg)

    t0 = time.time()
    graph = create_graph(primer_df, len(mutreg_nt), args)
    graph_time = time.time() - t0

    # ---- Step 2: Greedy solution ----
    t1 = time.time()
    greedy_solution, greedy_obj = run_greedy(graph, primer_df, args)
    greedy_time = time.time() - t1

    # ---- Step 3: ILP model ----
    ilp_res: ILPResult = ilp_model(graph, sequence_nt, mutreg_nt, args, cfg)

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

    print(f"Saved summary to: {csv_path}")

    # ---- Step 5: Save selected primers as CSVs (instead of JSON paths) ----
    def paths_to_primers_csv(paths, method_name: str):
        """
        Convert list-of-paths (each path is list of primer nodes) into a flat CSV:
        protein_name, variant_index, primer_index, start, stop, orientation, seq, efficiency
        """
        rows = []

        for variant_idx, path_nodes in enumerate(paths):
            if not path_nodes:
                continue

            path_nodes = [
                ast.literal_eval(n) if isinstance(n, str) else n
                for n in path_nodes
            ]

            # path_nodes should match primer_df index entries
            df_sel = primer_df.loc[path_nodes].copy().reset_index()

            # rename fr → orientation if exists
            if "fr" in df_sel.columns:
                df_sel = df_sel.rename(columns={"fr": "orientation"})

            # add identifiers
            df_sel.insert(0, "protein_name", protein_name)
            df_sel.insert(1, "variant_index", variant_idx)
            df_sel.insert(2, "primer_index", range(len(df_sel)))

            # pick column order (keep only those that exist)
            wanted = ["protein_name", "variant_index", "primer_index",
                      "start", "stop", "orientation", "seq", "efficiency"]
            cols = [c for c in wanted if c in df_sel.columns]
            df_sel = df_sel[cols]

            rows.append(df_sel)

        out_path = out_dir / f"var_ILP_selected_primers_{method_name}.csv"
        if rows:
            pd.concat(rows, ignore_index=True).to_csv(out_path, index=False)
        else:
            pd.DataFrame(columns=["protein_name","variant_index","primer_index",
                                  "start","stop","orientation","seq","efficiency"]).to_csv(out_path, index=False)

        print(f"Saved {method_name} selected primers to: {out_path}")

    # ILP primers
    paths_to_primers_csv(ilp_res.paths, "ILP")

    # Greedy primers
    paths_to_primers_csv(greedy_solution, "greedy")


    return results_df, ilp_res.paths, greedy_solution



