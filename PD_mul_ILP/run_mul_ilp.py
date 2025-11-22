import json
import time
from pathlib import Path
import pandas as pd

from PD_mul_ILP.create_graphs import *
from PD_mul_ILP.ilp_model import *

def run_mul_ilp(mutreg_regions, sequences_nt, protein_names, args):
    """
    Runs graph creation, forbidden-pairs discovery, greedy baseline, and ILP;
    appends a single-row summary to <args.output>/mul_ilp_results.csv and
    writes per-run paths to <args.output>/paths_N{num_proteins}.json.
    """
    # --------------------------------------------------------
    # OUTPUT SETUP
    # --------------------------------------------------------
    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    # --------------------------------------------------------
    # GRAPH CREATION
    # --------------------------------------------------------
    print(f"[STEP] Creating graphs for {len(protein_names)} proteins...")
    graphs, graph_time, graph_memory, primer_dfs = create_graphs(
        mutreg_regions, sequences_nt, protein_names, args
    )
    print(f"[DONE] Graphs created in {graph_time:.2f} sec (peak {graph_memory:.1f} MB).")

    # --------------------------------------------------------
    # FORBIDDEN PAIRS
    # --------------------------------------------------------
    print("[STEP] Finding forbidden pairs across proteins...")
    t_forbid = time.time()
    single_forbidden, multiple_forbidden, single_pair_cnt, multi_pairs_cnt = find_forbidden_pairs(
        protein_names, sequences_nt, args
    )
    forbidden_time = time.time() - t_forbid
    print(f"[DONE] Forbidden pairs in {forbidden_time:.2f} sec "
          f"(intra={single_pair_cnt}, inter={multi_pairs_cnt}).")

    # --------------------------------------------------------
    # ILP OPTIMIZATION
    # --------------------------------------------------------
    print("[STEP] Running ILP model...")
    ilp_res = run_ilp(single_forbidden, multiple_forbidden, protein_names, graphs)
    print(f"[DONE] ILP status={ilp_res.status}, objective={ilp_res.objective:.2f}, "
          f"time={ilp_res.ilp_time:.2f} sec.")

    # free big sets if needed
    del single_forbidden, multiple_forbidden

    # ============================================================
    # SAVE RESULTS
    # ============================================================
    # ---- CSV (append one row per run) ----
    csv_path = outdir / "mul_ilp_results.csv"
    row = {
        "num_proteins": len(protein_names),
        "graph_time_sec": graph_time,
        "ilp_num_vars": getattr(ilp_res, "num_vars", None),
        "ilp_num_constraints": getattr(ilp_res, "num_constraints", None),
        "ilp_intra_forbidden_cnt": single_pair_cnt,
        "ilp_inter_forbidden_cnt": multi_pairs_cnt,
        "forbidden_time_sec": forbidden_time,
        "ilp_setup_time_sec": getattr(ilp_res, "setup_time", None),
        "ilp_optimize_time_sec": getattr(ilp_res, "ilp_time", None),
        "ilp_feasibility": "FEASIBLE" if ilp_res.status == 2 else "INFEASIBLE",
        "total_primer_efficiency": getattr(ilp_res, "objective", None), 
        "num_primers": sum(len(path) for path in (ilp_res.protein_paths or {}).values()),
    }
    df = pd.DataFrame([row])
    header = not csv_path.exists()
    df.to_csv(csv_path, mode="a", index=False, header=header)
    print(f"[SAVE] Appended summary → {csv_path}")

    # ---- JSON (paths per run) ----
    def to_jsonable_paths(paths_dict):
        out = {}
        for name, path in (paths_dict or {}).items():
            out[name] = [node if isinstance(node, str) else list(node) for node in path]
        return out

    json_path = outdir / f"primers_N{len(protein_names)}.json"
    combined_paths = {
        "ILP_primers": to_jsonable_paths(getattr(ilp_res, "protein_paths", {})),
    }
    with open(json_path, "w") as f:
        json.dump(combined_paths, f, indent=2)
    print(f"[SAVE] Paths saved → {json_path}")

    return df,ilp_res.protein_paths