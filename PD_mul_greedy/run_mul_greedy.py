from __future__ import annotations

from pathlib import Path
import time
import json
import pandas as pd
import networkx as nx

from General.primer_graphs import create_primer_df, create_graph
from General.primer_data import *
from General.utils import *

MUTREG_START = len(UPSTREAM_NT)

def run_mul_greedy(
    sequences_nt: list[str],
    mutreg_regions: list[str],
    protein_names: list[str],
    args,
) -> None:
    """
    Greedy, multi–non-homologous mode:
      • Build per-protein graphs and iterate until a valid (conflict-free) path is found
      • Resolve cross-hybridizations by removing offending primer nodes
      • Save outputs under args.output/ as:
          - summary.csv (one row)
          - per_protein_metrics.csv (one row per protein)
          - paths.json (all accepted paths)
    """
    # Ensure output is a directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    (paths,
     total_efficiency,
     cross_cnt,
     protein_cnt,
     re_iters,
     unresolved,
     per_protein_rows) = run_greedy(sequences_nt, mutreg_regions, protein_names, args)
    greedy_time = time.time() - t0

    # ---------- overall summary ----------
    summary_row = {
        "num_proteins": len(protein_names),
        "total_primer_efficiency": total_efficiency,
        "greedy_time_sec": round(greedy_time, 3),
        "cross_hybridizations_cnt": cross_cnt,
        "proteins_with_reiterations_cnt": protein_cnt,
        "total_reiterations": re_iters,
        "unresolved_proteins_cnt": len(unresolved),
        "unresolved_proteins": ",".join(unresolved)
    }
    summary_df = pd.DataFrame([summary_row])
    summary_df.to_csv(output_dir / "PD_mul_Greedy_summary.csv", index=False)

    # ---------- per-protein metrics ----------
    # per_protein_rows already contains per-protein timings & metadata
    if per_protein_rows:
        pd.DataFrame(per_protein_rows).to_csv(output_dir / "per_protein_metrics.csv", index=False)

    # ---------- paths ----------
    paths_json_path = output_dir / "primers_per_protein.json"
    out = {
        "paths": paths
    }
    with open(paths_json_path, "w") as f:
        json.dump(out, f, indent=2)

    print(f" Saved summary to: {output_dir/'summary.csv'}")
    print(f" Saved per-protein metrics to: {output_dir/'per_protein_metrics.csv'}")
    print(f" Saved paths to: {paths_json_path}")

    return summary_df, paths


def run_greedy(
    sequences_nt: list[str],
    mutreg_regions: list[str],
    protein_names: list[str],
    args,
):
    """
    Returns:
        paths: dict[str, list[tuple[int,int,str]]]
        total_efficiency: float
        cross_cnt: int
        protein_cnt: int
        re_iterations: int
        unresolved_proteins: list[str]
        per_protein_rows: list[dict]  # timing & metadata per protein
    """
    paths: dict[str, list[tuple[int, int, str]]] = {}
    selected_primers: list[str] = []   # sequences accepted so far (S)
    total_efficiency = 0.0
    proteins_with_retries = 0
    cross_cnt = 0
    re_iterations = 0
    unresolved_proteins: set[str] = set()

    per_protein_rows: list[dict] = []

    n_total = len(sequences_nt)

    for i, (seq_nt, mutreg_nt, protein_name) in enumerate(
        zip(sequences_nt, mutreg_regions, protein_names)
    ):
        if i % 10 == 0:
            print(f"[INFO] Processing protein {i}/{n_total}: {protein_name}")

        # 1) Build primer DataFrame (time it)
        t_df0 = time.time()
        primer_df = create_primer_df(seq_nt, args)
        primer_df_time = time.time() - t_df0

        # 2) Build graph (time it)
        t_g0 = time.time()
        graph = create_graph(primer_df, len(mutreg_nt), args)
        graph_time = time.time() - t_g0

        # 3) Iterate until a non-conflicting path is found (or none exists)
        needs_retry = True
        iters_for_this_protein = 0
        path_found = False

        while needs_retry:
            if iters_for_this_protein > 0:
                re_iterations += 1
            iters_for_this_protein += 1

            try:
                lp = longest_path_dag(graph, 's', 'd')[1:-1]
            except nx.NetworkXNoPath:
                print(f"[WARN] No valid path found for {protein_name}. Skipping.")
                break

            primer_seqs = {
                node: seq_nt[node  [0] + MUTREG_START : node[1] + MUTREG_START] for node in lp
            }

            violating_nodes = set()
            for node, pseq in primer_seqs.items():
                for other_seq in selected_primers:
                    tm = calc_tm(other_seq, pseq)
                    if tm >= MAX_TM:
                        cross_cnt += 1
                        violating_nodes.add(node)
                        break

            if violating_nodes:
                graph.remove_nodes_from(violating_nodes)
                needs_retry = True
            else:
                selected_primers.extend(primer_seqs.values())
                paths[protein_name] = lp
                path_found = True
                path_efficiency = primer_df.loc[lp, 'efficiency'].sum()
                total_efficiency += float(path_efficiency)
                needs_retry = False

        if iters_for_this_protein > 1:
            proteins_with_retries += 1
        if not path_found:
            unresolved_proteins.add(protein_name)

        # Record per-protein metrics row
        per_protein_rows.append({
            "protein_name": protein_name,
            "seq_length": len(seq_nt),
            "mutreg_length": len(mutreg_nt),
            "primer_df_time_sec": round(primer_df_time, 4),
            "graph_time_sec": round(graph_time, 4),
            "iterations": iters_for_this_protein,
            "resolved": bool(path_found),
            "num_primers": (len(paths[protein_name]) if path_found else 0),
            "primer_efficiency": (float(primer_df.loc[paths[protein_name], 'efficiency'].sum())
                                if path_found else 0.0),
        })

    return (
        paths,
        total_efficiency,
        cross_cnt,
        proteins_with_retries,
        re_iterations,
        list(unresolved_proteins),
        per_protein_rows,
    )