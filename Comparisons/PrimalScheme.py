import time
import json
import tracemalloc
from pathlib import Path
import pandas as pd
import networkx as nx
from General.primer_graphs import create_primer_df, create_graph
from General.utils import *
from General.args import *
import sys

mutreg_nt = (
    "ATGCAAAGCCCAGCACCTGCCGCAGCGCCTGCCCCTGCGGCACGTTCCATCGCAGCTACGCCTCCTAAACTGATCGTGGCAATTAGCGTGGACCAGTTTAGTGCAGACTTGTTCTCGGAGTATCGTCAATATTACACCGGAGGTTTAAAGCGTCTTACATCCGAAGGAGCTGTGTTCCCACGTGGTTATCAGAGTCATGCGGCAACAGAAACGTGTCCTGGTCACTCAACGATCCTGACAGGATCACGTCCGTCACGTACGGGTATTATCGCTAATAACTGGTTCGACTTGGACGCAAAGCGTGAGGATAAAAATCTGTACTGTGCTGAGGATGAATCCCAACCCGGTAGTTCGTCTGACAAGTACGAAGCTTCGCCACTGCACTTAAAGGTACCCACCCTGGGGGGACGCATGAAAGCCGCCAATCCTGCGACTCGTGTCGTCTCTGTTGCCGGCAAGGATCGCGCGGCCATTATGATGGGTGGCGCCACAGCGGATCAGGTCTGGTGGTTAGGGGGGCCTCAGGGGTATGTTTCGTATAAGGGTGTAGCGCCAACTCCCCTTGTAACACAGGTCAATCAGGCCTTTGCACAGCGCTTAGCTCAGCCGAACCCGGGATTTGAGTTGCCTGCTCAGTGCGTCAGCAAGGACTTTCCTGTTCAAGCGGGAAATCGCACAGTGGGTACCGGCCGCTTCGCCCGTGATGCTGGTGACTACAAAGGTTTTCGCATTTCCCCGGAGCAGGATGCTATGACGCTTGCATTCGCTGCCGCGGCCATTGAAAATATGCAATTAGGGAAGCAGGCCCAGACCGATATTATTAGCATTGGACTGAGCGCTACGGATTACGTGGGACACACCTTCGGCACGGAGGGTACGGAGAGTTGCATCCAAGTGGATCGTTTAGACACGGAGCTTGGTGCATTCTTTGATAAACTGGATAAGGATGGGATTGACTACGTAGTAGTGCTGACTGCAGATCATGGAGGACACGATCTGCCCGAACGTCATCGTATGAATGCCATGCCGATGGAACAGCGCGTAGACATGGCCCTGACACCTAAAGCTCTGAATGCTACCATCGCTGAGAAAGCTGGCCTTCCGGGCAAAAAGGTTATTTGGTCAGATGGACCTTCTGGCGATATTTACTATGATAAGGGCCTTACAGCCGCTCAACGTGCCCGTGTTGAAACCGAGGCGTTAAAATACTTGCGCGCGCATCCCCAAGTACAGACTGTATTCACTAAGGCGGAAATCGCGGCTACCCCTTCTCCGTCGGGACCACCTGAGAGCTGGAGTTTGATCCAGGAAGCTCGCGCGTCATTTTACCCGTCGCGCTCCGGGGACCTGTTACTTTTATTGAAACCTCGTGTGATGAGCATTCCTGAGCAAGCAGTCATGGGCTCGGTTGCAACCCATGGATCTCCATGGGATACGGATCGCCGTGTGCCTATCCTGTTTTGGCGCAAAGGTATGCAGCATTTCGAACAACCCTTAGGAGTAGAGACTGTTGATATTTTGCCCTCCTTGGCTGCACTTATTAAGCTTCCTGTTCCTAAGGATCAGATCGACGGCCGCTGTCTGGACTTGGTCGCCGGCAAGGATGATTCCTGTGCTGGACAGGGA"
)
sequence_nt = UPSTREAM_NT + mutreg_nt + DOWNSTREAM_NT
protein_name = 'SPAP'

sys.argv = [
        sys.argv[0],
        "--file_path", "input_path",
        "--output", "output_path",
        ]

args = get_args()

t0 = time.time()

# ---- Build primer table--
primer_df = create_primer_df(sequence_nt, args)
primal_scheme_primers = []

with open("Comparisons/PrimalScheme_primer.bed") as f:
    for line in f:
        # Skip comment or empty lines
        if line.startswith("#") or not line.strip():
            continue

        fields = line.strip().split("\t")

        # Extract coordinates and direction
        start = int(fields[1])
        end = int(fields[2])
        strand = fields[6]

        # Determine primer direction
        direction = 'f' if strand == '+' else 'r'

        # Append as tuple
        primal_scheme_primers.append((start-len(UPSTREAM_NT), end-len(UPSTREAM_NT), direction))

primal_scheme_set = primer_df.loc[primal_scheme_primers].copy().reset_index()
primal_scheme_efficiency = float(primal_scheme_set['efficiency'].sum())

t_graph0 = time.time()
tracemalloc.start()
graph = create_graph(primer_df, len(mutreg_nt), args)
graph_time = time.time() - t_graph0
graph_peak_mb = tracemalloc.get_traced_memory()[1] / 1e6
tracemalloc.stop()

longest_path_t0 = time.time()

full_path = longest_path_dag(graph, 's', 'd')
# strip s/d for primer nodes only (they should match primer_df index)
PD_single_path = full_path[1:-1]
# cost = sum of node weights as you stored them on edges entering those nodes
# but you already stored the per-primer cost in primer_df, so:
primer_set = primer_df.loc[PD_single_path].copy().reset_index()
PD_single_efficiency = float(primer_set['efficiency'].sum())

longest_path_time = time.time() - longest_path_t0
total_time = time.time() - t0


# ---- CSV (no path in CSV) ----
row = {
    "protein_name": protein_name,
    "graph_nodes": len(graph.nodes),
    "graph_edges": len(graph.edges),
    "graph_time_sec": round(graph_time, 3),
    "graph_peak_mem_MB": round(graph_peak_mb, 1),
    "PD_single_efficiency": PD_single_efficiency,
    "PrimalScheme_efficiency": primal_scheme_efficiency,
    "PD_single_primers":len(PD_single_path),
    "PrimalScheme_primers":len(primal_scheme_primers),
    "longest_path_time": longest_path_time,
    "total_time_sec": round(total_time, 3)
}
df = pd.DataFrame([row])

out_base = Path("PrimalScheme_comparison")
out_base.parent.mkdir(parents=True, exist_ok=True)
csv_path = out_base.with_suffix(".csv")
df.to_csv(csv_path, index=False)



