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

# ---- Build primer table and graph ----
primer_df = create_primer_df(sequence_nt, args)

t_graph0 = time.time()
tracemalloc.start()
graph = create_graph(primer_df, len(mutreg_nt), args)
graph_time = time.time() - t_graph0
graph_peak_mb = tracemalloc.get_traced_memory()[1] / 1e6
tracemalloc.stop()

longest_path_t0 = time.time()

full_path = longest_path_dag(graph, 's', 'd')
# strip s/d for primer nodes only (they should match primer_df index)
primer_path_nodes = full_path[1:-1]
# cost = sum of node weights as you stored them on edges entering those nodes
# but you already stored the per-primer cost in primer_df, so:
primer_designer_set = primer_df.loc[primer_path_nodes].copy().reset_index()
primer_designer_efficiency = float(primer_designer_set['efficiency'].sum())

longest_path_time = time.time() - longest_path_t0

total_time = time.time() - t0

primer_locations = [(31,55,'f'),(211,230,'r'),(183,208, 'f'),(357,382,'r'),(322,343,'f'),(499,518,'r'),(471,492,'f'),
                (644,669,'r'),(612,638,'f'),(784,808,'r'),(752,773,'f'),(925,947,'r'),(892,919,'f'),(1066,1088,'r'),(988,1009,'f'),
                (1156,1180,'r'),(1132,1154,'f'),(1303,1324,'r'),(1261,1281,'f'),(1432,1459,'r'),(1406,1427,'f'),(1563,1586,'r'),
                (1529,1550,'f'),(1639,1662,'r'),(1613,1637,'f'),(1730,1751,'r')]

quick_primers = []

for primer in primer_locations:
    start, end, strand = primer
    quick_primers.append((start-len(UPSTREAM_NT), end-len(UPSTREAM_NT), strand))

quick_set = primer_df.loc[quick_primers].copy().reset_index()
quick_efficiency = float(quick_set['efficiency'].sum())

# ---- CSV (no path in CSV) ----
row = {
    "protein_name": protein_name,
    "graph_nodes": len(graph.nodes),
    "graph_edges": len(graph.edges),
    "graph_time_sec": round(graph_time, 3),
    "graph_peak_mem_MB": round(graph_peak_mb, 1),
    "PD_single_efficiency": primer_designer_efficiency,
    "QuickChange_efficiency": quick_efficiency,
    "PD_single_primers":len(primer_path_nodes),
    "QuickChange_primers":len(quick_primers),
    "shortest_path_time": longest_path_time,
    "total_time_sec": round(total_time, 3)
}
df = pd.DataFrame([row])

out_base = Path("QucikChange_comparison")
out_base.parent.mkdir(parents=True, exist_ok=True)
csv_path = out_base.with_suffix(".csv")
df.to_csv(csv_path, index=False)


