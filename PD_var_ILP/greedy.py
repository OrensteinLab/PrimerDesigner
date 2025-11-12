import pandas as pd
import itertools as it
import networkx as nx
from General.utils import *

def run_greedy(graph, primer_df,args):

  print("Running greedy algorithm")

  path_ls = [[]]
  graph_sub = graph.copy()
  for i in range(args.num_proteins):
    nodes_to_remove = []
    for p,n in it.product(path_ls[-1],graph_sub.nodes()):
      if n=='s' or n=='d':
        continue
      p_start, p_end, _ = p
      n_start, n_end, _ = n
      overlap = max(0, min(p_end, n_end) - max(p_start, n_start))
      pn_intersect = (
        overlap > args.allowed_overlap and
        (n[2] == p[2])  # same strand
      )
      if pn_intersect:
        nodes_to_remove.append(n)
    graph_sub.remove_nodes_from(set(nodes_to_remove))
    try:
      path_ls.append([primer for primer in longest_path_dag(graph_sub, 's', 'd')[1:-1]])
    except:
      print(f'WARNING: No feasible primer sequence for lib_{i}; reduce number of libraries or relax constraints.')

  path_ls = path_ls[1:]

  primer_set = []

  for primer_ls in path_ls:
    if not primer_ls:
      continue
    try:
      primers = primer_df.loc[primer_ls]
      primer_set.append(primers)
    except KeyError:
      print("WARNING: Some primers in the path were not found in primer_df.")
      continue

  if not primer_set:
    print("No valid primer paths were found.")
    return path_ls, float("inf")

  primer_set = pd.concat(primer_set)
  total_cost = primer_set["efficiency"].sum()

  return path_ls, total_cost