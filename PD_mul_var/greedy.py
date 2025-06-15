import pandas as pd
import itertools as it
import networkx as nx


def run_greedy(G, primer_df,args):

  print("Running greedy algorithm")

  path_ls = [[]]
  G_sub = G.copy()
  for i in range(args.num_proteins):
    nodes_to_remove = []
    for p,n in it.product(path_ls[-1],G_sub.nodes()):
      if n=='s' or n=='d':
        continue
      pn_intersect = n[1]-p[0] > args.allowed_overlap and p[1]-n[0] > args.allowed_overlap and n[2]==p[2]  ## check overlap
      if pn_intersect:
        nodes_to_remove.append(n)
    G_sub.remove_nodes_from(set(nodes_to_remove))
    try:
      path_ls.append([primer for primer in nx.algorithms.shortest_path(G_sub,'s','d', weight='weight')][1:-1])
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
  total_cost = primer_set["cost"].sum()

  return path_ls, total_cost