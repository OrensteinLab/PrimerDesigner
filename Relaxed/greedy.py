import pandas as pd
import itertools as it
import networkx as nx


def run_greedy(G, primer_df,args):

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

  # # https://pandas.pydata.org/pandas-docs/stable/user_guide/merging.html
  primer_set = pd.DataFrame()
  for i,primer_ls in enumerate(path_ls):
    primer_set1 = primer_df.loc[primer_ls].copy().reset_index()
    primer_set1['lib_i'] = i
    primer_set = pd.concat([primer_set, primer_set1])
  primer_set['tile_i'] = primer_set.groupby(['lib_i','fr']).cumcount()
  primer_set = primer_set[['lib_i','tile_i']+primer_set.columns[:-2].to_list()]

  #print(primer_set.groupby('lib_i').cost.sum())
  cost = primer_set.cost.sum()

  return path_ls, cost