import itertools as it
import networkx as nx
import pandas as pd

def run_greedy(graphs, primer_dfs,multiple_forbidden,protein_names):

  path_ls = {}

  for protein in protein_names:
    nodes_to_remove = []
    G = graphs[protein]
    G_sub = G.copy()
    for i,path in enumerate(path_ls.values()):
      other_protein = protein_names[i]
      for p,n in it.product(path,G_sub.nodes()):
        if n=='s' or n=='d':
          continue
        pn_forbidden = ((p[0],p[1]),(n[0],n[1])) in multiple_forbidden[(other_protein,protein)]
        if pn_forbidden:
          nodes_to_remove.append(n)
    G_sub.remove_nodes_from(set(nodes_to_remove))
    try:
      path_ls[protein]= [primer for primer in nx.algorithms.shortest_path(G_sub,'s','d', weight='weight')][1:-1]
    except:
      print(f'WARNING: No feasible primer sequence for lib_{i}; reduce number of libraries or relax constraints.')

  primer_set = pd.DataFrame()
  for i,primer_ls in enumerate(path_ls.values()):
    protein = protein_names[i]
    primer_df = primer_dfs[protein]
    primer_set1 = primer_df.loc[primer_ls].copy().reset_index()
    primer_set1['lib_i'] = protein
    primer_set = pd.concat([primer_set, primer_set1])
  primer_set['tile_i'] = primer_set.groupby(['lib_i','fr']).cumcount()
  primer_set = primer_set[['lib_i','tile_i']+primer_set.columns[:-2].to_list()]

  #print(primer_set.groupby('lib_i').cost.sum())
  cost = primer_set.cost.sum()

  return path_ls, cost
