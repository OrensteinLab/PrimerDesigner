import json
import gurobipy as gp
import pandas as pd
import primer3 as p3
from Bio.Seq import Seq
import networkx as nx

# defines global cross-hybridization thershold
MAX_TM = 45.0
PCR = p3.thermoanalysis.ThermoAnalysis()

# Define global upstream and downstream regions
upstream_nt = 'GCTAGTGGTGCTAGCCCCGCGAAATTAATACGACTCACTATAGGGTCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACAT'
downstream_nt = 'GGAGGGTCTGGGGGAGGAGGCAGTGGCATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCG'

def calc_tm(seq1, seq2):
    return PCR.calc_heterodimer(seq1, seq2,).tm

def get_model(file_path='gurobi.json'):

    with open(file_path, 'r') as json_file:
        params = json.load(json_file)

    env = gp.Env(params=params)

    # Create the model within the Gurobi environment
    model = gp.Model('max-sum', env=env)

    model.ModelSense = -1  # This  makes model maximize

    return model

def read_sequences(file_path):

    mutreg_regions = []
    protein_names=[]

    # read  protein coding-sequences from the file path
    with open(file_path) as file:
        for line in file.readlines():
            p_name, mutreg_region = line.strip().split('\t')
            mutreg_regions.append(mutreg_region)
            protein_names.append(p_name)

    full_sequences = []

    # add constant upstream and downstream regions to each sequence
    for mutreg_nt in mutreg_regions:

        sequence = upstream_nt + mutreg_nt + downstream_nt
        full_sequences.append(sequence)

    return mutreg_regions,full_sequences,protein_names

def revcomp(seq):
  return str(Seq(seq).reverse_complement())

def subsequences(sequence,primer_lmin,primer_lmax): #Generates all subsequences w/ all poss. start-stop pairs
  ls = []
  for j in range(primer_lmin, primer_lmax+1): #length
    for i in range(len(sequence)-j+1): #starting index
      start = i
      stop = i+j
      ls.append([sequence[start:stop], start, stop, stop-start])
  return ls


class NoPathError(Exception):
    """Raised when no path exists between source and target in the DAG."""
    pass

def longest_path_dag(G, source, target):  
    """
    Finds the maximum-weight path between source and target in a DAG.
    Raises NoPathError if no path exists.
    """
    # Step 1: Negate weights to convert max → min
    G_neg = G.copy()
    for u, v, data in G_neg.edges(data=True):
        data['weight'] = -data.get('weight', 1)
    
    # Step 2: Topological order
    topo_order = list(nx.topological_sort(G_neg))
    
    # Step 3: Initialize distances and parent pointers
    dist = {v: float('inf') for v in G_neg.nodes()}
    dist[source] = 0
    parent = {v: None for v in G_neg.nodes()}
    
    # Step 4: Relax edges in topological order
    for u in topo_order:
        if dist[u] != float('inf'):
            for v in G_neg.successors(u):
                w = G_neg[u][v]['weight']
                if dist[v] > dist[u] + w:
                    dist[v] = dist[u] + w
                    parent[v] = u
    
    # Step 5: Raise error if target not reachable
    if dist[target] == float('inf'):
        raise NoPathError(f"No path found from {source!r} to {target!r} in the DAG.")
    
    # Step 6: Reconstruct path
    path = []
    curr = target
    while curr is not None:
        path.append(curr)
        curr = parent[curr]
    path.reverse()
    
    return path
