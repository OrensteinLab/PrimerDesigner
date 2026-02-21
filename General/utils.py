import json
import gurobipy as gp
import pandas as pd
import primer3 as p3
from Bio.Seq import Seq
import networkx as nx
import random 
from Bio import SeqIO
from dataclasses import dataclass

@dataclass(frozen=True)
class Config:
    upstream: str
    downstream: str
    max_tm: float = 45.0

def load_config(path="config.json") -> Config:

    with open(path) as f:
        cfg = json.load(f)

    return Config(
        upstream=cfg["upstream_nt"],
        downstream=cfg["downstream_nt"],
        max_tm=float(cfg.get("max_tm", 45.0)),
    )

_PCR = None

def get_PCR():
    global _PCR
    if _PCR is None:
        _PCR = p3.thermoanalysis.ThermoAnalysis()
    return _PCR

def calc_tm(seq1, seq2):
    return get_PCR().calc_heterodimer(seq1, seq2).tm

def read_fasta(fasta_path):

    with open(fasta_path, "r") as handle:
        record = next(SeqIO.parse(handle, "fasta"))
        return str(record.seq)
    
def get_model(file_path='gurobi.json'):

    with open(file_path, 'r') as json_file:
        params = json.load(json_file)

    env = gp.Env(params=params)

    model = gp.Model('max-sum', env=env)
    model.ModelSense = -1  # maximize

    return model

def read_sequences(file_path, cfg):

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

        sequence = cfg.upstream + mutreg_nt + cfg.downstream
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

def sample_paths_dag_uniform(G, source, target, k=100, max_tries=10000, seed=42, verbose = False):
    rng = random.Random(seed)

    topo = list(nx.topological_sort(G))

    # DP: number of paths from node -> target
    ways = {n: 0 for n in G.nodes()}
    ways[target] = 1
    for u in reversed(topo):
        if u == target:
            continue
        s = 0
        for v in G.successors(u):
            s += ways.get(v, 0)
        ways[u] = s

    if ways[source] == 0:
        raise ValueError(f"source {source!r} cannot reach target {target!r}")

    paths, seen = [], set()
    tries = 0

    while len(paths) < k and tries < max_tries:
        if verbose and tries % 100 == 0:
            print("Number of paths found:", len(paths), "tries:", tries, end="\r")
        
        tries += 1
        curr = source
        path = [curr]

        while curr != target:
            candidates = [(v, ways[v]) for v in G.successors(curr) if ways[v] > 0]
            if not candidates:
                break
            vs, ws = zip(*candidates)
            # weighted random choice
            total = sum(ws)
            r = rng.uniform(0, total)
            acc = 0.0
            nxt = None
            for v, w in candidates:
                acc += w
                if r <= acc:
                    nxt = v
                    break
            curr = nxt
            path.append(curr)

        if path[-1] == target:
            t = tuple(path)
            if t not in seen:
                seen.add(t)
                paths.append(path)

    return paths 