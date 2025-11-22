import time
import tracemalloc
from General.utils import *
import gurobipy as gp

from dataclasses import dataclass
from typing import Dict, Any
import gurobipy as gp
import time, tracemalloc
from PD_mul_ILP.brute_force import *
from itertools import combinations

@dataclass
class ILPResult:
    num_vars: int
    num_constraints: int
    setup_time: float
    setup_memory: float
    ilp_time: float
    ilp_memory: float
    protein_paths: Dict[str, list]
    objective: float
    status: int

def find_forbidden_pairs(protein_names, sequences_nt, args):

    multiple_forbidden_cnt=0
    single_forbidden_cnt=0

    single_forbidden = {}

    print("Finding intra-protein forbidden pairs constraints")

    total_single = len(protein_names)

    for idx, (protein, sequence) in enumerate(zip(protein_names, sequences_nt), start=1):
        print(f"[{idx}/{total_single}] Processing {protein}...")

        forbidden_pairs = find_forbidden_pairs_intra(sequence,args)

        single_forbidden_cnt += len(forbidden_pairs)

        single_forbidden[protein] = forbidden_pairs

    print("Number of intra-protein forbidden pairs constraints: ",single_forbidden_cnt)

    # Generate all possible pairs of proteins
    protein_pairs = list(combinations(range(len(protein_names)), 2))
    total_pairs = len(protein_pairs)
    multiple_forbidden = {}

    print(f"Finding inter-protein forbidden pairs for {total_pairs} combinations...")

    for idx, (p1, p2) in enumerate(protein_pairs, start=1):

        print(f"[{idx}/{total_pairs}] Fractions of pairs considered")

        sequence1 = sequences_nt[p1]
        sequence2 = sequences_nt[p2]

        protein1 = protein_names[p1]
        protein2 = protein_names[p2]

        # Find forbidden pairs between two protein sequences
        forbidden_pairs = find_forbidden_pairs_inter(sequence1,sequence2,args)

        multiple_forbidden[(protein1,protein2)]= forbidden_pairs

        multiple_forbidden_cnt += len(forbidden_pairs)

    print("Number of inter-protein forbidden pairs constraints: ",multiple_forbidden_cnt)

    return single_forbidden,multiple_forbidden, single_forbidden_cnt,multiple_forbidden_cnt


def add_constraints(graphs, protein_names, single_forbidden, multiple_forbidden, model):

    model_vars = {}

    # create variables and add intra-forbidden pairs constraints within same protein
    for protein in protein_names:

        single_forbidden_pairs = single_forbidden[protein]

        graph = graphs[protein]

        # Extracting edges and nodes from the graph
        graph_edges = graph.edges(data=True)
        graph_nodes = [node for node in graph.nodes if node not in ('s', 'd')]  # Removing 's' & 'd' nodes

        # Converting graph edges to model variables
        ij = gp.tuplelist()
        w_ij = gp.tupledict()

        for edge in graph_edges:
            l = (str(edge[0]), str(edge[1]))  # i, j
            ij.append(l)
            w_ij[l] = edge[-1]['weight']

        # Adding variables to the model
        protein_vars = model.addVars(ij, obj=w_ij, vtype=gp.GRB.BINARY)

        # Adding  graph single path constraints
        model.addConstrs(
        (
            gp.quicksum(protein_vars[i, j] for i, j in ij.select(v, '*')) -
            gp.quicksum(protein_vars[j, i] for j, i in ij.select('*', v)) ==
            (1 if v == 's' else -1 if v == 'd' else 0)
            for v in map(str, graph_nodes + ['s', 'd'])
        ),
        name="flow_conservation"
        )

        model.addConstrs(
        (
            gp.quicksum([
                protein_vars.sum(str((p1[0], p1[1], 'f')), '*'),
                protein_vars.sum(str((p1[0], p1[1], 'r')), '*'),
                protein_vars.sum(str((p2[0], p2[1], 'f')), '*'),
                protein_vars.sum(str((p2[0], p2[1], 'r')), '*'),
            ]) <= 1
            for (p1, p2) in single_forbidden_pairs
        ),
        name="intra_forbidden"
        )

        model_vars[protein] = protein_vars

    # add inter-forbidden pairs constraints between multiple difference proteins
    # Generate all possible pairs of proteins
    protein_pairs = list(combinations(range(len(protein_names)), 2))

    for p1, p2 in protein_pairs:

        protein1 = protein_names[p1]
        protein2 = protein_names[p2]

        protein1_vars = model_vars[protein1]
        protein2_vars = model_vars[protein2]

        multi_forbidden_pairs = multiple_forbidden[(protein1,protein2)]

        # Add <= 1 constraint for every forbidden pair
        model.addConstrs(
        (gp.quicksum([
            protein1_vars.sum(str((p1[0], p1[1], 'f')), '*'),
            protein1_vars.sum(str((p1[0], p1[1], 'r')), '*'),
            protein2_vars.sum(str((p2[0], p2[1], 'f')), '*'),
            protein2_vars.sum(str((p2[0], p2[1], 'r')), '*')
        ]) <= 1
        for (p1, p2) in multi_forbidden_pairs),
        name="inter_protein_forbidden"
        )


    return model_vars


def run_ilp(single_forbidden, multiple_forbidden, protein_names, graphs):
    """
    Run the ILP (Integer Linear Programming) model to optimize primer design with intra- and inter-protein constraints.

    Args:
        single_forbidden (dict): Dictionary of forbidden pairs within each protein.
        multiple_forbidden (dict): Dictionary of forbidden pairs between different proteins.
        protein_names (list): List of protein names.
        graphs (dict): Dictionary of graphs representing the connections between nodes for each protein.

    Returns:
        tuple: Contains model statistics and paths for each protein.
    """

    # Measure setup time and memory usage
    setup_start = time.time()
    tracemalloc.start()

    print("Creating ILP model...")

    # Create ILP model
    model = get_model()

    print("Adding constraints to ILP model...")

    protein_vars = add_constraints(graphs, protein_names, single_forbidden, multiple_forbidden, model)

    # Capture setup time and memory usage
    setup_time = time.time() - setup_start
    setup_memory = tracemalloc.get_traced_memory()[1] / 10**6
    tracemalloc.stop()

    # Measure optimization time and memory usage
    tracemalloc.start()
    start_time = time.time()

    model.optimize()

    ILP_time = time.time() - start_time
    ILP_memory = tracemalloc.get_traced_memory()[1] / 10**6
    tracemalloc.stop()

    def post_processing(vars_tupledict):
        # collect chosen edges
        chosen = [(i, j) for (i, j), var in vars_tupledict.items() if var.X > 0.5]
        succ = {}
        for i, j in chosen:
            succ[i] = j
        path = ['s']
        while path[-1] in succ:
            path.append(succ[path[-1]])
            if len(path) > len(chosen) + 2:  # safety
                break
        return path[1:-1]

    # Post-process to extract paths for each protein
    protein_paths = {}

    if model.Status == gp.GRB.OPTIMAL:
        objective = model.objVal
        for name, vars in protein_vars.items():
            actual_values = post_processing(vars)
            protein_paths[name] = actual_values

        # Print results
        for name, vals in protein_paths.items():
            print(f"Protein #{name} ({len(vals)})")
            print(vals)
            print()
    else: 
        print(f"Model status: {model.Status}")
        objective = float('nan')


    return ILPResult(
    model.numVars,
    model.numConstrs,
    setup_time,
    setup_memory,
    ILP_time,
    ILP_memory,
    protein_paths,
    objective,
    model.Status
    )
