from itertools import combinations
from Non_relaxed.PairFinder import *
import gurobipy as gp

def multiple_pairs(protein_names, sequences_nt,args):
    """
    Identify inter-protein forbidden pairs for combinations of protein sequences.

    Args:
        protein_names (list): List of protein names.
        sequences_nt (list): List of nucleotide sequences corresponding to the protein names.

    Returns:
        dict: A dictionary where keys are tuples of protein name pairs and values are lists of forbidden pairs.
    """
    multiple_forbidden = {}

    # Generate all possible pairs of proteins
    protein_pairs = list(combinations(range(len(protein_names)), 2))

    for p1, p2 in protein_pairs:
        sequence1 = sequences_nt[p1]
        sequence2 = sequences_nt[p2]

        protein1 = protein_names[p1]
        protein2 = protein_names[p2]

        # Find forbidden pairs between two protein sequences
        pairs_finder = PairFinder(sequence1, sequence2)
        forbidden_pairs = pairs_finder.find_all_pairs(args)

        multiple_forbidden[(protein1, protein2)] = forbidden_pairs

    return multiple_forbidden

def multiple_constraints(model, model_vars, multiple_forbidden):
    """
    Add constraints to the optimization model based on inter-protein forbidden primer pairs between each pair of protein seqeunces.

    Args:
        model (gurobipy.Model): The optimization model.
        model_vars (dict): Dictionary of model variables indexed by protein names.
        multiple_forbidden (dict): Dictionary of forbidden pairs for protein pairs.

    Returns:
        None
    """

    forbidden_pairs_cnt=0

    for protein_pair, forbidden_pairs in multiple_forbidden.items():
        protein1, protein2 = protein_pair

        protein1_vars = model_vars[protein1]
        protein2_vars = model_vars[protein2]

        # Add <= 1 constraint for every forbidden pair
        for forbidden_pair in forbidden_pairs:
            all_edges = []
            pair1, pair2 = forbidden_pair

            # Collect outgoing edges for the forbidden pairs
            node1_fwd = (pair1[0], pair1[1], 'f')
            node1_rev = (pair1[0], pair1[1], 'r')
            all_edges.append(protein1_vars.sum(str(node1_fwd), '*'))
            all_edges.append(protein1_vars.sum(str(node1_rev), '*'))

            node2_fwd = (pair2[0], pair2[1], 'f')
            node2_rev = (pair2[0], pair2[1], 'r')
            all_edges.append(protein2_vars.sum(str(node2_fwd), '*'))
            all_edges.append(protein2_vars.sum(str(node2_rev), '*'))

            # Constraint: the sum of outgoing edges of forbidden pairs <= 1
            model.addConstr(gp.quicksum(all_edges) <= 1)
            forbidden_pairs_cnt+=1

    print("Number of inter-protein forbidden pairs constraints: ",forbidden_pairs_cnt)


