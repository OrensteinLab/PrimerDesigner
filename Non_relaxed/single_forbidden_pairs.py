from Non_relaxed.PairFinder import *
import gurobipy as gp

def single_pairs(protein_names, sequences_nt, args):
    """
    Identify intra-protein forbidden pairs for each protein sequence.

    Args:
        protein_names (list): List of protein names.
        sequences_nt (list): List of nucleotide sequences corresponding to the protein names.

    Returns:
        dict: A dictionary where keys are protein names and values are lists of forbidden pairs.
    """
    single_forbidden = {}

    for name, sequence in zip(protein_names, sequences_nt):
        # Find forbidden pairs in a single protein sequence
        pairs_finder = PairFinder(sequence)
        forbidden_pairs = pairs_finder.find_all_pairs(args)
        single_forbidden[name] = forbidden_pairs

    return single_forbidden

def single_constraints(graph, forbidden_pairs, model):
    """
    Add constraints to the optimization model based on the graph structure and intra-protein forbidden primer pairs.

    Args:
        graph (networkx.Graph): The graph representing connections between nodes.
        forbidden_pairs (list): List of forbidden pairs of primers.
        model (gurobipy.Model): The optimization model.

    Returns:
        gurobipy.Var: Variables representing the protein in the model.
    """
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
    for n in graph_nodes + ['s', 'd']:  # Adding 's' & 'd' back just here
        v = str(n)
        model.addConstr(
            sum(protein_vars[i, j] for i, j in ij.select(v, '*')) -
            sum(protein_vars[j, i] for j, i in ij.select('*', v)) ==
            (1 if v == 's' else -1 if v == 'd' else 0), v
        )

    # Adding intra-protein forbidden pairs constraints
    cnt=0
    for cnt, pairs in enumerate(forbidden_pairs):
        all_edges = []
        for pair in pairs:
            node1 = (pair[0], pair[1], 'f')
            node2 = (pair[0], pair[1], 'r')
            all_edges.append(protein_vars.sum(str(node1), '*'))
            all_edges.append(protein_vars.sum(str(node2), '*'))
        model.addConstr(gp.quicksum(all_edges) <= 1)

    print("Number of intra-protein forbidden pairs constraint: ",cnt)


    return protein_vars
