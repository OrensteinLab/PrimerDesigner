import time
import tracemalloc
from General.utils import *
from Non_relaxed.PairFinder import *
import gurobipy as gp


def single_constraints(graph,sequence, model,args):
    """
    Add constraints to the optimization model based on the graph structure and intra-protein forbidden primer pairs.

    Args:
        graph (networkx.Graph): The graph representing connections between nodes.
        forbidden_pairs (list): List of forbidden pairs of primers.
        model (gurobipy.Model): The optimization model.

    Returns:
        gurobipy.Var: Variables representing the protein in the model.
    """
    pairs_finder = PairFinder(sequence)
    forbidden_pairs = pairs_finder.find_all_pairs(args)

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
            protein_vars.sum(str((pair[0], pair[1], 'f')), '*'),
            protein_vars.sum(str((pair[0], pair[1], 'r')), '*')
        ]) <= 1
        for pairs in forbidden_pairs
        for pair in pairs
    ),
    name="intra_forbidden")

    num_constr = sum(len(pairs) for pairs in forbidden_pairs)

    return protein_vars, num_constr


def multiple_constraints(model, model_vars, protein_names, sequences_nt, args):
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

    # Generate all possible pairs of proteins
    protein_pairs = list(combinations(range(len(protein_names)), 2))

    multiple_forbidden = {}

    for p1, p2 in protein_pairs:

        sequence1 = sequences_nt[p1]
        sequence2 = sequences_nt[p2]

        protein1 = protein_names[p1]
        protein2 = protein_names[p2]

        # Find forbidden pairs between two protein sequences
        pairs_finder = PairFinder(sequence1, sequence2)
        forbidden_pairs = pairs_finder.find_all_pairs(args)

        multiple_forbidden[(protein1,protein2)]= forbidden_pairs

        protein1_vars = model_vars[protein1]
        protein2_vars = model_vars[protein2]

        num_constr = sum(len(pairs) for pairs in forbidden_pairs)

        forbidden_pairs_cnt += num_constr

        # Add <= 1 constraint for every forbidden pair
        model.addConstrs(
        (gp.quicksum([
            protein1_vars.sum(str((pair1[0], pair1[1], 'f')), '*'),
            protein1_vars.sum(str((pair1[0], pair1[1], 'r')), '*'),
            protein2_vars.sum(str((pair2[0], pair2[1], 'f')), '*'),
            protein2_vars.sum(str((pair2[0], pair2[1], 'r')), '*')
        ]) <= 1
        for (pair1, pair2) in forbidden_pairs),
        name="inter_protein_forbidden"
        )

    print("Number of inter-protein forbidden pairs constraints: ",forbidden_pairs_cnt)

    return multiple_forbidden, forbidden_pairs_cnt


def run_ilp(sequences_nt, protein_names, graphs, args):
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
    print("Running ILP")

    # Measure setup time and memory usage
    setup_start = time.time()
    tracemalloc.start()

    # Create ILP model
    model = get_model()

    # Dictionary to store variables for each protein
    protein_vars = {}

    single_pair_cnt = 0 
    # Add single path and intra-protein forbidden pair ILP constraints
    for name,sequence in zip(protein_names,sequences_nt):

        graph = graphs[name]
        # Returns group of protein graph variables with single-protein constraints
        vars,num_constraints = single_constraints(graph,sequence,model,args)

        single_pair_cnt += num_constraints
        # save protein vars in dictionary
        protein_vars[name] = vars

    # Add inter-protein forbidden pair constraints
    multiple_forbidden, multi_pairs_cnt = multiple_constraints(model, protein_vars, protein_names, sequences_nt, args)

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

    def post_processing(vars):
        """
        Extracts the path from the optimized ILP variables.

        Args:
            vars (dict): Dictionary of ILP variables.

        Returns:
            list: List representing the path.
        """
        path = ['s']
        true_edges = [index for index, var in vars.items() if var.X != 0]

        while true_edges:
            edge = true_edges.pop(0)
            added_edge = False
            if edge[0] == path[-1]:
                path.append(edge[1])
                added_edge = True
            if not added_edge:
                true_edges.append(edge)

        return path

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

    return multiple_forbidden, model.numVars, model.numConstrs, single_pair_cnt,multi_pairs_cnt, setup_time, setup_memory, ILP_time, ILP_memory, protein_paths, objective
