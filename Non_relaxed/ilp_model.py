import time
import tracemalloc
from General.get_model import *
from Non_relaxed.single_forbidden_pairs import *
from Non_relaxed.multi_forbidden_pairs import *

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

    # Create ILP model
    model = get_model()

    # Dictionary to store variables for each protein
    protein_vars = {}

    # Add single path and intra-protein forbidden pair ILP constraints
    for name in protein_names:
        graph = graphs[name]
        # find intra-sequence forbidden pairs in each protein
        forbidden_pairs = single_forbidden[name]
        # Returns group of protein graph variables with single-protein constraints
        vars = single_constraints(graph, forbidden_pairs, model)
        # save protein vars in dictionary
        protein_vars[name] = vars

    # Add inter-protein forbidden pair constraints
    multiple_constraints(model, protein_vars, multiple_forbidden)

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

    for name, vars in protein_vars.items():
        actual_values = post_processing(vars)
        protein_paths[name] = actual_values

    # Print results
    for name, vals in protein_paths.items():
        print(f"Protein #{name} ({len(vals)})")
        print(vals)
        print()

    return model.numVars, model.numConstrs, setup_time, setup_memory, ILP_time, ILP_memory, protein_paths, model.objVal
