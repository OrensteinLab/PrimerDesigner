import time
import tracemalloc
import gurobipy as gp
from General.utils import *
from collections import defaultdict
from dataclasses import dataclass
from typing import List, Any

@dataclass
class ILPResult:
    num_vars: int
    num_constraints: int
    setup_time: float
    setup_peak_mem_mb: float
    optimize_time: float
    optimize_peak_mem_mb: float
    paths: List[List[Any]]   # list of paths, one per protein; each path is a list of nodes
    objective: float
    status: int    

def ilp_model(graph, sequence_nt, mutreg_nt, args):
    
    print("Running ILP")
    setup_start = time.time()
    tracemalloc.start()
    # Create Gurobi Model
    model = get_model()

    # Getting Nodes/Edges
    graph_edges = graph.edges(data=True)
    graph_nodes = [node for node in graph.nodes if node != 's' and node != 'd']  # removing s & d nodes

    mutreg_start = len(UPSTREAM_NT)
    l_range = (args.allowed_overlap + 1, args.primer_lmax + 1)

    def create_bins(l_range):
        all_bins = defaultdict(list)

        for node in graph_nodes:
            node_start, node_end, direction = node
            for bin_start in range(node_start, node_end):
                for bin_length in range(*l_range):
                    bin_end = bin_start + bin_length
                    if bin_end > node_end:
                        break
                    all_bins[(bin_start, bin_end, direction)].append(node)

        return dict(all_bins) 

    all_bins = create_bins(l_range)
    print("Number of Constraints:", len(all_bins))
    if len(all_bins):
        print("Average Vars Per Constraint", 1 / len(all_bins) * sum(len(val) for _, val in all_bins.items()))

    def unite_bins(all_bins):
        """Merge all bins corresponding to identical sequences (plus direction)."""
        united_bins = {}

        for (start, end, fr) in all_bins.keys():
            seq_len = end - start
            bin_seq = sequence_nt[start + mutreg_start:end + mutreg_start]

            # Skip if already processed
            if (bin_seq, fr) in united_bins:
                continue

            union_nodes = []

            # Efficiently find all occurrences of bin_seq in sequence_nt
            search_start = 0
            while True:
                match_pos = sequence_nt.find(bin_seq, search_start)
                if match_pos == -1:
                    break
                bin_start = match_pos - mutreg_start
                key = (bin_start, bin_start + seq_len, fr)
                if key in all_bins:
                    union_nodes.extend(all_bins[key])
                search_start = match_pos + 1  # allow overlapping matches

            united_bins[(bin_seq, fr)] = union_nodes

        return united_bins


    # if the merge_bins flag is set to true, calls function to unite bins with identical sequences
    if args.merge_bins:
        united_bins = unite_bins(all_bins)
        print("Number of united bins constraints:", len(united_bins))
        if len(all_bins):
            print("Average Vars Per Constraint", 1 / len(united_bins) * sum(len(ls) for seq, ls in united_bins.items()))
        all_bins = united_bins

    # Converting Graphs to lists of model variables
    ij = gp.tuplelist()
    w_ij = gp.tupledict()

    # add dinary variable for each graph edge
    for edge in graph_edges:
        l = (str(edge[0]), str(edge[1]))
        ij.append(l)
        w_ij[l] = edge[-1]['weight']

    print("Graph edges: ",len(graph_edges))

    # Adding Variables to model
    model_vars = model.addVars(ij, obj=w_ij, vtype=gp.GRB.BINARY)

    print("Finished ILP Variable Creations")

    for cnt, nodes in enumerate(all_bins.values()):
        all_edges = []
        if cnt % (len(all_bins) // 25) == 0:
            print(int(cnt / len(all_bins) * 100))
        for node in nodes:
            all_edges.append(model_vars.sum(str(node), '*'))
        model.addConstr(gp.quicksum(all_edges) <= 1)
        
    print("Finished Intersection constraints!")

    # implement single Path Constraints
    model.addConstrs(
    (
        gp.quicksum(model_vars[i, j] for i, j in ij.select(v, '*')) -
        gp.quicksum(model_vars[j, i] for j, i in ij.select('*', v)) ==
        (args.num_proteins if v == 's' else -args.num_proteins if v == 'd' else 0)
        for v in map(str, graph_nodes + ['s', 'd'])
    ),
    name="single_path"
    )

    setup_time = time.time() - setup_start
    setup_memory = tracemalloc.get_traced_memory()[1] / 10 ** 6
    tracemalloc.stop()

    tracemalloc.start()
    start_time = time.time()
    model.optimize()
    ILP_time = time.time() - start_time
    ILP_memory = tracemalloc.get_traced_memory()[1] / 10 ** 6
    tracemalloc.stop()

    def post_processing(model_vars):
        # retrieve solution from graph
        all_proteins = [['s'] for _ in range(args.num_proteins)]
        true_edges = [index for index, var in model_vars.items() if var.X != 0]

        while true_edges:
            edge = true_edges.pop(0)
            added_edge = False
            for protein_list in all_proteins:
                if edge[0] == protein_list[-1]:
                    protein_list.append(edge[1])
                    added_edge = True
                    break
            if not added_edge:
                true_edges.append(edge)

        return [p[1:-1] for p in all_proteins if len(p) > 2]


    actual_values= []
    if model.Status == gp.GRB.OPTIMAL:
        objective = model.objVal
        actual_values = post_processing(model_vars)
        for cnt, vals in enumerate(actual_values):
            print(f"Protein #{cnt + 1} ({len(vals)})")
            print(vals)
            print()
    else:
        print(f"Model status: {model.Status}")
        objective = float('nan') 

    return ILPResult(
        num_vars=model.numVars,
        num_constraints=model.numConstrs,
        setup_time=setup_time,
        setup_peak_mem_mb=setup_memory,
        optimize_time=ILP_time,
        optimize_peak_mem_mb=ILP_memory,
        paths=actual_values,
        objective=objective,
        status=model.Status
    )
