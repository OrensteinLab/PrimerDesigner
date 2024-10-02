import time
import tracemalloc
from General.get_model import *
import gurobipy as gp
from General.constants import *


def ilp_model(graph, sequence_nt, mutreg_nt, args):
    print("Running ILP")
    setup_start = time.time()
    tracemalloc.start()
    # Create Model
    model = get_model()
    # Getting Nodes/Edges
    graph_edges = graph.edges(data=True)
    graph_nodes = [node for node in graph.nodes if node != 's' and node != 'd']  # removing s & d nodes

    mutreg_start = len(upstream_nt)
    nt_range = (-len(upstream_nt), len(mutreg_nt) + len(downstream_nt) + 1)  # range of nucleotides
    l_range = (args.allowed_overlap + 1, args.primer_lmax + 1)

    def create_bins():  # bins take the form (start, end, [])
        all_bins = {(start, start + len, c): [] for start in range(*nt_range) for len in range(*l_range) for c in
                    ('f', 'r')}

        for node in graph_nodes:  # for each node, add it into all the necessary bins
            for bin_start in range(node[0], node[1]):
                for bin_length in range(*l_range):
                    if bin_start + bin_length > node[1]:
                        break
                    else:
                        all_bins[(bin_start, bin_start + bin_length, node[2])].append(node)

        all_bins = {key: val for key, val in all_bins.items() if val}  # no empty bins

        return all_bins

    all_bins = create_bins()
    print("Number of Constraints:", len(all_bins))
    print("Average Vars Per Constraint", 1 / len(all_bins) * sum(len(val) for _, val in all_bins.items()))

    def unite_bins():
        # this function merges all bins corresponding to identical sequences
        united_bins = {}

        for bin in all_bins.keys():
            start = bin[0]
            end = bin[1]
            fr = bin[2]
            length = end - start
            bin_sequence = sequence_nt[start + mutreg_start:end + mutreg_start]
            if (bin_sequence, fr) in united_bins:
                continue
            bin_union = []
            # find subsequences that are identical to the bin sequence and add them the the union of bins
            for i in range(len(sequence_nt) - length + 1):
                if sequence_nt[i:i + length] == bin_sequence:
                    start = i - mutreg_start  # subtract from start to account for upstream region
                    bin_union.extend(all_bins[(start, start + length, fr)])  # add the bins with identical sequences to the union list

            # add bin union tp dictionary in the sequence entry
            united_bins[(bin_sequence, fr)] = bin_union

        return united_bins

    # if the merge_bins flag is set to true, calls function to unite bins with identical sequences
    if args.merge_bins:
        united_bins = unite_bins()
        print("Number of united bins constraints:", len(united_bins))
        print("Average Vars Per Constraint", 1 / len(united_bins) * sum(len(ls) for seq, ls in united_bins.items()))
        all_bins = united_bins

    # Converting Graphs to lists of model variables
    ij = gp.tuplelist()
    w_ij = gp.tupledict()

    for edge in graph_edges:
        l = (str(edge[0]), str(edge[1]))  # i, j
        ij.append(l)
        w_ij[l] = edge[-1]['weight']

    print("Finished Conversion")

    # Adding Variables to model
    model_vars = model.addVars(ij, obj=w_ij, vtype=gp.GRB.BINARY)
    print("Finished Variable Creations")

    # implements forbidden pairs constraints - bins
    for cnt, nodes in enumerate(all_bins.values()):
        all_edges = []
        if cnt % (len(all_bins) // 25) == 0:
            print(int(cnt / len(all_bins) * 100))
        for node in nodes:
            all_edges.append(model_vars.sum(str(node), '*'))
        model.addConstr(gp.quicksum(all_edges) <= 1)
    print("Finished Intersection constraints!")

    # implement single Path Constraints
    for n in graph_nodes + ['s', 'd']:  # adding s & d back just here
        v = str(n)
        model.addConstr(sum(model_vars[i, j] for i, j in ij.select(v, '*')) - sum(
            model_vars[j, i] for j, i in ij.select('*', v)) == (
                            args.num_proteins if v == 's' else -1 * args.num_proteins if v == 'd' else 0), v)

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

        return all_proteins

    actual_values = post_processing(model_vars)
    for cnt, vals in enumerate(actual_values):
        print(f"Protein #{cnt + 1} ({len(vals)})")
        print(vals)
        print()

    return model.numVars, model.numConstrs, setup_time, setup_memory, ILP_time, ILP_memory, actual_values, model.objVal
