from General.primer_graphs import *
from PD_mul_var.greedy import *
from PD_mul_var.ilp_model import *
import time
import tracemalloc
import pandas as pd


def run_relaxed_ilp(sequence_nt,mutreg_nt,protein_name,args):
    run_data = []

    primer_f, primer_df = create_primer_df(sequence_nt, args)


    # Creating the Graph
    start_time = time.time()
    tracemalloc.start()
    graph = create_graph(primer_df, primer_f, len(mutreg_nt), args)
    graph_time = int(time.time() - start_time)
    graph_memory = tracemalloc.get_traced_memory()[1] / 10 ** 6  # MB
    tracemalloc.stop()

    # Greedy Solution
    start_time = time.time()
    greedy_solution, greedy_obj = run_greedy(graph, primer_df, args)
    greedy_time = time.time() - start_time

    numVars, numConstrs, setup_time, setup_memory, ILP_time, ILP_memory, actual_values, objective = ilp_model(graph,
                                                                                                              sequence_nt,
                                                                                                              mutreg_nt, args)
    run_data.append({
                    "protein name":protein_name,
                    "Nodes": len(graph.nodes),
                     "Edges": len(graph.edges),
                     "Time (Graph)": graph_time,
                     "MP (Graph)": graph_memory,
                     "Vars": numVars,
                     "Constr": numConstrs,
                     "Time (Setup)": setup_time,
                     "MP (Setup)": setup_memory,
                     "Time (ILP)": ILP_time,
                     "MP (ILP)": ILP_memory,
                     "ILP Solution": actual_values,
                     "ILP Objective": objective,
                     "Greedy Solution": greedy_solution,
                     "Greedy Objective": greedy_obj,
                     "Greedy Time": greedy_time})
    # Assuming all_data is a list of dictionaries
    df = pd.DataFrame(run_data)

    # Write the DataFrame to a CSV file
    df.to_csv(f'{args.output}.csv', index=False)
