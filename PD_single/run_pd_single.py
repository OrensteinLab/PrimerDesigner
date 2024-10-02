from General.primer_graphs import *
import time
import tracemalloc
import pandas as pd


def run_single(sequence_nt,mutreg_nt,protein_name,args):

    global_start = time.time()

    run_data = []

    primer_f, primer_df = create_primer_df(sequence_nt, args)


    # Creating the Graph
    start_time = time.time()
    tracemalloc.start()
    graph = create_graph(primer_df, primer_f, len(mutreg_nt), args)
    graph_time = int(time.time() - start_time)
    graph_memory = tracemalloc.get_traced_memory()[1] / 10 ** 6  # MB
    tracemalloc.stop()

    shortest_path = nx.algorithms.shortest_path(graph, 's', 'd', weight='weight')[1:-1]

    print("Shortest path: ", shortest_path)

    primer_set = primer_df.loc[[p for p in shortest_path]].copy().reset_index()
    primer_cost = primer_set['cost'].sum()

    total_time = time.time()-global_start


    run_data.append({
                    "protein name":protein_name,
                    "Nodes": len(graph.nodes),
                     "Edges": len(graph.edges),
                     "Time (Graph)": graph_time,
                     "MP (Graph)": graph_memory,
                    "Shortest path": shortest_path,
                    "Cost": primer_cost,
                    "Total time":total_time
    })

    # Assuming all_data is a list of dictionaries
    df = pd.DataFrame(run_data)

    # Write the DataFrame to a CSV file
    df.to_csv(f'{args.output}.csv', index=False)
