from Non_relaxed.create_graphs import *
from Non_relaxed.single_forbidden_pairs import *
from Non_relaxed.multi_forbidden_pairs import *
from Non_relaxed.greedy import *
from Non_relaxed.ilp_model import *
import pandas as pd

def run_non_relaxed(mutreg_regions,sequences_nt,protein_names,args):

    # Initialize a list to store all data
    run_data = []

    # Create graphs and measure time and memory usage
    graphs, graph_time, graph_memory, primer_dfs = create_graphs(mutreg_regions, sequences_nt, protein_names,args)

    # Measure time for calculating intra-protein forbidden primer pairs
    start_time = time.time()
    single_forbidden = single_pairs(protein_names, sequences_nt,args)
    single_pairs_time = time.time() - start_time

    # Measure time for calculating inter-protein forbidden primer pairs
    start_time = time.time()
    multiple_forbidden = multiple_pairs(protein_names, sequences_nt,args)
    multi_pairs_time = time.time() - start_time

    # Run the greedy solution and measure time
    start_time = time.time()
    greedy_solution, greedy_obj = run_greedy(graphs, primer_dfs, multiple_forbidden, protein_names)
    greedy_time = time.time() - start_time

    # Run the ILP model and collect statistics
    numVars, numConstrs, setup_time, setup_memory, ILP_time, ILP_memory, protein_paths, objective = run_ilp(single_forbidden, multiple_forbidden, protein_names, graphs)

    # Append collected data to the all_data list
    run_data.append({
    "num proteins": len(protein_names),
    "Time (Graph)": graph_time,
    "MP (Graph)": graph_memory,
    "Vars": numVars,
    "Constr": numConstrs,
    "Time (Setup)": setup_time,
    "MP (Setup)": setup_memory,
    "multi pair time": multi_pairs_time,
    "single pair time": single_pairs_time,
    "Time (ILP)": ILP_time,
    "MP (ILP)": ILP_memory,
    "ILP Solution": protein_paths,
    "ILP Objective": objective,
    "Greedy Solution": greedy_solution,
    "Greedy Objective": greedy_obj,
    "Greedy Time": greedy_time
    })


    # Convert the data to a DataFrame and save it as a CSV file
    run_df = pd.DataFrame(run_data)
    run_df.to_csv(f'{args.output}.csv', index=False)
