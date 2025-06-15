from Non_relaxed.create_graphs import *
from Non_relaxed.greedy import *
from Non_relaxed.ilp_model import *


def run_non_relaxed(mutreg_regions,sequences_nt,protein_names,args):

    # Initialize a list to store all data
    run_data = []

    # Create graphs and measure time and memory usage
    graphs, graph_time, graph_memory, primer_dfs = create_graphs(mutreg_regions, sequences_nt, protein_names,args)

    # Run the ILP model and collect statistics
    multiple_forbidden, numVars, numConstrs, single_pair_cnt,multi_pairs_cnt, setup_time, setup_memory, ILP_time, ILP_memory, protein_paths, objective = run_ilp(sequences_nt, protein_names, graphs, args)

    # Run the greedy solution and measure time
    start_time = time.time()
    greedy_solution, greedy_obj = run_greedy(graphs, primer_dfs, multiple_forbidden, protein_names) 
    greedy_time = time.time() - start_time

    # Append collected data to the all_data list
    run_data.append({
    "num proteins": len(protein_names),
    "Time (Graph)": graph_time,
    "MP (Graph)": graph_memory,
    "Vars": numVars,
    "Constr": numConstrs,
    "Single_pair_constr":single_pair_cnt,
    "Multi_pair_constr":multi_pairs_cnt,
    "Time (Setup)": setup_time,
    "MP (Setup)": setup_memory,
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
