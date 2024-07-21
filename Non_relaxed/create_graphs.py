from General.primer_graphs import *
import time
import tracemalloc

def create_graphs(mutreg_regions, sequences_nt, protein_names,args):
    """
    Create directed graphs for each protein sequence and measure the time and memory usage.

    Args:
        mutreg_regions (list): List of mutation regions for each protein.
        sequences_nt (list): List of nucleotide sequences corresponding to the protein names.
        protein_names (list): List of protein names.

    Returns:
        tuple: Contains the graphs, graph creation time, graph memory usage, and primer dataframes.
    """

    # Measure the start time and memory usage
    start_time = time.time()
    tracemalloc.start()

    graphs = {}
    primer_dfs = {}

    for i, (mutreg_nt, sequence, protein_name) in enumerate(zip(mutreg_regions, sequences_nt, protein_names)):

        print("Creating Graph for protein: ",i)

        # Create primer dataframe for the given sequence
        primer_f, primer_df = create_primer_df(sequence,args)

        primer_dfs[protein_name] = primer_df

        mutreg_l = len(mutreg_nt)

        graphs[protein_name] = create_graph(primer_df,primer_f,mutreg_l,args)

    # Measure the graph creation time and memory usage
    graph_time = int(time.time() - start_time)
    graph_memory = tracemalloc.get_traced_memory()[1] / 10**6  # MB
    tracemalloc.stop()

    return graphs, graph_time, graph_memory, primer_dfs
