from General.primer_graphs import *
from General.primer_data import *
from General.constants import *
import networkx as nx
import primer3 as p3
import time
import pandas as pd

mutreg_start=len(upstream_nt)


def run_extension(sequences_nt, mutreg_regions,protein_names, args):
    start_time = time.time()
    greedy_solution, greedy_obj, cross_cnt, protein_cnt, re_iterations = run_greedy(sequences_nt, mutreg_regions,protein_names,args)
    greedy_time = time.time() - start_time

    run_data = []

    run_data.append({"Greedy Solution": greedy_solution,
                     "Greedy Objective": greedy_obj,
                     "Greedy Time": greedy_time,
                     "Cross hybridizations": cross_cnt,
                     "re-iterations": re_iterations,
                     "protein count": protein_cnt})

    run_df = pd.DataFrame(run_data)

    # Write the DataFrame to a CSV file
    run_df.to_csv(f'{args.output}.csv', index=False)



def run_greedy(sequences_nt, mutreg_regions,protein_names, args):

    path_ls={}
    selected_primers = []
    total_cost = 0
    proteins_cnt =0
    cross_cnt = 0
    re_iterations=0

    for i,(seq_nt, mutreg_nt,protein) in enumerate(zip(sequences_nt, mutreg_regions,protein_names)):

        print("Running protein number:",i)

        primer_f, primer_df = create_primer_df(seq_nt,args)

        G = create_graph(primer_df, primer_f,len(mutreg_nt),args)

        cross= True

        iterations = 0

        while cross:

            if iterations>0:
               re_iterations+=1

            iterations+=1

            nodes_to_remove=[]

            # Find the shortest path, excluding the start ('s') and end ('d') nodes
            shortest_path = nx.algorithms.shortest_path(G, 's', 'd', weight='weight')[1:-1]

            cross = False

            for primer in shortest_path:

                start, end, fr = primer

                # Calculate the primer sequence
                primer_seq = seq_nt[start+mutreg_start:end+mutreg_start]  # Assuming start and end are relative to the sequence

                for other_primer in selected_primers:

                    tm = p3.bindings.calc_heterodimer(other_primer, primer_seq, mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6,
                                                      dna_conc=50.0, temp_c=37.0, max_loop=30).tm
                    if tm >= 45:
                        cross_cnt+=1
                        cross = True
                        nodes_to_remove.append(primer)
                        break

            if cross:
                G.remove_nodes_from(nodes_to_remove)
            else:
              # add selected primers to forbidden primer list
              selected_primers.extend([seq_nt[start+mutreg_start:end+mutreg_start] for start,end,fr in shortest_path])
              path_ls[protein]=shortest_path # save shortest path
              # Calculate the cost for the found path
              primer_set = primer_df.loc[[p for p in shortest_path]].copy().reset_index()
              primer_cost = primer_set['cost'].sum()
              total_cost += primer_cost


        # if there was more than 1 iteration it means that there was cross-hybridization in this protein
        if iterations>1:
          proteins_cnt+=1


    return path_ls, total_cost, cross_cnt, proteins_cnt,re_iterations