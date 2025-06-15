from General.primer_graphs import *
from General.primer_data import *
from General.utils import *
import networkx as nx
import primer3 as p3
import time
import pandas as pd

mutreg_start= len(upstream_nt)

def run_pd_mul_nh(sequences_nt, mutreg_regions,protein_names, args):
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

        if i % 10 == 0:
            print(f"Running protein {i}/{len(sequences_nt)}")

        primer_df = create_primer_df(seq_nt,args)

        G = create_graph(primer_df, len(mutreg_nt),args)

        cross= True

        iterations = 0

        while cross:

            if iterations>0:
               re_iterations+=1

            iterations+=1

            nodes_to_remove=[]

            # Find the shortest path, excluding the start ('s') and end ('d') nodes
            try:
                shortest_path = nx.algorithms.shortest_path(G, 's', 'd', weight='weight')[1:-1]
            except nx.NetworkXNoPath:
                print(f"No valid path found for {protein}. Skipping.")
                break

            cross = False

            primer_seqs = {
                (start, end, fr): seq_nt[start + mutreg_start:end + mutreg_start]
                for (start, end, fr) in shortest_path
            }

            for primer, primer_seq in primer_seqs.items():
                for other_seq in selected_primers:
                    tm = calc_tm(other_seq, primer_seq)

                    if tm >= MAX_TM:
                        cross_cnt+=1
                        cross = True
                        nodes_to_remove.append(primer)
                        break

            if cross:
                G.remove_nodes_from(nodes_to_remove)
            else:
              # add selected primers to forbidden primer list
              selected_primers.extend(primer_seqs.values())

              path_ls[protein]=shortest_path # save shortest path
              # Calculate the cost for the found path and add to total cost
              primer_cost = primer_df.loc[shortest_path, 'cost'].sum()
              total_cost += primer_cost


        # if there was more than 1 iteration it means that there was cross-hybridization in this protein
        if iterations>1:
          proteins_cnt+=1


    return path_ls, total_cost, cross_cnt, proteins_cnt,re_iterations