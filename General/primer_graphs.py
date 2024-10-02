from General.Primer import *
from General.primer_data import *
import networkx as nx


def create_graph(primer_df, primer_f, mutreg_l,args):
    print("Creating graph")
    # initialize graph
    G = nx.DiGraph()

    ## all forward primers that end before mutreg and pass threshold
    primers_init = [Primer(primer_df, p.start, p.stop) for _, p in primer_f.query('stop<=0').iterrows() if
                    check_threshold(p.tm, p.gc,args)]
    for primer in primers_init:
        G.add_edge('s', primer.tup(), weight=primer.w)  # i ntializing the s-primer connection
        dfs(G, primer,primer_df,mutreg_l,args)  # create the rest of the graph

    print(G)

    return G


def dfs(G, primer, primer_df, mutreg_l,args):  # CREATING the graph using DFS

    tm = primer_df.at[primer.tup(), 'tm']
    gc = primer_df.at[primer.tup(), 'gc']

    # base case (end)
    if (primer.start >= mutreg_l) and primer.is_r and check_threshold(tm, gc,args):
        G.add_edge(primer.tup(), 'd', weight=0.)  # G is global variable defined in next section
        return

    for next_primer in actions(primer_df, primer,args):

        is_new = not G.has_node(next_primer.tup())  # check if primer node exists already

        tm = primer_df.at[next_primer.tup(), 'tm']
        gc = primer_df.at[next_primer.tup(), 'gc']

        passes_threshold = check_threshold(tm, gc,args)  # check if primer passes threshold

        # only add edge if primer passes threshold
        if passes_threshold:
            G.add_edge(primer.tup(), next_primer.tup(), weight=next_primer.w)  # weight is the cost of the new primer
        if passes_threshold and is_new:
            dfs(G, next_primer, primer_df, mutreg_l,args)  # recursive call to dfs


def paths_ct(G, u, d):  # count total # of paths between two points

    if u == d:
        return 1
    else:
        if not G.nodes[u]:  # npaths attribute is the # of paths out of u
            G.nodes[u]['npaths'] = sum(paths_ct(G, c, d) for c in G.successors(u))
        return G.nodes[u]['npaths']
