from General.Primer import *
from General.primer_data import *
import networkx as nx


def create_graph(primer_df, mutreg_l, args):

    print("Creating graph")

    # initialize graph
    G = nx.DiGraph()

    tm_dict = primer_df['tm'].to_dict()
    gc_dict = primer_df['gc'].to_dict()

    # all forward primers that end before mutreg and pass threshold
    primers_init = []
    for row in primer_df.reset_index().itertuples(index=False):
        if row.stop <= 0 and row.fr == "f" and check_threshold(row.tm, row.gc, args):
            try:
                primers_init.append(Primer(primer_df, row.start, row.stop))
            except Exception as e:
                print(f"Failed to construct primer for {key}: {e}")

    for primer in primers_init:
        key = primer.tup()
        G.add_edge('s', key, weight=primer.w)
        dfs(G, primer, primer_df, mutreg_l, args, tm_dict, gc_dict)

    return G


def dfs(G, primer, primer_df, mutreg_l, args, tm_dict, gc_dict):

    key = primer.tup()

    tm = tm_dict[key]
    gc = gc_dict[key]

    if (primer.start >= mutreg_l) and primer.is_r and check_threshold(tm, gc, args):

        G.add_edge(key, 'd', weight=0.0)
        return

    for next_primer in actions(primer_df, primer, args):

        next_key = next_primer.tup()

        is_new = not G.has_node(next_key)

        tm_next = tm_dict[next_key]
        gc_next = gc_dict[next_key]

        if check_threshold(tm_next, gc_next, args):
            G.add_edge(key, next_key, weight=next_primer.w)
            if is_new:
                dfs(G, next_primer, primer_df, mutreg_l, args, tm_dict, gc_dict)


def check_threshold(tm, gc, args):
    
    if not args.apply_threshold:
        return True

    # Convert thresholds
    gc_min = args.min_gc / 100.0
    gc_max = args.max_gc / 100.0
    tm_min = float(args.min_tm)
    tm_max = float(args.max_tm)

    passed = gc_min <= gc <= gc_max and tm_min <= tm <= tm_max
    
    return passed
