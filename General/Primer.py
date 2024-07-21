import itertools as it
from General.args import ARGS


class Primer:
    def __init__(self, primer_df, start, stop, is_r=False):
        assert start < stop
        self.start = start
        self.stop = stop
        self.is_r = is_r  # forward or reverse
        self.l = stop - start  # length
        self.primer_df = primer_df
        self.w = primer_df.at[self.tup(), 'cost']  # total cost value; the fancy notation is b/c
        # of the hierarchal lookup system in panda.df

    def __str__(self):
        return ' '.join(map(str, (self.start, self.stop, self.is_r)))

    def __repr__(self):
        return f'{("r" if self.is_r else "f")}({self.start},{self.stop})'

    def tup(self):
        return (self.start, self.stop, ("r" if self.is_r else "f"))


def actions(primer_df, primer,args):  # returning possible counterparts (forward -> reverse; reverse -> forward)
    # i.e. this method gets the "neighbors"
    if not primer.is_r:  # fwd
        for oligo_l, primer_l in it.product(reversed(range(args.oligo_lmin, args.oligo_lmax + 1)),
                                            range(args.primer_lmin, args.primer_lmax + 1)):

            stop = primer.start + oligo_l
            start = stop - primer_l

            if args.apply_threshold:
                # finds tm of forward and reverse primers
                tm_f = primer_df.at[primer.tup(), 'tm']
                tm_r = primer_df.at[(start, stop, "r"), 'tm']

                # if tm difference is larger then max_difference threshold do not add primer to graph
                if abs(tm_f - tm_r) > args.max_difference:
                    continue

            yield Primer(primer_df, start, stop, is_r=True)

    elif primer.is_r:  # rev
        for overlap_l, primer_l in it.product(reversed(range(args.overlap_lmin, args.overlap_lmax + 1)),
                                              range(args.primer_lmin, args.primer_lmax + 1)):

            start = primer.stop - overlap_l
            stop = start + primer_l

            # filter
            no_split = (primer.start - stop) >= primer.start % 3
            if (stop > primer.start) or (not no_split):
                continue

            yield Primer(primer_df, start, stop)
