import primer3 as p3
import multiprocessing
from General.utils import * 
from itertools import combinations




class PairFinder:

    def __init__(self, sequence_nt1, sequence_nt2=None, upstream_nt_length=len(upstream_nt)):
        self.sequence_nt1 = sequence_nt1
        self.sequence_nt2 = sequence_nt2 if sequence_nt2 else sequence_nt1
        self.memo = {}
        self.upstream_len = upstream_nt_length

    def forbidden_pairs(self, primer, start, sequence_nt, args):

        def search(start, end):
            results = []
            if (end - start) <= 2 * args.primer_lmax:
                tm = calc_tm(sequence_nt[start:end], primer)
                return self.handle_special_case(start, end, sequence_nt, primer, args) if tm >= MAX_TM else []

            mid = (start + end) // 2
            left_end = mid + args.primer_lmax - 1
            right_start = max(mid - (args.primer_lmax - 1), 0)

            if calc_tm(sequence_nt[start:left_end], primer) >= MAX_TM:
                results.extend(search(start, left_end))
            if calc_tm(sequence_nt[right_start:end], primer) >= MAX_TM:
                results.extend(search(right_start, end))

            return results

        return search(max(start - args.allowed_overlap, 0), len(sequence_nt))

    def handle_special_case(self, start, end, sequence_nt, primer, args):
        if (start, end) in self.memo:
            return self.memo[(start, end)]

        if end - start == args.primer_lmax:
            tm = calc_tm(sequence_nt[start:end], primer)
            result = [(start, end)] if tm >= MAX_TM else []
            self.memo[(start, end)] = result
            return result

        shorten_left = self.handle_special_case(start + 1, end, sequence_nt, primer, args)
        shorten_right = self.handle_special_case(start, end - 1, sequence_nt, primer, args)
        result = list(set(shorten_left + shorten_right))
        self.memo[(start, end)] = result
        return result

    def find_all_pairs(self, args):
        forbidden_pairs = set()

        # Determine whether we're checking within a single sequence or between two distinct sequences
        single_sequence = self.sequence_nt1 == self.sequence_nt2

        for p_start in range(len(self.sequence_nt1) - args.primer_lmax + 1):
            p_end = p_start + args.primer_lmax
            primer = self.sequence_nt1[p_start:p_end]

            if single_sequence:
                # Search within the same sequence for potential cross-hybridizations
                off_targets = self.forbidden_pairs(primer, p_end, self.sequence_nt1, args)
            else:
                # Search the second sequence for cross-hybridizations with the primer from the first
                off_targets = self.forbidden_pairs(primer, 0, self.sequence_nt2, args)

            forbidden_pairs.update([((p_start, p_end), target) for target in off_targets])

        return self.sub_pairs_parallel(forbidden_pairs, args)

    def process_pair(self, pair, args):
        sub_forbidden_pairs = []
        (p1_start, p1_end), (p2_start, p2_end) = pair
        primer1 = self.sequence_nt1[p1_start:p1_end]
        primer2 = self.sequence_nt2[p2_start:p2_end]

        for len1 in range(args.primer_lmin, args.primer_lmax + 1):
            for len2 in range(args.primer_lmin, args.primer_lmax + 1):
                for s1 in range(args.primer_lmax - len1 + 1):
                    for s2 in range(args.primer_lmax - len2 + 1):
                        sub1 = primer1[s1:s1 + len1]
                        sub2 = primer2[s2:s2 + len2]
                        if calc_tm(sub1, sub2) >= MAX_TM:
                            sub_p1 = (p1_start + s1 - self.upstream_len, p1_start + s1 + len1 - self.upstream_len)
                            sub_p2 = (p2_start + s2 - self.upstream_len, p2_start + s2 + len2 - self.upstream_len)
                            sub_forbidden_pairs.append((sub_p1, sub_p2))
        return sub_forbidden_pairs

    def sub_pairs_parallel(self, forbidden_pairs, args):
        tasks = [(pair, args) for pair in forbidden_pairs]
        with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
            result = pool.starmap(self.process_pair, tasks)
        return set(item for sublist in result for item in sublist)

