import argparse

def get_args():

    parser = argparse.ArgumentParser(description='Process some sequences.')
    parser.add_argument('--version', type=str, default='PD-mul-var', help='Version to use. Options: PD-mul-var, Non_relaxed,PD-mul-nh, PD-single', choices=['PD-mul-var', 'Non_relaxed', 'PD-mul-nh','PD-single'])
    parser.add_argument('--file_path', type=str, required=True, help='Path of the input file')
    parser.add_argument('--output', type=str, required=True, help='Path of the csv program output file')
    parser.add_argument('--primer_lmin', type=int, default=18, help='Minimum primer length')
    parser.add_argument('--primer_lmax', type=int, default=30, help='Maximum primer length')
    parser.add_argument('--overlap_lmin', type=int, default=45, help='Minimum oligo overlap length')
    parser.add_argument('--overlap_lmax', type=int, default=50, help='Maximum oligo overlap length')
    parser.add_argument('--oligo_lmin', type=int, default=195, help='Minimum oligo length')
    parser.add_argument('--oligo_lmax', type=int, default=205, help='Maximum oligo length')
    parser.add_argument('--allowed_overlap', type=int, default=6, help='Maximum allowed overlap between selected primers')
    parser.add_argument('--apply_threshold', type=bool, default=False, help='Choose whether to apply efficiency thresholds')
    parser.add_argument('--min_gc', type=int, default=40, help='Minimum GC content')
    parser.add_argument('--max_gc', type=int, default=60, help='Maximum GC content')
    parser.add_argument('--min_tm', type=int, default=58, help='Minimum tm')
    parser.add_argument('--max_tm', type=int, default=65, help='Maximum tm')
    parser.add_argument('--max_difference', type=int, default=3, help='Maximum tm difference between forward and reverse primer')
    parser.add_argument('--merge_bins', type=bool, default=False,help='Choose whether to merge bins (for the relaxed version)')
    parser.add_argument('--num_proteins', type=int, default=3, help='Number of variants of the same protein (for the relaxed version)')

    args = parser.parse_args()
    return args

ARGS = get_args()