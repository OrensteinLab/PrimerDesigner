from PD_mul_greedy.run_mul_greedy import *
from General.args import *
import sys

sys.argv = [
        sys.argv[0],
        "--file_path", "data/100_ccds_protein_sequences.txt",
        "--output", "Experiment_results/mul_greedy_ccds",
        ]

args = get_args()

# Create output directory if not exists
output_dir = Path(args.output)
output_dir.mkdir(parents=True, exist_ok=True)
print(f"[INFO] Output directory: {output_dir.resolve()}")

# ============================================================
# LOAD SEQUENCES
# ============================================================
print(f"[INFO] Reading protein coding sequences from: {args.file_path}")
all_mutreg_regions, all_full_sequences, all_protein_names = read_sequences(args.file_path)
print(f"[INFO] Total proteins loaded: {len(all_protein_names)}")

# all_mutreg_regions, all_full_sequences, all_protein_names = all_mutreg_regions[:2], all_full_sequences[:2], all_protein_names[:2] 

run_mul_greedy(all_full_sequences,all_mutreg_regions,all_protein_names,args)