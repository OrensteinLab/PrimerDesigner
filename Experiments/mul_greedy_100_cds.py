from PD_mul_greedy.run_mul_greedy import *
from General.args import get_args
import General.utils as GU

def main():

    args = get_args()

    args.output = "Results"

    args.file_path = "data/100_ccds_protein_sequences.txt"

    cfg = GU.load_config("configs/SPAP_experiment.json")

    # Create output directory if not exists
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"[INFO] Output directory: {output_dir.resolve()}")

    # ============================================================
    # LOAD SEQUENCES
    # ============================================================
    print(f"[INFO] Reading protein coding sequences from: {args.file_path}")
    all_mutreg_regions, all_full_sequences, all_protein_names = read_sequences(args.file_path,cfg)
    print(f"[INFO] Total proteins loaded: {len(all_protein_names)}")

    run_mul_greedy(all_full_sequences,all_mutreg_regions,all_protein_names,args,cfg)

if __name__ == '__main__':
    main()
