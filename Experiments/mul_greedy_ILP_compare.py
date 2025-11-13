import pandas as pd
from pathlib import Path
from PD_mul_greedy.run_mul_greedy import *
from General.args import *
import sys

def main():

    sys.argv = [
        sys.argv[0],
        "--file_path", "data/10_protein_coding_sequences_example.txt",
        "--output", "Experiment_results/mul_greedy_compare",
        ]

    args = get_args()


    # Create output directory if not exists
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"[INFO] Output directory: {output_dir.resolve()}")

    print(f"[INFO] Reading protein coding sequences from: {args.file_path}")
    all_mutreg_regions, all_full_sequences, all_protein_names = read_sequences(args.file_path)
    print(f"[INFO] Total proteins loaded: {len(all_protein_names)}")


    summary_rows = []  # collect rows for a final aggregated CSV

    for i in range(2, len(all_protein_names)):
        print(f"\n[INFO] Processing {i} protein(s)...")

        # Use the first i proteins (make sure to pass the SLICED lists!)
        mutreg_regions = all_mutreg_regions[:i]
        sequences_nt   = all_full_sequences[:i]
        protein_names  = all_protein_names[:i]

        # run and collect the summary row
        summary_row = run_mul_greedy(sequences_nt, mutreg_regions, protein_names, args)
        summary_rows.append(summary_row)

    # write the aggregated results once
    greedy_results_csv = output_dir / "greedy_results.csv"
    pd.DataFrame(summary_rows).to_csv(greedy_results_csv, index=False)
    print(f"[INFO] Saved aggregate greedy results to: {greedy_results_csv}")


if __name__ == '__main__':
    main()