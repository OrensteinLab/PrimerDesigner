# Reproducing PrimerDesigner Experiments

## Reproducing comparison on SpAP and CVB3 proteins

From the main directory, run:

```bash
python -m Comparisons.compare_CVB3
python -m Comparisons.compare_SpAP
```

This will reproduce the comparisons between PrimerDesigner and competing methods on the SpAP and CVB3 genes.

The scripts generate the following output files in `Comparisons/Results/`:

- `*_comparison.csv` — summary of runtime, memory usage, and average PCR efficiency for each method  
- `*_primers.csv` — selected primer sets for each method  
- `null_paths_primers_*.csv` — 1000 randomly sampled primer paths from the primer graph (used for null distribution analysis)
