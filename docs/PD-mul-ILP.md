# PD-mul-ILP  
**ILP-based PrimerDesigner version for multiple non-homologous proteins**

## Overview
`PD-mul-ILP` designs an efficient primer set for **multiple non-homologous protein-coding sequences**.  
This version first identifies primer pairs that may cross-hybridize (intra- and inter-protein) using a brute-force search, then encodes these as forbidden-pair constraints in an ILP. 

Requires a Gurobi license (`gurobi.json` in the repository root). See the main README.

## Input Format (Required)
Create a text file where **each line** contains a protein name and its DNA coding sequence, separated by a **tab**:

```text
SHP2	ATGACATCGCGGAGATGGTTTCACCCAAATATCACTGGTGTGGAGGCAGAAAACCTACTGTTGACAAGAGGAGT...
CXAR	ATGGCGCTCCTGCTGTGCTTCGTGCTCCTGTGCGGAGTAGTGGATTTCGCCAGAAGTTTGAGTATCACTACTCC...
```

Upstream and downstream flanks come from `--config` (see the main README). 

## Parameters

All parameters match the global program configuration described in the main README and are supplied as command-line arguments. Cross-hybridization uses `max_tm` from the config file (heterodimer cutoff), not `--min_tm` / `--max_tm`. `--allowed_overlap` is used when enumerating forbidden pairs.

## Output

Results are written to the directory given by `--output`:

### 1. Summary file

**`mul_ilp_results.csv`**

One row per run (appended if the file already exists):

| Column | Meaning |
|---|---|
| `num_proteins` | Number of input sequences |
| `graph_time_sec` | Graph construction time |
| `ilp_num_vars`, `ilp_num_constraints` | ILP model size |
| `ilp_intra_forbidden_cnt`, `ilp_inter_forbidden_cnt` | Forbidden-pair counts |
| `forbidden_time_sec` | Cross-hybridization search time |
| `ilp_setup_time_sec`, `ilp_optimize_time_sec` | ILP setup and solve times |
| `ilp_feasibility` | `FEASIBLE` if Gurobi status is 2, else `INFEASIBLE` |
| `total_primer_efficiency` | ILP objective (sum of selected primer efficiencies) |
| `num_primers` | Total primers across all proteins |

Peak memory and per-protein primer counts are not written to this CSV.

### 2. Primer selection file

**`mul_ilp_selected_primers.csv`**

One row per selected primer (not a single cell listing the whole path):

| Column | Meaning |
|---|---|
| `protein_name` | Protein |
| `primer_index` | Order along that protein’s path (0-based) |
| `start`, `end` | Coordinates relative to the mutagenized CDS start (upstream flank is negative) |
| `orientation` | `f` or `r` |
| `length` | `end - start` |
| `seq` | Primer sequence (5′→3′) |
