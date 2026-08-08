# PD-var-ILP  
**ILP-based PrimerDesigner version for multiple variants of the same protein**

## Overview
`PD-var-ILP` designs primer sets for **multiple variants of the same protein-coding sequence**.  
It uses an ILP with forbidden-pair constraints so primers from different variants do not overlap beyond a user-defined threshold. A greedy baseline is run on the same graph and written alongside the ILP solution.

Requires a Gurobi license (`gurobi.json` in the repository root). See the main README.

## Input Format (Required)
Create a text file with **one line**: protein name and the **consensus DNA coding sequence** (typically after aligning the variants), separated by a **tab**. If the file has extra lines, only the first sequence is used.

```text
SHP2_reference	ATGACATCGCGGAGATGGTTTCACCCAAAT...
```

Upstream and downstream flanks come from `--config` (see the main README).

## Parameters
All parameters follow the global program settings described in the main README. PD-var-ILP-specific flags:

- **`--num_proteins`**  
  Number of variants to design primers for.  
  *Default:* 3

- **`--merge_bins`**  
  Merge bins that correspond to identical non-overlapping sequences.  
  *Default:* off 

- **`--allowed_overlap`**  
  Maximum allowed overlap between primer pairs from overlapping regions.  
  *Default:* 6

## Output

Results are written to the directory given by `--output`:

### 1. Summary file

**`var_ILP_results.csv`**

One row with:

| Column | Meaning |
|---|---|
| `protein_name` | Input sequence name |
| `num_variants` | Value of `--num_proteins` |
| `seq_length` | Mutagenized CDS length |
| `graph_nodes`, `graph_edges` | Primer-graph size |
| `graph_time_sec` | Graph construction time |
| `ilp_num_vars`, `ilp_num_constraints` | ILP model size |
| `ilp_setup_time_sec`, `ilp_optimize_time_sec` | ILP setup and solve times |
| `ilp_path_length`, `greedy_path_length` | Total primers selected by ILP / greedy |
| `ILP_primer_efficiency`, `greedy_primer_efficiency` | Objective values |
| `ILP_solution_feasibility`, `greedy_solution_feasibility` | `FEASIBLE` or `INFEASIBLE` |
| `greedy_time_sec` | Greedy runtime |

Peak memory is not written.

### 2. Primer selection files

**`var_ILP_selected_primers_ILP.csv`**  
**`var_ILP_selected_primers_greedy.csv`**

One row per selected primer:

| Column | Meaning |
|---|---|
| `protein_name` | Input sequence name |
| `variant_index` | Variant (0-based; one path per variant) |
| `primer_index` | Order along that variant’s path (0-based) |
| `start`, `stop` | Coordinates relative to the consensus CDS start (upstream flank is negative) |
| `orientation` | `f` or `r` |
| `seq` | Primer sequence (5′→3′) |
| `efficiency` | Predicted primer efficiency |

There is no primer-length column.
