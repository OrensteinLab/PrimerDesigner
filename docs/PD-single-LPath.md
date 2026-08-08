# PD-single-LPath  
**Longest-Path PrimerDesigner version for a single protein**

## Overview
`PD-single-LPath` computes the most efficient primer set for **one protein-coding sequence** by constructing the primer graph and selecting the longest path.  
This avoids ILP and provides a fast and deterministic solution. No Gurobi license is required.

## Input Format (Required)
Create a text file with a **single line**: protein name and DNA coding sequence, separated by a **tab**. If the file has extra lines, only the first sequence is used.

```text
SHP2	ATGACATCGCGGAGATGGTTTCACCCAAATATCACTGGTGTGGAGGCAGAAAACCTACTGTTGACAAGAGGAGT...
```

Upstream and downstream flanks come from `--config` (see the main README). Example input: `data/SPAP_reference.txt`.

## Parameters
The parameters are identical to the global program parameters described in the main README and are provided through command-line arguments.

## Output

Results are written to the directory given by `--output`:

### 1. Summary file

**`PD_single_LPath_results.csv`**

One row with:

| Column | Meaning |
|---|---|
| `protein_name` | Input sequence name |
| `graph_nodes`, `graph_edges` | Primer-graph size |
| `graph_time_sec` | Graph construction time |
| `longest_path_efficiency` | Sum of selected primer efficiencies |
| `total_time_sec` | End-to-end runtime |
| `num_primers` | Number of selected primers (not primer pairs) |

Peak memory is not written.

### 2. Primer selection file

**`PD_single_LPath_selected_primers.csv`**

One row per selected primer:

| Column | Meaning |
|---|---|
| `protein_name` | Input sequence name |
| `primer_index` | Order along the assembly path (0-based) |
| `start`, `stop` | Coordinates relative to the mutagenized CDS start (upstream flank is negative) |
| `orientation` | `f` (forward) or `r` (reverse) |
| `seq` | Primer sequence (5′→3′) |
| `efficiency` | Predicted primer efficiency |

There is no primer-length column.
