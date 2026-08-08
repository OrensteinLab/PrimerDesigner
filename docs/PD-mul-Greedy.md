# PD-mul-Greedy  
**Greedy PrimerDesigner version for multiple non-homologous proteins**

## Overview
`PD-mul-Greedy` computes primer sets for **multiple non-homologous proteins** using an iterative greedy strategy.  
For each protein, the longest-path algorithm is applied, and if the selected primers cross-hybridize with any previously selected primers, they are removed from the graph and the process repeats. 

## Input Format (Required)
Create a text file where **each line** contains a protein name and its DNA coding sequence, separated by a **tab**:

```text
SHP2	ATGACATCGCGGAGATGGTTTCACCCAAAT...
CXAR	ATGGCGCTCCTGCTGTGCTTCGTGCTCC...
```

Upstream and downstream flanks come from `--config` (see the main README). 

## Parameters

All parameters match the global program configuration described in the main README and are supplied as command-line arguments. Cross-hybridization uses `max_tm` from the config file (heterodimer cutoff).

## Output

Results are written to the directory given by `--output`:

### 1. Summary file

**`PD_mul_Greedy_summary.csv`**

One row for the whole run:

| Column | Meaning |
|---|---|
| `num_proteins` | Number of input sequences |
| `total_primer_efficiency` | Sum of efficiencies over all accepted paths |
| `greedy_time_sec` | End-to-end runtime |
| `cross_hybridizations_cnt` | Heterodimer hits that triggered a retry |
| `proteins_with_reiterations_cnt` | Proteins that needed more than one longest-path attempt |
| `total_reiterations` | Extra longest-path attempts across proteins |
| `unresolved_proteins_cnt`, `unresolved_proteins` | Proteins with no valid path |
| `total_primers` | Total primers across all accepted paths |


### 2. Per-protein metrics

**`PD_mul_Greedy_per_protein_metrics.csv`**

One row per protein: `protein_name`, sequence / mutreg lengths, `primer_df_time_sec`, `graph_time_sec`, `iterations`, `resolved`, `num_primers`, `primer_efficiency`.

### 3. Primer selection file

**`PD_mul_Greedy_selected_primers.csv`**

One row per selected primer:

| Column | Meaning |
|---|---|
| `protein_name` | Protein |
| `primer_index` | Order along that protein’s path (0-based) |
| `start`, `stop` | Coordinates relative to the mutagenized CDS start (upstream flank is negative) |
| `orientation` | `f` or `r` |
| `seq` | Primer sequence (5′→3′) |

