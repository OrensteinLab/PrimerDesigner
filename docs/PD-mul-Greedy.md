
# PD-mul-Greedy  
**Greedy PrimerDesigner version for multiple non-homologous proteins**

## Overview
`PD-mul-Greedy` computes primer sets for **multiple non-homologous proteins** using an iterative greedy strategy.  
For each protein, the longest-path algorithm is applied, and if the slected primer cross-hybridize with any previously selected primers, they are removed from the graph and the process repeats.

## Input Format (Required)
Create a text file where each line contains a protein name and its DNA coding sequence, separated by a tab:

```text
SHP2   ATGACATCGCGGAGATGGTTTCACCCAAAT...
CXAR   ATGGCGCTCCTGCTGTGCTTCGTGCTCC...
```

## Parameters

All parameters match the global program configuration described in the main README and are supplied as command-line arguments.

## Output

Results are written to the specified output directory and include the following files:

### 1. Summary file

**PD_mul_Greedy_summary.csv**

This file contains summary statistics for the run, including:

- total runtime (seconds)  
- peak memory usage  
- number of selected primers per protein  
- Total efficiency of the selected primer sets  

---

### 2. Primer selection file

**PD_mul_Greedy_selected_primers.csv**

This file contains the selected primers for each protein, including:

- protein name  
- primer ID  
- start and end positions  
- strand (`f` for forward, `r` for reverse)  
- primer sequence (5′→3′)  
- primer length  

Each row corresponds to a single primer in the selected primer set.