# PD-single-LPath  
**Longest-Path PrimerDesigner version for a single protein**

## Overview
`PD-single-LPath` computes the most efficient primer set for **one protein-coding sequence** by constructing the primer graph and selecting the longest path.  
This avoids ILP and provides a fast and deterministic solution.

## Input Format (Required)
Create a text file with a single line containing the protein name and its DNA coding sequence, separated by a tab:

```text
SHP2  ATGACATCGCGGAGATGGTTTCACCCAAATATCACTGGTGTGGAGGCAGAAAACCTACTGTTGACAAGAGGAGT..
```

## Parameters
The parameters are identical to the global program parameters described in the main README and are provided through command-line arguments.

## Output

Results are written to the specified output directory and include the following files:

### 1. Summary file

**PD_single_LPath_results.csv**

This file contains summary statistics for the run, including:

- total runtime  
- peak memory usage  
- number of selected primer pairs  
- total predicted efficiency of the selected primer set  

---

### 2. Primer selection file

**PD_single_LPath_selected_primers.csv**

This file contains the optimal primer-selection path for the input protein.

Each row corresponds to a selected primer and includes:

- protein name  
- primer index in the assembly path  
- start and end positions relative to the coding sequence  
- strand (`f` for forward, `r` for reverse)  
- primer sequence (5′→3′)  
- primer length  
