# PD-mul-ILP  
**ILP-based PrimerDesigner version for multiple non-homologous proteins**

## Overview
`PD-mul-ILP` designs an efficient primer set for **multiple non-homologous protein-coding sequences**.  
This version first identifies all primer pairs that may cross-hybridize between different proteins using a brute-force search, and then encodes these as forbidden-pair constraints in an ILP formulation. The ILP ensures that no two selected primers form a high-risk cross-hybridizing pair, while maximizing the overall primer efficiency.

## Input Format (Required)
Create a text file where **each line** contains a protein name and its DNA coding sequence, separated by a tab character:

```text
SHP2   ATGACATCGCGGAGATGGTTTCACCCAAATATCACTGGTGTGGAGGCAGAAAACCTACTGTTGACAAGAGGAGT...
CXAR   ATGGCGCTCCTGCTGTGCTTCGTGCTCCTGTGCGGAGTAGTGGATTTCGCCAGAAGTTTGAGTATCACTACTCC...
```

## Parameters

All parameters match the global program configuration described in the main README and are supplied as command-line arguments.

## Output

Results are written to the specified output directory and include the following files:

### 1. Summary file

**mul_ilp_results.csv**

This file contains summary statistics for the run, including:

- total runtime and runtime of individual program stages (graph construction, cross-hybridization detection, and ILP solving)  
- number of detected intra-protein and inter-protein cross-hybridization risks  
- peak memory usage  
- number of selected primers per protein  
- total predicted efficiency of the selected primer sets  

---

### 2. Primer selection file

**mul_ilp_selected_primers.csv**

This file contains the optimal primer-selection paths for each protein.

For each protein, the file records:

- ordered list of selected primers  
- start and end positions  
- strand (`f` for forward, `r` for reverse)  
- primer sequence (5′→3′)  
- primer length  
