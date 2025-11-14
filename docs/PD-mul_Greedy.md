
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

Results are written to the specified output directory.
The output includes a summary CSV containing:
	•	total runtime
	•	memory usage statistics
	•	efficiency of the primer sets
	•	selected primer set for each protein

A JSON file containing the optimal primer-selection paths for each protein is also produced.
