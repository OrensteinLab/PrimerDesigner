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
All results are written to the specified output directory.  
The output includes a summary CSV containing:

- total runtime  
- memory usage statistics  
- the cost (efficiency) of the selected primer set  
