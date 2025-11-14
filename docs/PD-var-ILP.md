
# PD-var-ILP  
**ILP-based PrimerDesigner version for multiple variants of the same protein**

## Overview
`PD-var-ILP` designs an efficient primer set for **multiple variants of the same protein-coding sequence**.  
This version uses an ILP formulation with forbidden-pair constraints to prevent selecting primers that overlap beyond a user-defined threshold across variants.

## Input Format (Required)
Create a text file with one line contains the protein name and a reference DNA coding sequence of all variants (can be computed after alignment)

```text
SHP2_reference   ATGACATCGCGGAGATGGTTTCACCCAAAT...
```
## Parameters
All parameters follow the global program settings described in the main README and are provided through command-line arguments.  
Specific parameters for PD-var-ILP include

- **num_proteins**  
The number of variants for which you want to design specific primers for

- **merge_bins**  
  Boolean flag for merging bins corresponding to identical non-overlapping sequences.  

## Output
All results are saved to the specified output directory.  
The output includes:

- total runtime  
- memory usage statistics  
- ILP objective value - the primer efficiency

A JSON file containing the optimal primer-selection paths for each variant is also produced.
