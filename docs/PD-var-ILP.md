
# PD-var-ILP  
**ILP-based PrimerDesigner version for multiple variants of the same protein**

## Overview
`PD-var-ILP` designs an efficient primer set for **multiple variants of the same protein-coding sequence**.  
This version uses an ILP formulation with forbidden-pair constraints to prevent selecting primers that overlap beyond a user-defined threshold across variants.

## Input Format (Required)
Create a text file containing **one line** with the protein name and the **consensus DNA coding sequence** representing all variants (typically obtained after aligning the variants):

```text
SHP2_reference   ATGACATCGCGGAGATGGTTTCACCCAAAT...
```
## Parameters
All parameters follow the global program settings described in the main README and are provided through command-line arguments.  
Specific parameters for PD-var-ILP include:

- **num_proteins**  
  The number of variants for which you want to design specific primers for.  
  *Default:* 3

- **merge_bins**  
  Boolean flag for merging bins corresponding to identical non-overlapping sequences.  
  *Default:* False

- **allowed_overlap**  
  Maximum allowed overlap between primer pairs from overlapping regions.  
  *Default:* 6

## Output

Results are written to the specified output directory and include the following files:

### 1. Summary file

**var_ILP_results.csv**

This file contains summary statistics for the run, including:

- total runtime and runtime of major program stages  
- peak memory usage  
- total number of selected primers by greedy and ILP algorithms
- ILP and greedy objective values (efficiencies) 

---
### 2. Primer selection file

**PD_var_ILP_selected_primers_greedy.csv & PD_var_ILP_selected_primers_ILP.csv**

These file contains the optimal primer-selection paths for each protein variant for the freeedy and ILP approaches.

For each variant, the file records:

- ordered list of selected primers  
- start and end positions relative to the reference sequence  
- strand (`f` for forward, `r` for reverse)  
- primer sequence (5′→3′)  
- primer length  

