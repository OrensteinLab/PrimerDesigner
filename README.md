
# PrimerDesigner

## Introduction

This repository contains the source code for PrimerDesigner, a tool designed to find the most efficient primer set with complete coverage and no cross hybridization risk for protein synthesis using assembly PCR. It accompanies the paper titled "PrimerDesigner: Designing efficient primers for synthesizing large protein libraries without cross-hybridization".

![alt text](https://github.com/OrensteinLab/PrimerDesigner/blob/main/primer_design_illustration.png)


## Creating the Conda Environment

This project uses **Conda** for environment and dependency management.  
You can easily recreate the exact environment used in this project using the provided `environment.yml` file.

- **1. Create the environment**

  ```bash
  conda env create -f environment.yml

- **2. Activate the environment**

   ```bash
  conda activate primer_env

## The Different Versions

**PrimerDesigner** has 4 different version options:

- **[docs/PD-single-LPath.md](docs/PD-single-LPath.md)**
  - Finds the most efficient primer set for a single protein by computing the longest path in the primer graph.

- **PD-var-ILP**
  - Finds the most efficient primer set  for a number of variants of the same protein-coding sequence by imposing ILP constraints to prevent the selections of primers from overlapping subsequences.

- **PD-mul-Greedy**
  - Finds the optimal primer set for multiple non-homologous proteins by applying an iterative greedy apporoach using the longest path algorithm in the primer graph.

- **PD-mul-ILP**
  - Finds the optimal primer set for multiple non-homologous proteins.
  - Identifies all cross-hybridize primers pairs using a brute-force algorithm and adds them as forbidden-pair ILP constraints.


## Setting Up the Tool

First, create a file named `gurobi.json` containing the details for the Gurobi license and put it in the main directory.

```json
{
  "WLSACCESSID": "XXXXX",
  "WLSSECRET": "XXXXX",
  "LICENSEID": 12345
}
```

Create a text file containing the protein names and their DNA coding-sequences. Each line should contain a protein's name and its DNA coding-sequence, separated by a tab. For example: 

```text
SHP2  ATGACATCGCGGAGATGGTTTCACCCAAATATCACTGGTGTGGAGGCAGAAAACCTACTGTTGACAAGAGGAGT....
CXAR  ATGGCGCTCCTGCTGTGCTTCGTGCTCCTGTGCGGAGTAGTGGATTTCGCCAGAAGTTTGAGTATCACTACTCC....
```

## Running PrimerDesginer

To execute PrimerDesginer, use the following command:

```bash
python ./tool.py --version <version> --file_path <file-path>  --output <output-file>
```
- **file_path**: The file path of the protein coding-sequences
- **version**: Specifies which version of the algorithm to run. The options are: `PD-single-LPath`, `PD-mul-ILP`, `PD-mul-Greedy` and `PD-var-ILP` (default:  `PD-single-LPath`)
- **output**: The path of the folder that the run's output file will be saved to.
  
The other arguments are optional and include the algorithm parameters:

### Parameters

- **primer_lmin**, **primer_lmax**  
  Minimum and maximum primer lengths.  
  *Default:* 18, 30  

- **oligo_lmin**, **oligo_lmax**  
  Minimum and maximum oligonucleotide lengths.  
  *Default:* 195, 205  

- **overlap_lmin**, **overlap_lmax**  
  Minimum and maximum overlap length between oligonucleotides.  
  *Default:* 45, 50  

- **allowed_overlap**  
  Allowed overlap between primer pairs.  
  *Default:* 6  

- **num_proteins**  
  Number of variants of the same sequence.  
  *Used for:* `PD-var-ILP` version only.  
  *Default:* 3  

- **apply_threshold**  
  Boolean flag for applying primer quality threshold.  
  *Default:* False  

- **min_gc**, **max_gc**  
  Minimum and maximum GC content threshold (in %).  
  *Default:* 40, 60  

- **min_tm**, **max_tm**  
  Minimum and maximum melting temperature (Tm) thresholds (in °C).  
  *Default:* 58, 65  

- **max_difference**  
  Maximum allowed Tm difference between forward and reverse primers in a pair.  
  *Default:* 3  

- **merge_bins**  
  Boolean flag for merging bins corresponding to identical non-overlapping sequences.  
  *Used for:* `PD-var-ILP` version.  
  *Default:* False

Example command:
```bash
python ./tool.py --version PD-mul-ILP --file_path example_proteins.txt  --output run_output --primer_lmin 20 --primer_lmax 26 --oligo_lmin 180 --oligo_lmax 200
```










