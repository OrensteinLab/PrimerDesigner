
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

**PrimerDesigner** offers four variants, each supporting a different primer-design setting.
Detailed input/output specifications are provided in the linked documentation.

- **PD-single-LPath**  
  Longest-path algorithm for primer design on a single protein.  
  Documentation: [docs/PD-single-LPath.md](docs/PD-single-LPath.md)

- **PD-var-ILP**  
  ILP-based design for variants of the same protein, preventing selection of primers that overlap beyond a specified threshold.  
  Documentation: [docs/PD-var-ILP.md](docs/PD-var-ILP.md)

- **PD-mul-Greedy**  
  Greedy iterative approach for multiple non-homologous proteins, using repeated longest-path computations.  
  Documentation: [docs/PD-mul-Greedy.md](docs/PD-mul-Greedy.md)

- **PD-mul-ILP**  
  ILP formulation for multiple non-homologous proteins, with brute-force cross-hybridization detection and forbidden-pair constraints.  
  Documentation: [docs/PD-mul-ILP.md](docs/PD-mul-ILP.md)

## Setting Up the Tool

First, create a file named `gurobi.json` containing the details for the Gurobi license and put it in the main directory.

```json
{
  "WLSACCESSID": "XXXXX",
  "WLSSECRET": "XXXXX",
  "LICENSEID": 12345
}
```

Create a text file containing the protein names and their DNA coding sequences.  
Input formats that are specific to each PrimerDesigner variant are described in their corresponding documentation files.

Example:
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













