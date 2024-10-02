
# PrimerDesigner

## Introduction

This repository contains the source code for PrimerDesigner, a tool designed to find the most efficient primer set with complete coverage and no cross hybridization risk for protein synthesis using assembly PCR. It accompanies the paper titled "PrimerDesigner: Designing efficient primers for protein synthesis without cross-hybridization."

## Requirements

The tool has been tested with the following configuration on a Linux machine:
- Python 3.9.18
- GurobiPy 11.0.2
- NetworkX 3.2.1
- NumPy 1.24.3
- Pandas 2.1.3
- Biopython 1.75
- Primer3 2.0.1

## The Different Versions

**PrimerDesigner** has 4 different version options:

- **PD-single**
  - Finds the optimal primer set for a single protein by computing the shortest path in the primer graph.
  - Runs only on the first protein-coding sequence in the sequence file.

- **PD-mul-var**
  - Finds the optimal primer set for a number of variants of the same protein-coding sequence.
  - Runs only on the first protein-coding sequence in the sequence file.

- **PD-mul-nh**
  - Finds the optimal primer set for multiple non-homologous proteins by computing the shortest path in each protein primer graph.
  - If there are primers in the shortest path that cross-hybridize with previously selected primers, it removes these primers from the graph and computes the shortest path again.

- **Non relaxed**
  - Finds the optimal primer set for multiple non-homologous proteins.
  - Identifies all primers pairs that may cross-hybridize and adds them as forbidden pair ILP constraints.


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
python ./tool.py --file_path <file-path> --version <version> --output <output-file>
```
- **file_path**: The file path of the protein coding-sequences
- **version**: Specifies which version of the algorithm to run. The options are: PD-single, PD-mul-var, PD-mul-nh and non-ralxed (default:  PD-mul-var)
- **output**: The file path that the program output will be saved to.
  
The other arguments are optional and include the algorithm parameters:

- **primer_lmin, primer_lmax**: Minimum and maximum primer lengths (default: 18,30).
- **oligo_lmin, oligo_lmax**: Minimum and maximum oligonucleotides lengths (default: 195,205)
- **overlap_lmin, overlap_lmax**: Minimum and maximum overlap length between oligonucleotides  (default: 45,50)
- **allowed_overlap**: Allowed overlap between primer_pairs (default: 6)
- **num_proteins**: Number of variants of the same sequence - used for PD-mul-var version only (default: 3)
- **apply_threshold**: Boolean flag for applying primer quality threshold (default: False)
- **min_gc, max_gc**: Minimum and maximum threshold on gc content (default: 40,60)
- **min_tm, max_tm**: Minimum and maximum threshold on melting temperature tm (default: 58,65)
- **max_difference**: Maximum difference threshold between the forward and reverse primer in each pair (default: 3)
- **merge_bins**: Boolean flag for merging bins corresponding to identical non-overlapping sequences - used for the PD-mul-var version (default: False)


Example command:
```bash
python ./tool.py --file_path example_proteins.txt --version Non_relaxed --output run_output --primer_lmin 20 --primer_lmax 26 --oligo_lmin 180 --oligo_lmax 200
```


