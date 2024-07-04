
# PrimerDesigner

## Introduction
This repository contains the source code for the PrimerDesigner project. 
The main goal of this project is to create an algorithm for primer design for protein synthesis using microarray oligonucleotide probes. 
The algorithm aims to find the most efficient primer set with complete coverage and no cross hybridization risk.
The program uses a primer graph to represent all valid forward and reverse primer combinations.
It also uses Integer Linear Programming (ILP) with specific forbidden pair and single path constraints to find the best path in the primer graph.


## Requirments

- PrimerDesigner requires Python 3.7 or higher.<br>
- To run PrimerDesigner, you will need a local installation of Gurobi with an appropriate license (academic licenses are provided for free direct from Gurobi).<br>
- You will also need to make sure that you are able to import Gurobi from within your local environment.<br>
- The installs are located at the top of the colab file. Note also that a large RAM may be required due to the space complexity of ILP solvers. <br>
- For the full sequence provided, a machine with 70 GB of RAM is sufficient. <br>

## Getting Started
<br>
1. Clone the Repository:

   ```
   git clone https://github.com/OrensteinLab/PrimerDesigner.git
   cd PrimerDesigner
   ```
<br>
2. Open the PrimerDeisgner_relaxed.ipynb or PrimerDeisgner_non_relaxed.ipynb notebook using Jupyter or any compatible environment.
<br>
<br>
3. Run the Notebook:
   Execute the cells in the notebook to track timings and analyze the primer design process.


## Instructions

- Input your gurobi ID and key in the beginning of the Gurobi Setup section <br><br>
```
# input gurobi license id and key
  params = {
  "WLSACCESSID":"",
  "WLSSECRET":"",
  "LICENSEID":
  }
  env = gp.Env(params=params)

```
- Define upstream, mutreg and downstream regions based on your protein coding sequence in the "Full Sequence" section <br><br>

```
# define sequences 
upstream_nt = "ATTTGAATGTATTTAGAAAAATAAACAAATAGGGGTTCCGCGCACATTTCCCCGAAAAGTGCTAGTGGTGCTAGCCCCGCGAAATTAAT..."
mutreg_nt_full = "CAAAGCCCAGCACCTGCCGCAGCGCCTGCCCCTGCGGCACGTTCCATCGCAGCTACGCCTCCTAAACTGATCGTGGCAATTAGCGT..."
downstream_nt = "GGAGGAGGGTCTGGGGGAGGAGGCAGTGGCATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTG..."

```
- You can adjust the algorithm parameters in the "Parameters" section. The parameters include the primer length range, overlap length range, oligonucleotide length range, number of proteins and maximum allowed overlap between primers <br><br>

```
# algorithm parameters
# adjust algorithm parameters
primer_lmin, primer_lmax = 18, 30
overlap_lmin,overlap_lmax = 45,50
oligo_lmin,oligo_lmax = 195,205
allowed_overlap = 6
num_proteins=3

```
- You can choose to apply efficiency thresholds in the "Parameters" section by setting the apply_threshold flag to True. You can adjust the min and max thresholds on melting temeprature (tm) and gc content. In addition, you can set the maximum allowed tm difference between the forward and reverse primers in every pair .<br> 
<br>

```
# threshold params
apply_threshold= True # apply threshold flag.
min_gc=40
max_gc=60
min_tm=58
max_tm=65
max_difference=3 # max tm difference between forward & reverse primers

```
<br>


