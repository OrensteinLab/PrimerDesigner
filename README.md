
# PrimerDesigner

## Introduction
 
This repository contains the source code for the PrimerDesigner project. 
The main goal of this project is to create an algorithm for primer design for protein synthesis using microarray oligonucleotide probes. 
Protein synthesis using microarray oligonucleotide probes is an efficient protein synthesis technique which enables high throughput protein mutagensis experiments.
The algorithm aims to find the most efficient primer set with complete coverage and no cross hybridization risk.
The program uses a primer graph to represent all valid forward and reverse primer combinations.
It also uses Integer Linear Programming (ILP) with specific  overlap and single path constraints to find the best path in the primer graph.


## Requirments

- PrimerDesinger requires Python 3.7 or higher.<br>
- To run PrimerDesinger, you will need a local installation of Gurobi with an appropriate license (academic licenses are provided for free direct from Gurobi).<br>
- You will also need to install gurobipy and make sure that you are able to import it from within your local environment.<br>
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
2. Open the PrimerDeisgner.ipynb notebook using Jupyter or any compatible environment.
<br>
<br>
3. Run the Notebook:
   Execute the cells in the notebook to track timings and analyze the primer creation process.


## Instructions

- Input your gurobi ID and key in the "Creating Gurobipy Solver" section. <br>
 <br>
<img width="341" alt="image" src="https://github.com/OrensteinLab/PrimerDesigner/assets/118123854/a266ad50-bb01-4f6e-a8b0-c2e236df913b"><br>
<br>
<br>

- Define upstream, mutreg and downstream regions based on the on your protein coding sequence <br>
 <br>
<img width="656" alt="image" src="https://github.com/OrensteinLab/PrimerDesigner/assets/118123854/9d418c7b-425d-48ca-a510-d527e344f04f"><br>
<br>
<br>

-You can adjust the algorithm's parameters. The parameters include the primer length range, overlap length range, oligonucleotide length range, number of proteins and maximum allowed overlap between primers <br> <br>
<img width="475" alt="image" src="https://github.com/OrensteinLab/PrimerDesigner/assets/118123854/a645bac2-9054-47c6-9327-824245db8048"><br>
<br>
<br>

- You can choose to apply threshold by setting the apply_threshold flag to True. You can adjust the min and max threshold on melting temeprature (tm) and gc contecnt. In addition, you can set the maxmimum allowed difference between forward and reverse primer melting temperatures.<br> <br>
<img width="474" alt="image" src="https://github.com/OrensteinLab/PrimerDesigner/assets/118123854/85ab62cf-c243-4067-9e5d-370bb7dc5b0f">
<br>
<br>



