
# PrimerDesigner

## Introduction
 
This repository contains the source code for the PrimerDesigner project. <br> <br>
The main goal of this project is to create an algorithm for primer design for protein synthesis using microarray oligonucleotide probes. <br><br>
Protein synthesis using microarray oligonucleotide probes is an efficient protein synthesis technique which enables high throughput protein mutagensis experiments.<br><br>
The algorithm aims to find the most efficient primer set with complete coverage and no cross hybridization risk.<br><br>
The program uses a primer graph to represent all valid forward and reverse primer combinations.<br><br>
It also uses Integer Linear Programming (ILP) with specific  overlap and single path constraints to find the best path in the primer graph.<br><br>


## Requirments

PrimerDesinger requires Python 3.7 or higher.<br><br>
To run PrimerDesinger, you will need a local installation of Gurobi with an appropriate license (academic licenses are provided for free direct from Gurobi).<br><br>
You will also need to install gurobipy and make sure that you are able to import it from within your local environment.<br><br>
The installs are located at the top of the colab file. Note also that a large RAM may be required due to the space complexity of ILP solvers. <br><br>
For the full sequence provided, a machine with 70 GB of RAM is sufficient. <br><br>

## Getting Started
<br>
1. Clone the Repository:

   ```
   git clone https://github.com/OrensteinLab/PrimerDesigner.git
   cd PrimerDesigner
   ```
<br>
<br>
2. Open the Tracking_Timings_Creation.ipynb notebook using Jupyter or any compatible environment.
<br>
<br>
3. Run the Notebook:
   Execute the cells in the notebook to track timings and analyze the primer creation process.



# Programs

**Mutagensis.ipynb**
<br>
* The import section installs and imports the necessary libraries. <br>
* The setup subsection establishes the data structures required to carry out primer design (data_frame for primer costs, the actual DNA sequence).  <br>
* Under Algorithm: first a primer graph is created, then the user must input their gurobi ID and key in the "Creating Gurobipy Solver" section. <br>
* Next, the model hyperparameters( including primer length, maximum overlap and oligonucleotide length) are initalized.  <br>
* The program construct a few data structures necessary for the ILP, such as a dictionary of bins which store primers with overlapping indexes.<br>
* Finally, ILP is executed using gurobi . 
<br>

**Tracking_Timings_Creation.ipynb**
<br>
* This program has similar inital setup to Mutagenesis, but it also iterates through many different coding sequence lengths & provides primer design outputs for each one.
* After it finishes iterating, it will save all of the relevant output data to a list and save it to a csv file.
* In addition, the program keeps track of time and memory usage for every part and records it for the user


