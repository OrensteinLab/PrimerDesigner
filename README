
# PrimerDesigner

## Overview
 
This repository contains the source code for the PrimerDesigner project. The main focus of this project is to create an algorithm that designs the most efficient primer set for pcr mutagensis experiments with complete coverage and no cross hybridization risk. The algorithm uses a primer graph and Integer Linear Programming (ILP)


# Requirments

PrimerDesinger requires Python 3.7 or higher. To run PrimerDesinger, you will need a local installation of Gurobi with an appropriate license (academic licenses are provided for free direct from Gurobi). You will also need to install gurobipy and make sure that you are able to import it from within your local environment. The installs are located at the top of the colab file. Note also that a large RAM may be required due to the space complexity of ILP solvers. For the full sequence provided, a machine with 70 GB of RAM is sufficient. 

## Mutagensis.ipynb

The import section installs and imports the necessary libraries. The setup subsection establishes the data structures required to carry out primer design (data_frame for primer costs, the actual DNA sequence). Under Algorithm: first a primer graph is created, then the user must input their gurobipy ID and key in the "Creating Gurobipy Solver" section. Next, a few hyperparameters are initalized. Finally, we construct a few data structures necessary for the ILP, and then the ILP is executed. 

## Tracking_Timings_Creation.ipynb

This program has similar inital setup to Mutagenesis, but theit also iterates through many different coding sequence sequence lengths & provides primer design outputs for each one. After it finishes iterating, it will save all of the relevant output data to a list and save it to a csv file.  In addition, the program keeps track of time and memory usage for every part and records it for the user

## Getting Started

1. **Clone the Repository:**
   ```bash
   git clone https://github.com/OrensteinLab/PrimerDesigner.git
   cd PrimerDesigner

2. Open the Tracking_Timings_Creation.ipynb notebook using Jupyter or any compatible environment.

3. Run the Notebook:
Execute the cells in the notebook to track timings and analyze the primer creation process.


