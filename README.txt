Requirements: 

PrimerDesinger requires Python 3.7 or higher.

To run PrimerDesinger, you will need a local installation of Gurobi with an appropriate lisence (academic licenses are provided for free direct from Gurobi). You will also need to install gurobipy and make sure that you are able to import it from within your local environment. The installs are located at the top of the colab file. 

Mutagenesis: The import section installs and imports the necessary libraries. The setup subsection establishes the data structures required to carry out primer design (data_frame for primer costs, the actual DNA sequence. Under Algorithm: first a primer graph is created, then the user must input their gurobipy ID and key in the "Creating Gurobipy Solver" section. Next, a few hyperparameters are initalized. Finally, we construct a few data structures necessary for the ILP, and then the ILP is executed. 

