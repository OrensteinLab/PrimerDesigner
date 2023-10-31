Requirements: 

PrimerDesinger requires Python 3.7 or higher.

To run PrimerDesinger, you will need a local installation of Gurobi with an appropriate lisence (academic licenses are provided for free direct from Gurobi). You will also need to install gurobipy and make sure that you are able to import it from within your local environment. The installs are located at the top of the colab file. Note also that a large RAM may be required due to the space complexity of ILP solvers. For the full sequence provided, a machine with 70 GB of RAM is sufficient. 

Mutagenesis: The import section installs and imports the necessary libraries. The setup subsection establishes the data structures required to carry out primer design (data_frame for primer costs, the actual DNA sequence. Under Algorithm: first a primer graph is created, then the user must input their gurobipy ID and key in the "Creating Gurobipy Solver" section. Next, a few hyperparameters are initalized. Finally, we construct a few data structures necessary for the ILP, and then the ILP is executed. 

Tracking Timings: Similar inital setup to Mutagenesis, but the program will iterate through many different sequence lengths & provide outputs for each one. After it finishes iterating, it will save all of the relevant output data to a list and save it to a csv file. 

