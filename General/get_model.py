import json
import gurobipy as gp

def get_model():
    with open('gurobi.json', 'r') as json_file:
        params = json.load(json_file)

    env = gp.Env(params=params)

    # Create the model within the Gurobi environment
    model = gp.Model('min-sum', env=env)

    return model