from Nozzle import *
from scipy.optimize import minimize, Bounds, LinearConstraint
import plotly.graph_objects as go

# Define fixed parameters
fixed_params = {'z': 0.3,
                'b': 2.3935}


# Define initial solution
x0 = array([geometry['c'], geometry['d']])

# Define bounds
bounds = Bounds()


# Define objective function
def objective_function(solution):
    return solution[2] / (2 * solution[3] * sqrt(fixed_params['b'] + 1000 * solution[2] * fixed_params['z']))


# Perform minimization
res = minimize(objective_function, x0, method='nelder-mead', options={'xatol': 1e-5, 'disp': True})


# z = 0
# z_list = []
# data = []
#
# while z < fixed_params['z']:
#     derivative = objective_function(geometry, z)
#     data.append(derivative)
#     z_list.append(z)
#
#     z += 0.0001
#
# fig = go.Figure(data=go.Scatter(x=z_list, y=data))
# fig.show()