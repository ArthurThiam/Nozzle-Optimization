from Nozzle import *
import numpy as np
import scipy.optimize as sp_optimize

# Define fixed parameters
fixed_params = {'z': 0.3,
                'b': 2.3935}



# Define constraints
ineq_constraints = {
    'type': 'ineq',
    'fun': lambda x: np.array([sqrt(2.3935 + 300*x[0]) - sqrt(2.3935 + 195.5*x[0]) - 102*x[1]],
                              [2.3935 + 300 * x[1]]),

    'jac': lambda x: np.array([[150 / sqrt(2.3935 + 300 * x[0]) - 99.25 / sqrt(2.3935 + 198.5 * x[0]), -102],
                               [0, 300]])
    }

eq_constraints = {
    'type': 'eq',
    'fun': lambda x: np.array([sqrt(2.3935 + 300*x[0]) - 103.336*x[1],
                               sqrt(2.3935 + 198.5*x[0]) - 73.486*x[1]]),

    'jac': lambda x: np.array([[150 / sqrt(300 * x[0] + 2.3935), -103.336],
                              [99.25 / sqrt(198.5 * x[0] + 2.3935), -73.486]])
}

bounds = ((0, 100), (0, 100))


# Define objective function
def objective_function(x):
    f = x[0] / (2 * x[1] * sqrt(2.3935 + 1000 * x[0] * 0.3))
    return float(f)


# Derivate of objective function
def objective_function_der(x):
    dx1 = (0.25*(x[0] + 0.0159567))/(x[1]*(x[0] + 0.00797833) * sqrt(300*x[0] + 2.3935))
    dx2 = - x[0] / (2 * sqrt(300*x[0] + 2.3935) * x[1] ** 2)

    func_derivative = [dx1, dx2]

    return func_derivative


# Optimize
max_count = 3
x0 = np.array([0.5, 0.02])
solutions = []
iterator = 0

while iterator != max_count:
    res = sp_optimize.minimize(objective_function, x0, args=(), jac=objective_function_der, bounds=bounds, method='SLSQP',
                           constraints=eq_constraints, options={'disp': False, 'maxiter': 10000}, )

    solutions.append([[x0[0], x0[1]], res.x])
    x0[0] *= 1.1
    x0[1] *= 1.1


    iterator += 1

print("Solutions:")
for solution in solutions:
    print(solution)

solution = [solutions[0][1][0], solutions[0][1][1]]

a = geometry['a']
b = geometry['b']
c = solution[0]
d = solution[1]

R_t = radius_function(geometry['a'], geometry['b'], solution[0], solution[1], (0.3 - geometry['dz']))
R_e = radius_function(geometry['a'], geometry['b'], solution[0], solution[1], 0.3)

print('\nTransition Radius: ', 1000*R_t, ' mm')
print('Exit Radius: ', 1000*R_e, ' mm')

ineq_cons_1 = False
ineq_cons_2 = False
eq_cons_1 = False
eq_cons_2 = False

# c = 0.5
# d = 0.02


if sqrt(2.3935 + 300*c) - sqrt(2.3935 + 195.5*c) - 102*d > 0:
    ineq_cons_1 = True

else:
    print(sqrt(2.3935 + 300*c) - sqrt(2.3935 + 195.5*c) - 102*d)

if 2.3935 + 300 * c > 0:
    ineq_cons_2 = True

print('\nInequality constraint 1 satisfied: ', ineq_cons_1)
print('Inequality constraint 2 satisfied: ', ineq_cons_2)
