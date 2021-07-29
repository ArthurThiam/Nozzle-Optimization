from ObjectiveFunctions import *
import configparser
import plotly.express as px
import random


# ============================= IMPORT SETTINGS ===============================================


config = configparser.ConfigParser()
config.read('settings.ini')

geometry = [config.getfloat('Geometry', 'a'),
            config.getfloat('Geometry', 'b'),
            config.getfloat('Geometry', 'c'),
            config.getfloat('Geometry', 'd'),
            config.getfloat('Geometry', 'dz')]

engine_properties = [config.getfloat('Engine Properties', 'gamma'),
                     config.getfloat('Engine Properties', 'm'),
                     config.getfloat('Engine Properties', 't_thrust'),
                     config.getfloat('Engine Properties', 'A_t'),
                     config.getfloat('Engine Properties', 'p_c')]

simulation_settings = [config.getfloat('Simulation Settings', 'dt'),
                       config.getfloat('Simulation Settings', 'c_d')]


# ============================== GENETIC ALGORITHM FUNCTIONS ==================================


# Randomization function for generating initial population
def randomize_solution(input_solution):

    new_solution = []
    for i in range(len(input_solution)):
        new_entry = input_solution[i] + (- 0.5 + 1*random.random())*input_solution[i]  # up to 50% variation
        new_solution.append(new_entry)

    return new_solution


# Initialization: generate an initial population based on the initial solution.
def initialize(input_solution, population_size):
    initial_population = []
    iterator = 0

    # vary each parameter within the initial solution by a random percentage to form initial population
    while iterator < population_size:
        initial_population.append(randomize_solution(input_solution))
        iterator += 1

    return initial_population


# Selection: Randomly selecting the solutions to carry over from the solution space to the next generation.
def select_pair(evaluation):

    selected = []

    return selected

# Cross-over: Cutting genomes of different solutions at a random spot and combining them to for a new solution.

# Elitism: Adding the X best solutions of this generation to the next generation.

# Mutation: Introducing a random variation in the genomes of the new generation.






# ============================================== UNIT TESTING ==========================================================


# Create genome instance
chromosome_1 = Chromosome(geometry, engine_properties)

# Run simulation test
simulation = objective_function_2(chromosome_1, simulation_settings)
print(simulation['apogee'])
