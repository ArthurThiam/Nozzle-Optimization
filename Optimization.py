#from ObjectiveFunctions import *
from Nozzle import *
import time
import plotly.express as px
from numpy import zeros


# ============================== GENETIC ALGORITHM FUNCTIONS ==================================


# Randomization function, used for generating initial population
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


# ============================================== MAIN CODE =============================================================


# Initialize population
solution_space = initialize(geometry, GA_settings[3])

# Create population of chromosomes from initial solution space
gene_set = []
for solution in solution_space:
    gene_set.append(Chromosome(solution, engine_properties))

population = Population(gene_set)

# Perform reproduction (and time it)
start_time = time.time()

population.reproduce()

runtime = time.time() - start_time

print('New population: ', len(population.population))
print('')
print('total runtime [s]: ', runtime)
print('average runtime per chromosome [s/chromosome]:', runtime/GA_settings[3])


# ============================================== UNIT TESTING ==========================================================

# Create offspring
#solution_1 = population[0].genes()
#print(solution_1)

#solution_2 = population[1].genes()
#print(solution_2)

#print(single_point_crossover(solution_1, solution_2))



# Create genome instance
#chromosome_1 = Chromosome(geometry, engine_properties)

# Run simulation test
#simulation = objective_function_2(chromosome_1, simulation_settings)
#print(simulation['apogee'])
