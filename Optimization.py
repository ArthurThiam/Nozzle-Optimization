#from ObjectiveFunctions import *
from Nozzle import *
import plotly.express as px
import random
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


# Selection: Randomly selecting the solutions to carry over from the solution space to the next generation.
def roulette_wheel_selection(population):

    # generate list of apogees
    apogee_list = []
    for solution in population:

        apogee = solution.objective_function_2(simulation_settings)['apogee']
        apogee_list.append(apogee)

    # calculate fitness and selection probabilities
    total = sum(apogee_list)
    fitness_probabilities = []
    selection_probabilities = []

    for apogee in apogee_list:
        fitness_probability = apogee/total
        selection_probability = random.randint(0,100)/100

        fitness_probabilities.append(fitness_probability)
        selection_probabilities.append(selection_probability)

    # make list of cumulative fitness probabilities
    cumulative_fitness = zeros(len(fitness_probabilities))
    for i in range(len(fitness_probabilities)):
        cumulative_fitness[i] = cumulative_fitness[i-1] + fitness_probabilities[i]

    # apply selection probabilities to cumulative fitness probabilities
    found = False
    selection = []

    for iterator in range(len(selection_probabilities)):
        if selection_probabilities[iterator] < fitness_probabilities[iterator]:
            iterator += 1

        else:
            found = True
            selection.append(population[iterator])



    #print(selection)

    return selection


# Cross-over: Cutting genomes of different solutions at a random spot and combining them for a new solution.
def single_point_crossover(parent_1, parent_2):

    cross_over_point = random.randint(0, len(parent_1))

    for gene in range(cross_over_point, len(parent_1)):
        parent_1[gene], parent_2[gene] = parent_2[gene], parent_1[gene]

    offspring = [parent_1, parent_2]

    return offspring

# Elitism: Adding the X best solutions of this generation to the next generation.

# Mutation: Introducing a random variation in the genomes of the new generation.


# ============================================== MAIN CODE =============================================================

# Initialize population
input_solution = geometry
solution_space = initialize(geometry, 5)

print('Initial solution space: ', solution_space)


# Create population of chromosomes from initial solution space
population = []
for solution in solution_space:
    population.append(Chromosome(solution, engine_properties))

# Select fittest chromosomes
roulette_wheel_selection(population)

# Print apogee of every solution
#for solution in population:
#    print(solution.objective_function_2(simulation_settings)['apogee'])




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
