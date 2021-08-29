from math import *
from isa import determine_atmosphere
from scipy.optimize import minimize
from numpy import zeros, mod, array, inf, arange
import random
from Settings_import import *
import sys

atmospheric_data = {'base_altitude': [0, 11000, 20000, 32000, 47000, 51000, 71000, 100000],
                    'base_temperature': [288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65],
                    'base_density': [1.225, 0.36391, 0.08803, 0.01322, 0.00143, 0.00086, 0.000064],
                    'base_pressure': [101325, 22632.1, 5474, 868.02, 110.91, 66.94, 3.96],
                    'lapse': [-0.0065, 0, 0.001, 0.0028, 0, -0.0028, -0.002]}


# ========================================= OPERATORS ============================================================

# Function to build geometry dictionary from gene values
def build_geometry(solution):
    solution_geometry = {'a': solution[0],
                         'b': solution[1],
                         'c': solution[2],
                         'd': solution[3],
                         'dz': solution[4],
                         'D_t': solution[5]}

    return solution_geometry

# Nozzle radius from geometric parameters + desired z_position
def radius_function(a, b, c, d, z):
    return (a + ((b + 1000 * c * z)**0.5 / d)) * 0.001


# Randomization function, used for generating initial population
def randomize_solution(input_solution):

    new_solution = []
    for i in range(GA_settings[6]):
        new_entry = input_solution[geometry_keys[i]] + (- 0.5 + 1*random.random())*input_solution[geometry_keys[i]]  # up to 50% variation
        new_solution.append(new_entry)

    if len(new_solution) < len(geometry):
        missing_genes = len(geometry) - len(new_solution)

        while missing_genes > 0:
            new_solution.append(geometry[geometry_keys[len(new_solution)]])
            missing_genes = len(geometry) - len(new_solution)

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


# The single point cross-over takes two chromosome instances and returns their two offspring
# NOTE: DOES NOT SUPPORT VARIABLE NUMBER OF UNSTABLE GENES
def single_point_crossover(chromosome_1, chromosome_2):

    # pull parent genes from input chromosomes
    parent_1 = chromosome_1.genes
    parent_2 = chromosome_2.genes

    # determine cross over point
    #cross_over_point = random.randint(0, len(parent_1))
    cross_over_point = 3

    # perform cross over
    for gene in range(cross_over_point, len(parent_1)):
        parent_1[gene], parent_2[gene] = parent_2[gene], parent_1[gene]

    # generate offspring chromosomes
    offspring = [Chromosome(parent_1, engine_properties), Chromosome(parent_2, engine_properties)]

    return offspring


# Uniform cross-over
def uniform_crossover(chromosome_1, chromosome_2):

    # pull parent genes from input chromosomes
    parent_1 = chromosome_1.genes
    parent_2 = chromosome_2.genes

    gene = 0

    while gene < GA_settings[6]:
        if mod(gene, 2) == 0:
            parent_1[gene], parent_2[gene] = parent_2[gene], parent_1[gene]

        gene += 1

    offspring = [Chromosome(build_geometry(parent_1), engine_properties), Chromosome(build_geometry(parent_2), engine_properties)]

    return offspring


# Function to evaluate fitness of population
def evaluate_fitness(population):

    # Generate list of apogees
    apogee_list = []
    for solution in population:
        # print('Calculating fitness of Chromosome ', solution)
        apogee = solution.objective_function_2(simulation_settings)['apogee']
        apogee_list.append(apogee)

    return apogee_list


# Roulette selection
def roulette_selection(population, evaluation, target_count):

    selection = []

    # Remove best performer and add it to the selection
    elite_chromosome = population[0]
    for chromosome in population:
        if chromosome.apogee > elite_chromosome.apogee:
            elite_chromosome = chromosome

    # Reset all elite statuses
    for chromosome in population:
        chromosome.is_elite = False

    elite_chromosome.is_elite = True

    selection.append(elite_chromosome)

    #evaluation.pop(elite_chromosome_index)
    #population.pop(elite_chromosome_index)

    # Substract lowest performer to increase relative differences
    min_performance = min(evaluation)
    for performance in evaluation:
        performance -= min_performance

    # Calculate fitness probabilities
    fitness_probabilities = []
    for apogee in evaluation:
        fitness_probabilities.append(apogee / sum(evaluation))

    # Calculate selection probabilities
    selection_probabilities = []
    iterator = 0

    # Make list of selection probabilities (selecting one less, since best performer was passed through earlier)
    while iterator < target_count - 1:
        selection_probabilities.append(random.randint(0, 100) / 100)
        iterator += 1

    # Make list of cumulative fitness probabilities
    cumulative_fitness = zeros(len(fitness_probabilities))

    for i in range(len(fitness_probabilities)):
        cumulative_fitness[i] = cumulative_fitness[i - 1] + fitness_probabilities[i]

    # Apply selection probabilities to cumulative fitness probabilities

    for selector in range(len(selection_probabilities)):  # selector running through each selection probability

        selected = False  # Boolean to indicate whether a chromosome was selected in current loop
        iterator = 0  # iterator running through each chromosome

        while not selected:

            if selection_probabilities[selector] <= cumulative_fitness[iterator]:

                # selection probability falls in the selection range of the current chromosome
                selection.append(population[iterator])
                selected = True

            elif selection_probabilities[selector] == 1.0:
                selection.append(population[-1])
                selected = True

            else:
                # selection probability does not fall in selection range of current chromosome, move on to next
                iterator += 1

    return selection


# Ranked selection
# TODO: unit list ranked selection
def ranked_selection(population, evaluation, target_count):

    selection = []
    iterator = 0
    temp_evaluation = []
    temp_evaluation += evaluation

    while iterator < target_count:
        index = temp_evaluation.index(max(temp_evaluation))         # find index of best performer
        selection.append(population[index])                         # add best performer to selection
        temp_evaluation.pop(index)                                  # remove previous max from evaluation list

        iterator += 1

    return selection


# Boundary condition checker
def boundary_conditions_ok(chromosome):

    # Maximum exit radius
    if chromosome.R_e >= design_constraints['Re_max']:
        return False

    # Maximum throat diameter
    elif chromosome.genes[-1] >= design_constraints['Rt_max']:
        return False

    # positive throat
    elif chromosome.genes[-1] <= 0.06:
        return False

    # Maximum length
    elif chromosome.genes[4] >= design_constraints['dz_max']:
        return False

    # Minimum expansion ratio
    elif chromosome.epsilon < design_constraints['epsilon_min']:
        return False

    # All conditions met
    else:
        return True


# ================================================= CLASSES ======================================================

# TODO: constrain nozzle divergent length dz
class Chromosome:

    def __init__(self, geometry, engine_properties):

        # initial conditions: -1.3360, 2.3935, 0.0888, 0.05, 0.1015
        self.genes = [geometry['a'], geometry['b'], geometry['c'], geometry['d'], geometry['dz'], geometry['D_t']]

        self.gamma = engine_properties[0]
        self.t_thrust = engine_properties[2]
        self.A_t = pi * self.genes[5] ** 2 / 4
        self.p_c = engine_properties[4]
        self.performance_loss_factor = 1
        self.pressure_ratio_value = 0
        self.apogee = 0
        self.z_max = 0
        self.z_min = 0
        self.epsilon = 8
        self.R_e = 0
        self.evaluated = False
        self.time = []
        self.data = []
        self.is_elite = False

    # NOZZLE CHARACTERISTICS METHODS

    # Determine graphite - zirconium transition point
    def transition(self):
        r = 0
        R_transition = 72.15 / 1000.
        #z_min = 0
        self.z_max = 0

        # Determine z_min(z_value for which the graphite - zirconium transition radius is reached)
        # while abs(r - R_transition) > 0.5 / 1000.:
        #     r = radius_function(self.genes[0], self.genes[1], self.genes[2], self.genes[3], z_min)
        #     z_min += 0.005 / 1000.

        z_min = ((R_transition - self.genes[0])**2 * self.genes[3]**2 - self.genes[1])/(1000*self.genes[2])


        self.z_min = z_min
        self.z_max = self.z_min + self.genes[4]

    # Vanderkerckhove function
    def vanderkerckhove(self):
        return sqrt(self.gamma)*(2/(self.gamma + 1))**((self.gamma + 1)/(2*(self.gamma - 1)))

    # Determine mass flow rate
    def massflow(self):
        Gamma = self.vanderkerckhove()

        #return Gamma * self.p_c * self.A_t / sqrt(287 * 2170)
        return engine_properties[1]

    # Determine exit radius
    def exit_radius(self):
        # TODO: fix relation between throat diameter and added exit radius
        self.R_e = radius_function(self.genes[0], self.genes[1], self.genes[2], self.genes[3], self.z_max) + geometry['D_t']/2
        return self.R_e

    # Determine expansion ratio
    def expansion_ratio(self):
        return (pi * self.exit_radius() ** 2)/self.A_t

    # Determine the pressure ratio based on geometric and flow properties
    def pressure_ratio(self):

        Gamma = sqrt(self.gamma) * (2 / (self.gamma + 1)) ** ((self.gamma + 1) / (2 * (self.gamma - 1)))
        pressure_ratio_calculated = 0.001
        pressure_stepsize = 0.00001

        epsilon_calculated = 0
        self.epsilon = self.expansion_ratio()

        error = abs(self.epsilon - epsilon_calculated)
        threshold = 1

        while error > threshold:

            epsilon_calculated = Gamma / sqrt((2 * self.gamma / (self.gamma - 1)) * pressure_ratio_calculated ** (2 / self.gamma) * (
                    1 - pressure_ratio_calculated ** ((self.gamma - 1) / self.gamma)))
            error = abs(self.epsilon - epsilon_calculated)

            pressure_ratio_calculated = pressure_ratio_calculated + pressure_stepsize

        return pressure_ratio_calculated

    # Determine the performance lost due to radial thrust force component
    def performance_loss(self):

        # self.transition()

        # Determine the nozzle exit angle
        delta_z = 0.005  # nozzle length over which straight wall is assumed to calculate exit angle (i.e. equivalent tangent)
        slope = (radius_function(self.genes[0], self.genes[1], self.genes[2], self.genes[3], self.z_max) - radius_function(self.genes[0], self.genes[1], self.genes[2], self.genes[3], self.z_max - delta_z)) / delta_z
        exit_angle = atan(slope)

        loss_factor = 1 - (1 - cos(exit_angle)) / 2

        return loss_factor

    # Determine nozzle exhaust velocity
    def exhaust_velocity(self):
        # Determine flow exit velocity based on geometric and flow conditions

        R = 640.91
        T_c = 2170

        U_e = sqrt(2 * (self.gamma * R * T_c) / (self.gamma - 1) * (1 - self.pressure_ratio_value ** ((self.gamma - 1) / self.gamma))) * self.performance_loss_factor

        return U_e

    # Objective function for optimization. Runs simulation based on input vector x to determine performance.
    def objective_function_2(self, simulation_settings):

        # check boundary conditions
        if not boundary_conditions_ok(self):
            self.apogee = 1
            self.evaluated = True


        # if boundary conditions ok, and chromosome not yet evaluated, run simulation
        if not self.evaluated:

            prev_altitude = 0  # m (used to determine if rocket is still ascending or not)
            altitude = 1  # m
            t = 0  # s
            dt = simulation_settings[0]
            c_d = simulation_settings[1]
            v = 0  # m / s

            # Determine nozzle parameters from input

            # TODO: note that nozzle mass changes are neglected to speed up the simulation speed
            self.transition()
            #nozzle_mass = self.objective_function_1()  # nozzle_properties returns [nozzle_mass, R_e].
            mass = 350 #+ nozzle_mass  # kg

            A_e = self.A_t * self.expansion_ratio()
            self.pressure_ratio_value = self.pressure_ratio()
            p_e = self.p_c * self.pressure_ratio_value
            m = self.massflow()


            engine_on = True
            ascending = True

            # Start simulation loop

            while ascending:  # or while t_simulation

                # Determine current atmospheric conditions
                atm = determine_atmosphere(altitude, atmospheric_data)
                p_a = atm['pressure']
                rho = atm['density']

                # Calculate net force and acceleration
                if t < self.t_thrust:
                    self.performance_loss_factor = self.performance_loss()
                    F_t = self.performance_loss_factor * m * self.exhaust_velocity() + (p_e - p_a) * A_e

                else:
                    engine_on = False
                    F_t = 0

                F_d = c_d * 0.5 * rho * v ** 2 * ((100. / 1000.) ** 2 * pi)
                F_w = mass * 9.80665

                F = F_t - F_d - F_w
                a = F / mass

                # Calculate distance travelled
                d = v * dt + 0.5 * a * dt ** 2
                altitude = altitude + d

                if altitude < prev_altitude:
                    ascending = False

                else:
                    prev_altitude = altitude

                # Update variables
                if t < self.t_thrust:
                    mass = mass - m * dt


                t = t + dt
                v = v + a * dt
                # time.append(t)
                self.data.append(altitude)

            self.apogee = max(self.data)

            self.evaluated = True       # Record that this chromosome has been evaluated. Whenever its genes are modified
                                        # this is set to false again. And whenever a population is being evaluated,
                                        # chromosomes that have already been evaluated can be skipped.

        return {'time': self.time,
                'data': self.data,
                'apogee': self.apogee}

    # Mutation operator
    def mutate(self):

        mutation_percentage = GA_settings[5]  # percentage by which a gene may be mutated
        random_number = random.random()

        # perform mutation
        if not self.is_elite:
            if random_number < GA_settings[2]:
                modified_gene = random.randint(0, len(self.genes) - 1)  # Determine what gene to modify
                modification = - mutation_percentage + 2 * mutation_percentage * random.random()  # Determine how much to modify the gene


                if modified_gene < GA_settings[6]:
                    self.genes[modified_gene] += modification * self.genes[modified_gene]
                    self.evaluated = False  # set the evaluation status to false

                else:
                    self.mutate()


class Population:

    def __init__(self, population):

        self.population = population
        self.population_steadiness = GA_settings[0]
        self.crossover_prob = GA_settings[1]
        self.mutation_prob = GA_settings[2]
        self.population_size = GA_settings[3]
        self.genes = []
        self.offspring_count = 0
        self.offspring = []
        self.evaluation = []


    # Derive how many offspring will be replacing population members (and so how many parents are needed)
    # The higher population steadiness is, the fewer members will be replaced by offspring
    def population_replacement(self):

        # Determine size of population
        population_size = GA_settings[3]

        # Calculate exact amount of population members that will stay
        steady_chromosome_count = round(population_size*self.population_steadiness, 0)

        # Derive number of offspring required
        offspring_count = population_size - steady_chromosome_count

        # ensure even amount of offspring is calculated
        if mod(offspring_count, 2) != 0:
            offspring_count += 1

        self.offspring_count = offspring_count

        return offspring_count

    # Selection: Fitness proportional parent selection.
    def parent_selection(self):

        # Evaluate population
        #print('Evaluating fitness...')
        self.evaluation = evaluate_fitness(self.population)

        # Perform roulette selection to select parents
        parent_selection = roulette_selection(self.population, self.evaluation, self.offspring_count)

        return parent_selection

    # Reproduce and add offspring to population
    def reproduce(self):

        offspring_target = self.population_replacement()    #   calling population replacement will also update self.offspring_count
                                                            #   so it doesn't have to be recalled in the future
        offspring_counter = 0
        offspring = []
        parent_list = self.parent_selection()               # Also runs evaluation of performance

        while offspring_counter < offspring_target:

            index_1 = offspring_counter
            index_2 = index_1 + 1

            offspring = offspring + uniform_crossover(parent_list[index_1], parent_list[index_2])

            offspring_counter += 2

        # Add random solutions
        base_chromosome_index = random.randint(0, len(self.population) - 1)
        base_geometry = {'a': self.population[base_chromosome_index].genes[0],
                           'b': self.population[base_chromosome_index].genes[1],
                           'c': self.population[base_chromosome_index].genes[2],
                           'd': self.population[base_chromosome_index].genes[3],
                           'dz': self.population[base_chromosome_index].genes[4],
                           'D_t': self.population[base_chromosome_index].genes[5]}        # base geometry for initialization

        new_solution_space = initialize(base_geometry, GA_settings[7])

        for solution in new_solution_space:
            solution_geometry = build_geometry(solution)
            offspring.append(Chromosome(solution_geometry, engine_properties))  # add chromosome with new geometry to offspring

        # evaluate offspring
        self.evaluation += evaluate_fitness(offspring)

        # add offspring to current population
        self.population += offspring


    # Select survivors of this generation
    def select_survivors(self):

        # Generate population offspring
        self.reproduce()

        for chromosome in self.population:
            chromosome.mutate()

        # Select members of the next generation
        self.population = roulette_selection(self.population, self.evaluation, self.population_size)

        # Update evaluation (they are now stored in individual chromosomes, so no need to re-evaluate everything)
        iterator = 0
        while iterator < self.population_size:
            self.evaluation[iterator] = self.population[iterator].apogee

            iterator += 1




