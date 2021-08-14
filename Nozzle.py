from math import *
from isa import determine_atmosphere
import sympy as sy
from numpy import zeros
from Operators import *



def radius_function(a, b, c, d, z):
    return (a + ((b + 1000 * c * z)**0.5 / d)) * 0.001


atmospheric_data = {'base_altitude': [0, 11000, 20000, 32000, 47000, 51000, 71000, 100000],
                    'base_temperature': [288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65],
                    'base_density': [1.225, 0.36391, 0.08803, 0.01322, 0.00143, 0.00086, 0.000064],
                    'base_pressure': [101325, 22632.1, 5474, 868.02, 110.91, 66.94, 3.96],
                    'lapse': [-0.0065, 0, 0.001, 0.0028, 0, -0.0028, -0.002]}


class Chromosome:

    def __init__(self, geometry, engine_properties):
        self.a = geometry[0]
        self.b = geometry[1]
        self.c = geometry[2]
        self.d = geometry[3]
        self.dz = geometry[4]

        self.gamma = engine_properties[0]
        self.m = engine_properties[1]
        self.t_thrust = engine_properties[2]
        self.A_t = engine_properties[3]
        self.p_c = engine_properties[4]

    # NOZZLE CHARACTERISTICS METHODS

    def transition(self):
        z_min = 0
        r = 0
        R_transition = 72.15 / 1000.

        # Determine z_min(z_value for which the graphite - zirconium transition radius is reached)
        while abs(r - R_transition) > 0.1 / 1000.:
            r = radius_function(self.a, self.b, self.c, self.d, z_min)
            z_min = z_min + 0.001 / 1000.

        z_max = z_min + self.dz

        return [z_min, z_max]

    def exit_radius(self):
        return radius_function(self.a, self.b, self.c, self.d, self.transition()[1])

    def expansion_ratio(self):
        return (pi * self.exit_radius() ** 2)/self.A_t

    def pressure_ratio(self):
        # Determine the pressure ratio based on geometric and flow properties

        Gamma = sqrt(self.gamma) * (2 / (self.gamma + 1)) ** ((self.gamma + 1) / (2 * (self.gamma - 1)))
        pressure_ratio_calculated = 0.001
        pressure_stepsize = 0.00005
        epsilon_calculated = 0

        epsilon_true = self.expansion_ratio()
        error = abs(epsilon_true - epsilon_calculated)
        threshold = 1

        while error > threshold:
            epsilon_calculated = Gamma / sqrt((2 * self.gamma / (self.gamma - 1)) * pressure_ratio_calculated ** (2 / self.gamma) * (
                    1 - pressure_ratio_calculated ** ((self.gamma - 1) / self.gamma)))
            error = abs(epsilon_true - epsilon_calculated)

            pressure_ratio_calculated = pressure_ratio_calculated + pressure_stepsize

        return pressure_ratio_calculated

    def performance_loss(self):
        # Determine the performance lost due to radial thrust force component

        # Determine z_min (z_value for which the graphite - zirconium transition radius is reached). Note that the transition
        # point is currently fixed at 72.15/1000 by design. This is considered a fixed constraint, so only the curve past
        # transition can be changed.

        z_min = 0
        r = 0
        r_transition = 72.15 / 1000.

        # gradually increase z_min until the resulting radius matches the transition radius
        while abs(r - r_transition) > 0.01:
            r = radius_function(self.a, self.b, self.c, self.d, z_min)
            z_min = z_min + (0.01 / 1000.)

        z_max = z_min + self.dz

        # Determine the nozzle exit angle
        delta_z = 0.005
        slope = (radius_function(self.a, self.b, self.c, self.d, z_max) - radius_function(self.a, self.b, self.c, self.d, z_max - delta_z)) / delta_z
        exit_angle = atan(slope)

        loss_factor = 1 - (1 - cos(exit_angle)) / 2

        return loss_factor

    def exhaust_velocity(self):
        # Determine flow exit velocity based on geometric and flow conditions

        R = 640.91
        T_c = 2170

        # TODO: Figure out what Gamma was needed for. It was also not used in the function-based version of the code
        Gamma = sqrt(self.gamma) * (2 / (self.gamma + 1)) ** ((self.gamma + 1) / (2 * (self.gamma - 1)))

        p_ratio = self.pressure_ratio()

        U_e = sqrt(2 * (self.gamma * R * T_c) / (self.gamma - 1) * (1 - p_ratio ** ((self.gamma - 1) / self.gamma))) * self.performance_loss()

        return U_e

    # OBJECTIVE FUNCTION METHODS

    def objective_function_1(self):
        # Determine nozzle volumetric parameters

        # Determine z_min(z_value for which the graphite - zirconium transition radius is reached)
        z_min = self.transition()[0]
        z_max = self.transition()[1]

        # Determine results through integration
        z = sy.Symbol("z")

        v_zirconium = (pi * 1. / 1000. * (1. / 1000. * self.dz + 2 * sy.integrate(radius_function(self.a,
                                                                                                        self.b,
                                                                                                        self.c,
                                                                                                        self.d,
                                                                                                        z),
                                                                                        (z, z_min, z_max))))  # cm3

        v_titanium = (pi * 4. / 1000. * (4. / 1000. * self.dz + 2 * sy.integrate(radius_function(self.a,
                                                                                                       self.b,
                                                                                                       self.c,
                                                                                                       self.d,
                                                                                                       z),
                                                                                       (z, z_min, z_max))))  # cm3

        # Return nozzle mass
        return v_zirconium * 6400 + v_titanium * 4420

    def objective_function_2(self, simulation_settings):
        # Objective function for optimization. Runs simulation based on input vector x to determine performance.
        # TODO: Assess input parameters for vector x (= Optimization parameters)

        # TODO: derive mass flow from throat area + pressure? adds a new geometric parameter to vary
        prev_altitude = 0  # m (used to determine if rocket is still ascending or not)
        altitude = 1  # m
        t = 0  # s
        dt = simulation_settings[0]
        p_c = self.p_c  # 15 * 10 ^ 5  # Pa
        c_d = simulation_settings[1]
        v = 0  # m / s

        # Determine nozzle parameters from input

        nozzle_mass = self.objective_function_1()  # nozzle_properties returns [nozzle_mass, R_e].
        mass = 350 + nozzle_mass  # kg
        # TODO: Change function outputs to dictionaries for clarity
        A_e = self.A_t * self.expansion_ratio()
        p_e = p_c * self.pressure_ratio()

        data = []
        time = []
        ascending = True
        # Start simulation loop

        while ascending:  # or while t_simulation

            # Determine current atmospheric conditions
            atm = determine_atmosphere(altitude, atmospheric_data)
            p_a = atm['pressure']
            rho = atm['density']

            # Calculate net force and acceleration
            if t < self.t_thrust:
                F_t = self.performance_loss() * self.m * self.exhaust_velocity() + (p_e - p_a) * A_e

            else:
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
                mass = mass - self.m * dt


            t = t + dt
            v = v + a * dt
            time.append(t)
            data.append(altitude)

        # plot(time, data)
        # fprintf('performance evaluated.\n');
        apogee = max(data)

        simulation = {'time': time,
                      'data': data,
                      'apogee': apogee,
                      'nozzle mass': nozzle_mass}

        return simulation

    # Chromosome information

    def genes(self):
        return [self.a, self.b, self.c, self.d, self.dz]





class Population:

    def __init__(self, population):

        self.population = population
        self.population_steadiness = GA_settings[0]
        self.crossover_prob = GA_settings[1]
        self.mutation_prob = GA_settings[2]
        self.genes = []
        self.offspring_count = 0


    # Function to derive how many offspring will be replacing population members (and so how many parents are needed)
    # The higher population steadiness is, the less members will be replaced by offspring
    def population_replacement(self):

        # Determine size of population
        population_size = GA_settings[3]

        # Calculate exact amount of population members that will stay
        steady_chromosome_count = round(population_size*self.population_steadiness, 0)

        # Derive number of offspring required
        offspring_count = population_size - steady_chromosome_count

        print('')
        print('population size: ', population_size)
        print('population steadiness: ', self.population_steadiness)
        print('number of offspring to be calculated: ', offspring_count)
        print('')

        self.offspring_count = offspring_count

        return offspring_count

    # Selection: Fitness proportional parent selection.
    def parent_selection_roulette(self):

        # Generate list of apogees
        apogee_list = []
        for solution in self.population:
            #print('Calculating fitness of Chromosome ', solution)
            apogee = solution.objective_function_2(simulation_settings)['apogee']
            apogee_list.append(apogee)

        # Calculate fitness probabilities
        total = sum(apogee_list)
        fitness_probabilities = []


        for apogee in apogee_list:

            fitness_probability = apogee / total
            fitness_probabilities.append(fitness_probability)


        # Calculate selection probabilities
        selection_probabilities = []
        iterator = 0


        # Select the amount of required offspring as parents (a couple makes two offspring)
        while iterator < self.offspring_count:

            selection_probability = random.randint(0, 100) / 100
            selection_probabilities.append(selection_probability)

            iterator += 1


        # Make list of cumulative fitness probabilities
        cumulative_fitness = zeros(len(fitness_probabilities))
        #print("selec. probs.: ", selection_probabilities)

        for i in range(len(fitness_probabilities)):

            cumulative_fitness[i] = cumulative_fitness[i - 1] + fitness_probabilities[i]

        #print("cumulative fitness: ", cumulative_fitness)
        # Apply selection probabilities to cumulative fitness probabilities
        parent_selection = []


        for selector in range(len(selection_probabilities)):  # selector running through each selection probability

            selected = False  # Boolean to indicated whether a chromosome was selected in current loop
            iterator = 0  # iterator running through each chromosome

            while not selected:

                if selection_probabilities[selector] < cumulative_fitness[iterator]:

                    # selection probability falls in the selection range of the current chromosome
                    parent_selection.append(self.population[iterator])
                    selected = True

                else:
                    # selection probability does not fall in selecion range of current chromosome, move on to next
                    iterator += 1




        print('Selected parents:', len(parent_selection))

        return parent_selection


    # Cross-over: Cutting genomes of different solutions at a random spot and combining them for a new solution.
    # TODO: check if method is still static

    # Reproduce and add offspring to population
    def reproduce(self):

        offspring_target = self.population_replacement()    #   calling population replacement will also update self.offspring_count
                                                            #   so it doesn't have to be recalled in the future
        offspring_counter = 0
        offspring = []
        parent_list = self.parent_selection_roulette()

        while offspring_counter < offspring_target:

            index_1 = offspring_counter
            index_2 = index_1 + 1

            offspring = offspring + single_point_crossover(parent_list[index_1], parent_list[index_2])

            offspring_counter += 2

        # add offspring to current population
        self.population = self.population + offspring


    # TODO: add mutation operator

    # TODO: add selection of survivors (# of selected survivors = population size)
#class Generational_Population:
#
#    def __init__(self):
