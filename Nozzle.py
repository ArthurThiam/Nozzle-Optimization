from math import *
from isa import determine_atmosphere
import sympy as sy
from Settings_import import *


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

        # Start simulation loop

        while altitude > 0:  # or while t_simulation

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

            if altitude < 0:
                altitude = 0

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
