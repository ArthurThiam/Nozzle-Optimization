from math import *


def radius_function(a, b, c, d, z):
    return (a + ((b + 1000 * c * z)**0.5 / d)) * 0.001


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


