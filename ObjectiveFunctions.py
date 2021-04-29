from Nozzle import *
from isa import determine_atmosphere
import sympy as sy

atmospheric_data = {'base_altitude': [0, 11000, 20000, 32000, 47000, 51000, 71000, 100000],
                    'base_temperature': [288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65],
                    'base_density': [1.225, 0.36391, 0.08803, 0.01322, 0.00143, 0.00086, 0.000064],
                    'base_pressure': [101325, 22632.1, 5474, 868.02, 110.91, 66.94, 3.96],
                    'lapse': [-0.0065, 0, 0.001, 0.0028, 0, -0.0028, -0.002]}

# Objective function 1 calculates the mass of the nozzle.
# TODO: individual material volume can be used to derive a nozzle cost function
def objective_function_1(chromosome):
    # Determine nozzle volumetric parameters

    # Determine z_min(z_value for which the graphite - zirconium transition radius is reached)
    z_min = chromosome.transition()[0]
    z_max = chromosome.transition()[1]

    # Determine results through integration
    z = sy.Symbol("z")

    v_zirconium = (pi * 1. / 1000. * (1. / 1000. * chromosome.dz + 2 * sy.integrate(radius_function(chromosome.a,
                                                                                                    chromosome.b,
                                                                                                    chromosome.c,
                                                                                                    chromosome.d,
                                                                                                    z), (z, z_min, z_max))))  # cm3

    v_titanium = (pi * 4. / 1000. * (4. / 1000. * chromosome.dz + 2 * sy.integrate(radius_function(chromosome.a,
                                                                                                   chromosome.b,
                                                                                                   chromosome.c,
                                                                                                   chromosome.d,
                                                                                                   z), (z, z_min, z_max))))  # cm3

    # Return nozzle mass
    return v_zirconium * 6400 + v_titanium * 4420


def objective_function_2(chromosome):
    # Objective function for optimization. Runs simulation based on input vector x to determine performance.
    # TODO: Assess input parameters for vector x (= Optimization parameters)



    # TODO: derive mass flow from throat area + pressure? adds a new geometric parameter to vary
    altitude = 1  # m
    t = 0  # s
    dt = 0.1
    p_c = 15 * 10 ^ 5  # Pa
    c_d = 1.5
    v = 0  # m / s

    # Determine nozzle parameters from input

    nozzle_mass = objective_function_1(chromosome)  # nozzle_properties returns [nozzle_mass, R_e].
    mass = 350 + nozzle_mass  # kg
    # TODO: Change function outputs to dictionaries for clarity
    A_e = chromosome.A_t * chromosome.expansion_ratio()
    p_e = p_c * chromosome.pressure_ratio()

    data = []
    time = []

    # Start simulation loop

    while altitude > 0: # or while t_simulation

        # Determine current atmospheric conditions
        atm = determine_atmosphere(altitude, atmospheric_data)
        p_a = atm['pressure']
        rho = atm['density']

        # Calculate net force and acceleration
        if t < chromosome.t_thrust:
            F_t = chromosome.performance_loss() * chromosome.m * chromosome.exhaust_velocity() + (p_e - p_a) * A_e

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
        if t < chromosome.t_thrust:
            mass = mass - chromosome.m * dt

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
