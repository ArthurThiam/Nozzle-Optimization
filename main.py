from math import *
import sympy as sy

# Define atmospheric data
atmospheric_data = {'base_altitude': [0, 11000, 20000, 32000, 47000, 51000, 71000, 100000],
                    'base_temperature': [288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65],
                    'base_density': [1.225, 0.36391, 0.08803, 0.01322, 0.00143, 0.00086, 0.000064],
                    'base_pressure': [101325, 22632.1, 5474, 868.02, 110.91, 66.94, 3.96],
                    'lapse': [-0.0065, 0, 0.001, 0.0028, 0, -0.0028, -0.002]}


def determine_atmosphere(altitude, atmosphere):
    # Determine pressure and density at the current altitude. Works up to and excluding 100,000 m

    data_found = False
    i = 0

    while not data_found:  # for i in range(len(atmosphere['base_altitude'])):

        if atmosphere['base_altitude'][i] <= altitude < atmosphere['base_altitude'][i+1]:

            base_altitude = atmosphere['base_altitude'][i]
            base_temperature = atmosphere['base_temperature'][i]
            base_density = atmosphere['base_density'][i]
            base_pressure = atmosphere['base_pressure'][i]
            lapse = atmosphere['lapse'][i]

            data_found = True

        i += 1


    if lapse == 0:
        density = base_density*exp((-9.80665*0.0289644*(altitude-base_altitude))/(8.3144598*base_temperature))
        pressure = base_pressure*exp((-9.80665*0.0289644*(altitude-base_altitude))/(8.3144598*base_temperature))

    elif lapse != 0:
        density = base_density * (base_temperature / (base_temperature + lapse * (altitude - base_altitude))) ** (1 + (9.80665 * 0.0289644) / (8.3144598 * lapse))
        pressure = base_pressure * (base_temperature / (base_temperature + lapse * (altitude - base_altitude))) ** ((9.80665 * 0.0289644) / (8.3144598 * lapse))

    atmospheric_properties = {'base_altitude': base_altitude,
                              'altitude': altitude,
                              'pressure': pressure,
                              'density': density}

    return atmospheric_properties


def pressure_ratio(epsilon_true, gamma):
    # Determine the pressure ratio based on geometric and flow properties


    Gamma = sqrt(gamma) * (2 / (gamma + 1)) ** ((gamma + 1) / (2 * (gamma - 1)))
    pressure_ratio_calculated = 0.001
    pressure_stepsize = 0.00005
    epsilon_calculated = 0

    error = abs(epsilon_true - epsilon_calculated)
    threshold = 1

    while error > threshold:
        epsilon_calculated = Gamma / sqrt((2 * gamma / (gamma - 1)) * pressure_ratio_calculated ** (2 / gamma) * (
                    1 - pressure_ratio_calculated ** ((gamma - 1) / gamma)))
        error = abs(epsilon_true - epsilon_calculated)

        pressure_ratio_calculated = pressure_ratio_calculated + pressure_stepsize

    return pressure_ratio_calculated


def radius_function(a, b, c, d, z):
    return (a + (sqrt(b + 1000 * c * z) / d)) * 0.001


def performance_loss(c, dz, d):
    # Determine the performance lost due to radial thrust force component


    a = -1.3360
    b = 2.3935
    # d = 0.0454


    # Determine z_min (z_value for which the graphite - zirconium transition radius is reached). Note that the transition
    # point is currently fixed at 72.15/1000 by design. This is considered a fixed constraint, so only the curve past
    # transition can be changed.

    z_min = 0
    r = 0
    r_transition = 72.15 / 1000.

    # gradually increase z_min until the resulting radius matches the transition radius
    while abs(r - r_transition) > 0.01:
        r = radius_function(a, b, c, d, z_min)
        z_min = z_min + (0.01 / 1000.)

    z_max = z_min + dz

    # Determine the nozzle exit angle. CHANGE: currently the radius function is extrapolated all the way to the end.
    # Shouldn't the function change past the transition point? Otherwise what is the point of determining the transition
    # point...
    delta_z = 0.005
    slope = (radius_function(a, b, c, d, z_max) - radius_function(a, b, c, d, z_max - delta_z)) / delta_z
    exit_angle = atan(slope)

    loss_factor = 1 - (1 - cos(exit_angle)) / 2

    return loss_factor


def exhaust_velocity(R_e, gamma, c, d, dz):
    # Determine flow exit velocity based on geometric and flow conditions


    # Seems like too many parameters here are hard-coded
    throat_area = pi * (40 * 0.001) ** 2
    R = 640.91
    T_c = 2170

    Gamma = sqrt(gamma) * (2 / (gamma + 1)) ^ ((gamma + 1) / (2 * (gamma - 1)))

    expansion_ratio = (pi * R_e ^ 2) / throat_area
    p_ratio = pressure_ratio(expansion_ratio, gamma)

    U_e = sqrt(2 * (gamma * R * T_c) / (gamma - 1) * (1 - p_ratio ^ ((gamma - 1) / gamma))) * performance_loss(c, dz,
                                                                                                               d)
    return U_e


def nozzle_properties(c, d, dz):
    # Determine nozzle volumetric parameters


    # Initialize parameters
    # dz = 100 / 1000; # m

    a = -1.3360
    b = 2.3935
    # c is called from the function input
    # d = 0.0454;

    # Determine z_min(z_value for which the graphite - zirconium transition radius is reached)
    z_min = 0
    r = 0
    R_transition = 72.15 / 1000.

    while abs(r - R_transition) > 0.1 / 1000.:
        r = radius_function(a, b, c, d, z_min)
        z_min = z_min + 0.001 / 1000.

    z_max = z_min + dz # m

    # Determine results through integration
    z = sy.Symbol("z")

    V_zirconium = (pi * 1. / 1000. * (1. / 1000. * dz + 2 * sy.integrate(radius_function(a, b, c, d, z),
                                                                         (z, z_min, z_max)))) * 1000000  # cm3
    V_titanium = (pi * 4. / 1000. * (4. / 1000. * dz + 2 * sy.integrate(radius_function(a, b, c, d, z),
                                                                        (z, z_min, z_max)))) * 1000000  # cm3


    #V_zirconium = (pi * 1. / 1000. * (1. / 1000. * dz + 2 * integrate(R, z_min, z_max))) * 1000000 # cm3
    #V_titanium = (pi * 4. / 1000. * (4. / 1000. * dz + 2 * integrate(R, z_min, z_max))) * 1000000 # cm3

    R_e = radius_function(a, b, c, d, z_max)

    # Return Results
    return [V_titanium, V_zirconium, R_e]


def objective_function(x):
    # Objective function for optimization. Runs simulation based on input vector x to determine performance.
    # TODO: Assess input parameters for vector x (= Optimization parameters)


    gamma = 1.15
    m = 10  # kg / s
    altitude = 0  # m
    t = 0  # s
    dt = 0.1
    t_simulation = 400
    t_thrust = 22.5
    p_c = 15 * 10 ^ 5  # Pa
    A_t = pi * (40. / 1000.) ** 2
    c_d = 1.5
    v = 0  # m / s
    mass = 350  # kg

    # Determine nozzle parameters from input

    nozzle = nozzle_properties(x[1], x[2], x[3])  # nozzle_properties returns [V_titanium, V_zirconium, R_e].
    # TODO: Change function outputs to dictionaries for clarity
    A_e = (pi * nozzle[3] ^ 2)

    u_e = exhaust_velocity(nozzle[3], gamma, x(1), x(2), x(3))
    loss_factor = performance_loss(x[1], x[2], x[3])

    epsilon = A_e / A_t
    p_e = p_c * pressure_ratio(epsilon, gamma)

    data = []
    time = []

    # Start simulation loop

    while t < t_simulation:

        # Determine current atmospheric conditions
        atm = determine_atmosphere(altitude, atmospheric_data)
        p_a = atm['pressure']
        rho = atm['density']

        # Calculate net force and acceleration
        if (t < t_thrust):
            F_t = loss_factor * m * u_e + (p_e - p_a) * A_e

        else:
            F_t = 0

        F_d = c_d * 0.5 * rho * v ^ 2 * ((100. / 1000.) ** 2 * pi)
        F_w = mass * 9.80665

        F = F_t - F_d - F_w
        a = F / mass

        # Calculate distance travelled
        d = v * dt + 0.5 * a * dt ^ 2
        altitude = altitude + d

        if altitude < 0:
            altitude = 0

        # Update variables
        if t < t_thrust:
            mass = mass - m * dt

        t = t + dt
        v = v + a * dt
        time = time.append(t)
        data = data.append(altitude)

    # plot(time, data)
    # fprintf('performance evaluated.\n');
    performance = (1 / max(data)) * 1e6

    return performance


# =================================================== OPTIMIZATION =====================================================


# Main optimization code here

