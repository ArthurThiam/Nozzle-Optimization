from math import *

# Define atmospheric data
atmospheric_data = {'base_altitude': [0, 11000, 20000, 32000, 47000, 51000, 71000, 100000],
                    'base_temperature': [288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65],
                    'base_density': [1.225, 0.36391, 0.08803, 0.01322, 0.00143, 0.00086, 0.000064],
                    'base_pressure': [101325, 22632.1, 5474, 868.02, 110.91, 66.94, 3.96],
                    'lapse': [-0.0065, 0, 0.001, 0.0028, 0, -0.0028, -0.002]}

# Determine pressure and density at the current altitude. Works up to and excluding 100,000 m
def determine_atmosphere(altitude, atmosphere):

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
            print('case one i:', i)

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


def pressure_ratio(expansion_ratio, gamma):
    return 0


def performance_loss(c, dz, d):
    return 0


def exhaust_velocity(R_e, gamma, c, d, dz):

    throat_area = pi * (40 * 0.001) ** 2
    R = 640.91
    T_c = 2170

    Gamma = sqrt(gamma) * (2 / (gamma + 1)) ^ ((gamma + 1) / (2 * (gamma - 1)))

    expansion_ratio = (pi * R_e ^ 2) / throat_area
    p_ratio = pressure_ratio(expansion_ratio, gamma)

    U_e = sqrt(2 * (gamma * R * T_c) / (gamma - 1) * (1 - (p_ratio) ^ ((gamma - 1) / gamma))) * performance_loss(c, dz,
                                                                                                                 d)

    return U_e