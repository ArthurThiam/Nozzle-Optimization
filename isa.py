from math import *




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

        elif altitude > atmosphere['base_altitude'][7]:
            lapse = 0
            base_altitude = 100000
            base_temperature = 214.65
            base_density = 0
            base_pressure = 0

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