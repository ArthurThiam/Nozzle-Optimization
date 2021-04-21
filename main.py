# Define atmospheric data
atmospheric_data = {'base_altitude': [0, 11000, 20000, 32000, 47000, 51000, 71000],
                    'base_temperature': [288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65],
                    'base_density': [1.225, 0.36391, 0.08803, 0.01322, 0.00143, 0.00086, 0.000064],
                    'base_pressure': [101325, 22632.1, 5474, 868.02, 110.91, 66.94, 3.96],
                    'lapse': [-0.0065, 0, 0.001, 0.0028, 0, -0.0028, -0.002]}


def atmosphere(altitude, atmosphere):

    for i in range(len(atmosphere['base_altitude'])):
        if atmosphere['base_altitude'][i] <= altitude < atmosphere['base_altitude'][i+1]:

            base_altitude = atmosphere['base_altitude'][i]
            base_temperature = atmosphere['base_temperature'][i]
            base_density = atmosphere['base_density'][i]
            base_pressure = atmosphere['base_pressure'][i]
            lapse = atmosphere['lapse'][i]

    data = [base_altitude, base_temperature, base_density, base_pressure, lapse]
    print(data)
    return data


atmosphere(71000, atmospheric_data)
# fix edge case of altitude>71000



