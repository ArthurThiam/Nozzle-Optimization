import configparser

config = configparser.ConfigParser()
config.read('settings.ini')

geometry = [config.getfloat('Geometry', 'a'),
            config.getfloat('Geometry', 'b'),
            config.getfloat('Geometry', 'c'),
            config.getfloat('Geometry', 'd'),
            config.getfloat('Geometry', 'dz')]

engine_properties = [config.getfloat('Engine Properties', 'gamma'),
                     config.getfloat('Engine Properties', 'm'),
                     config.getfloat('Engine Properties', 't_thrust'),
                     config.getfloat('Engine Properties', 'A_t'),
                     config.getfloat('Engine Properties', 'p_c')]

simulation_settings = [config.getfloat('Simulation Settings', 'dt'),
                       config.getfloat('Simulation Settings', 'c_d')]