import configparser

config = configparser.ConfigParser()
config.read('settings.ini')

# geometry = [config.getfloat('Geometry', 'a'),
#             config.getfloat('Geometry', 'b'),
#             config.getfloat('Geometry', 'c'),
#             config.getfloat('Geometry', 'd'),
#             config.getfloat('Geometry', 'dz'),
#             config.getfloat('Geometry', 'D_t')]

geometry_keys = ['a', 'b', 'c', 'd', 'dz', 'D_t']

geometry = {'a': config.getfloat('Geometry', 'a'),
            'b': config.getfloat('Geometry', 'b'),
            'c': config.getfloat('Geometry', 'c'),
            'd': config.getfloat('Geometry', 'd'),
            'dz': config.getfloat('Geometry', 'dz'),
            'D_t': config.getfloat('Geometry', 'D_t')}

design_constraints = {'Re_max': config.getfloat('Design Constraints', 'Re_max'),
                      'Rt_max': config.getfloat('Design Constraints', 'Rt_max'),
                      'dz_max': config.getfloat('Design Constraints', 'dz_max'),
                      'epsilon_min': config.getfloat('Design Constraints', 'epsilon_min')}

engine_properties = [config.getfloat('Engine Properties', 'gamma'),
                     config.getfloat('Engine Properties', 'm'),
                     config.getfloat('Engine Properties', 't_thrust'),
                     config.getfloat('Engine Properties', 'A_t'),
                     config.getfloat('Engine Properties', 'p_c')]

simulation_settings = [config.getfloat('Simulation Settings', 'dt'),
                       config.getfloat('Simulation Settings', 'c_d')]

GA_settings = [config.getfloat('GA Settings 1', 'population_steadiness'),
               config.getfloat('GA Settings 1', 'crossover_prob'),
               config.getfloat('GA Settings 1', 'mutation_prob'),
               config.getint('GA Settings 1', 'population_size'),
               config.getint('GA Settings 1', 'generation_target'),
               config.getfloat('GA Settings 1', 'mutation_percentage'),
               config.getint('GA Settings 1', 'unstable_genes'),
               config.getint('GA Settings 1', 'random_solutions')]
