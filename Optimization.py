from ObjectiveFunctions import *
import configparser
import plotly.express as px
import random


# Read initial parameters
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
                     config.getfloat('Engine Properties', 'A_t')]

# Create genome instance
chromosome_1 = Chromosome(geometry, engine_properties)

# Run simulation test
simulation = objective_function_2(chromosome_1)
print(simulation['apogee'])
