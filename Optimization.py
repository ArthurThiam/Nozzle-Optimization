from Nozzle import Genome
import sympy as sy
import plotly.express as px
import random
from math import *


# Define initial properties
a = -1.3360
b = 2.3935
c = 0.0888  # 0.0909        [-] c value for nozzle curve definition
d = 0.05  # 0.0454
dz = 0.1015  # 0.1          [m]   skirt length

geometry = [a, b, c, d, dz]

gamma = 1.15
m = 10 # kg/s
t_thrust = 22.5 # s
A_t = pi * (40. / 1000.) ** 2

engine_properties = [gamma, m, t_thrust, A_t]

# Create genome instance
genome_1 = Genome(geometry, engine_properties)

print(genome_1.exit_radius())