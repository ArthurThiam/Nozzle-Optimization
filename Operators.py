import random
from Nozzle import Chromosome
from Settings_import import *

# This file contains all operators required by the Population class


# The single point cross-over takes two chromosome instances and returns their two offspring
def single_point_crossover(chromosome_1, chromosome_2):

    # pull parent genes from input chromosomes
    parent_1 = chromosome_1.genes()
    parent_2 = chromosome_2.genes()

    # determine cross over point
    cross_over_point = random.randint(0, len(parent_1))

    # perform cross over
    for gene in range(cross_over_point, len(parent_1)):
        parent_1[gene], parent_2[gene] = parent_2[gene], parent_1[gene]

    # generate offspring chromosomes
    offspring = [Chromosome(parent_1, engine_properties), Chromosome(parent_2, engine_properties)]

    return offspring

# The survivor selection choses what population members get replaced with offspring and which ones continue to the next
# generation


def survivor_selection():
    return 0

