from Nozzle import *
import matplotlib

# ============================================ INITIALIZATION ==========================================================
solution_space = initialize(geometry, GA_settings[3])
gene_set = []

for solution in solution_space:
    gene_set.append(Chromosome(solution, engine_properties))

population = Population(gene_set)


# ============================================== MAIN CODE =============================================================
generation_list = [population]
generation_count = 0
generation_target = GA_settings[4]

counter = []
performance = []

while generation_count < generation_target:

    population.select_survivors()                               # Determine survivors of current generation
    print('Population', generation_count, ': ', population.evaluation)
    print('Best performance: ', max(population.evaluation), ' m')
    print('----------------------------------------------------')

    performance.append(int(max(population.evaluation)))
    counter.append(generation_count)

    population = Population(population.population)              # Generate new population with surviving members
    generation_list.append(population)                          # Add new population to list of generations

    generation_count += 1




