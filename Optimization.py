from Nozzle import *

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

while generation_count < generation_target:

    population.select_survivors()                               # Determine survivors of current generation
    print('Best performance: ', max(population.evaluation), ' m')
    print('Population', generation_count, ': ', population.evaluation)
    print('----------------------------------------------------')
    population = Population(population.population)              # Generate new population with surviving members
    generation_list.append(population)                          # Add new population to list of generations


    generation_count += 1

# ===============================================================================

# start_time = time.time()
# runtime = time.time() - start_time

# print('New population: ', len(population.population))
# print('')
# print('total runtime [s]: ', runtime)
# print('average runtime per chromosome [s/chromosome]:', runtime/GA_settings[3])

