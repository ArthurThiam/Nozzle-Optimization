from Nozzle import *
import plotly.express as px
import plotly.graph_objects as go

# ============================================ INITIALIZATION ==========================================================
solution_space = initialize(geometry, GA_settings[3])
gene_set = []

for solution in solution_space:
    solution_geometry = build_geometry(solution)
    gene_set.append(Chromosome(solution_geometry, engine_properties))

population = Population(gene_set)


# ============================================== MAIN CODE =============================================================
generation_list = [population]
generation_count = 0
generation_target = GA_settings[4]

counter = []
performance = []
average_performance = []

while generation_count < generation_target:

    population.select_survivors()                               # Determine survivors of current generation
    average = sum(population.evaluation)/len(population.evaluation)

    elite_chromosome = population.population[0]
    for chromosome in population.population:
        if chromosome.apogee > elite_chromosome.apogee:
            elite_chromosome = chromosome

    best = elite_chromosome.apogee
    best_genes = elite_chromosome.genes

    print('Population', generation_count, ': ', population.evaluation)
    print('Best chromosome: ', best_genes)
    print('Exit radius: ', elite_chromosome.R_e)
    print('Expansion ratio: ', elite_chromosome.epsilon)
    print('')
    print('Best performance: ', best, ' m')
    print('Average performance: ', average, ' m')
    print('----------------------------------------------------')

    performance.append(float(best))
    average_performance.append(float(average))
    counter.append(generation_count)

    population = Population(population.population)              # Generate new population with surviving members
    #generation_list.append([population])                          # Add new population to list of generations



    generation_count += 1

fig = go.Figure(data=go.Scatter(x=counter, y=performance))
fig.show()

fig = go.Figure(data=go.Scatter(x=counter, y=average_performance))
fig.show()



