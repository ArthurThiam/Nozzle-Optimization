from Nozzle import *
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots


# GA Function
def nozzle_ga():
    print('\nRunning GA.')
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
    elites = []

    while generation_count < generation_target:

        population.select_survivors()                               # Determine survivors of current generation
        average = sum(population.evaluation)/len(population.evaluation)

        elite_chromosome = population.population[0]
        for chromosome in population.population:
            if chromosome.apogee > elite_chromosome.apogee:
                elite_chromosome = chromosome

        best = elite_chromosome.apogee

        # print('Population', generation_count, ': ', population.evaluation)
        # print('Best chromosome: ', best_genes)
        # print('Exit radius: ', elite_chromosome.R_e)
        # print('Expansion ratio: ', elite_chromosome.epsilon)
        # print('')
        # print('Best performance: ', best, ' m')
        # print('Average performance: ', average, ' m')
        # print('----------------------------------------------------')

        elites.append(elite_chromosome)
        performance.append(float(best))
        average_performance.append(float(average))
        counter.append(generation_count)

        population = Population(population.population)              # Generate new population with surviving members
        generation_count += 1

    # fig = go.Figure()
    # fig.add_trace(go.Scatter(x=counter, y=performance, name='Best Performer'))
    # fig.add_trace(go.Scatter(x=counter, y=average_performance, name='Average Performance'))
    #
    # fig.update_layout(title='GA Nozzle Optimization for Flight Apogee',
    #                    xaxis_title='Iteration',
    #                    yaxis_title='Apogee [m]')
    #
    #fig.show()

    print('GA complete.')

    return [elites[-1], performance[-1]]


# GA Solution instance - runs nozzle_ga 10 times
def ga_run(settings_chapter):
    run_count = 0
    counter = []
    performance_tracker = []
    settings_log = []
    solutions = []
    runtimes = []

    while run_count < 10:
        begin = time.time()
        ga = nozzle_ga()
        end = time.time()

        performance_tracker.append(ga[1])
        solutions.append(ga[0].genes)
        settings_log.append(settings_chapter)

        runtimes.append(end-begin)
        counter.append(run_count)
        run_count += 1

    # fig = go.Figure()
    # fig.add_trace(go.Scatter(x=counter, y=performance_tracker, name='Best Performance'))
    #
    # fig.update_layout(title='Best performance variation',
    #                        xaxis_title='Iteration',
    #                        yaxis_title='Apogee [m]')
    #
    # fig.show()
    avg = sum(runtimes)/len(runtimes)

    print('=========================================')
    print('\nAverage performance: ', sum(performance_tracker)/len(performance_tracker))
    print('Average generation run time: ', avg)
    #print('\nSolutions: ', solutions)

    return performance_tracker, settings_log, avg


# Generate stats for various settings profiles
def ga_stats():
    settings_selector = 1
    x = []
    time_x = []
    time_y = []
    fig = make_subplots(specs=[[{"secondary_y": True}]])

    while settings_selector <= 5:
        settings_chapter = 'GA Settings ' + str(settings_selector)

        GA_settings = [config.getfloat(settings_chapter, 'population_steadiness'),
                       config.getfloat(settings_chapter, 'crossover_prob'),
                       config.getfloat(settings_chapter, 'mutation_prob'),
                       config.getint(settings_chapter, 'population_size'),
                       config.getint(settings_chapter, 'generation_target'),
                       config.getfloat(settings_chapter, 'mutation_percentage'),
                       config.getint(settings_chapter, 'unstable_genes'),
                       config.getint(settings_chapter, 'random_solutions')]

        ga_run_data = ga_run(settings_chapter)
        fig.add_trace(go.Box(x=ga_run_data[1], y=ga_run_data[0]), secondary_y=False)

        time_y.append(ga_run_data[2])
        time_x.append([ga_run_data[1][0]])


        settings_selector += 1

    fig.add_trace(go.Scatter(x=['GA Settings 1', 'GA Settings 2', 'GA Settings 3', 'GA Settings 4', 'GA Settings 5'], y=time_y), secondary_y=True)
    fig.update_yaxes(title_text="GA Performance [m]", secondary_y=False)
    fig.update_yaxes(title_text="Average Computation Time per Generation [s]", secondary_y=True)
    fig.update_xaxes(title_text="Settings")

    fig.show()


# ========================================= MAIN CODE =======================================

# generate one solution, print its genes, exit radius, expansion ratio and performance
solution = nozzle_ga()
print(solution[0].genes)
print(solution[0].R_e)
print(solution[0].epsilon)
print(solution[1])

# If you want to run the stats for all parameter sets
# ga_stats()



