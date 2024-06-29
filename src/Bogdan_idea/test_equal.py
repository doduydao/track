from A_QUBM import *
from A_QUBM_parallel import *



if __name__ == '__main__':
    src_path = '../../src/data_selected'
    folder = '/2hits/known_track/'
    out_path = '/Users/doduydao/daodd/PycharmProjects/track/src/Bogdan_idea/results'
    data_path = src_path + folder + 'hits.csv'
    figure_path_out = out_path + folder + "result_D_AQUBM.PNG"
    check_path(out_path + folder)

    # read data
    print("Loading data...")
    start = time.time()
    hits_by_layers, hits = read_hits(data_path)
    end = time.time()
    print("Loaded data! Execution time: ", end - start)
    # weight
    A = 2
    beta_max = math.pi / A
    m = 1

    re_calculate = False
    costs_path_out = out_path + folder + "pi_" + str(A) + "costs.json"
    if os.path.exists(costs_path_out) == False:
        re_calculate = True
    if re_calculate:
        # calculate costs
        print("\n----Compute cost----")
        start = time.time()
        costs = get_costs(hits_by_layers, beta_max)
        end = time.time()
        print('Complete!. Execution time: ', end - start, 's')

        print("\n---Write cost out---")
        print("Path: ", costs_path_out)
        write_costs(costs, costs_path_out, m)

    print("---Load cost---")
    costs = load_costs(costs_path_out)

    print("\n ----Create Hamiltonian----")
    model_path_out = out_path + folder + "model_A_QUBM.lp"
    alpha = 100
    gamma = 1


    H_1 = create_hamiltonian(costs, hits, hits_by_layers, alpha, gamma, model_path_out)
    H_2 = create_objective_function(costs, hits, hits_by_layers, alpha, gamma)

    print(H_1 == H_2)

    print("----Created Hamiltonian----")

    cal_expected_value(hits, m)
    print("\n----Simulating H_1 ----")
    num_reads = 100
    ob_value, result = simulate_annealing(H_1, num_reads)
    print("Objective value:", ob_value)
    display(hits, result, out=figure_path_out)

    print("\n----Simulating H_2----")
    ob_value, result = simulate_annealing(H_2, num_reads)
    print("Objective value:", ob_value)
    display(hits, result, out=figure_path_out)
