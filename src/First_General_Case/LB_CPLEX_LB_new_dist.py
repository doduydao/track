from data import *
from read_data import *
from docplex.mp.model import Model
import json
import matplotlib.pyplot as plt


def create_variables(model, hits):
    layers = sorted(list(hits.keys()))
    K = len(layers) + 1
    layers = [0] + layers + [0]
    v = []
    for p_1 in range(1, K - 1):
        n_p_1 = len(hits[layers[p_1]]) + 1
        for p_2 in range(p_1 + 1, K):
            n_p_2 = len(hits[layers[p_2]]) + 1
            for i in range(1, n_p_1):
                for j in range(1, n_p_2):
                    v.append((p_1, p_2, i, j))

    phi = model.binary_var_dict(v, name="phi")

    v = []
    for p_1 in range(1, K - 2):
        n_p_1 = len(hits[layers[p_1]]) + 1
        for i in range(1, n_p_1):
            v.append((p_1, i))

    c = model.continuous_var_dict(v, name="c", lb=0)
    ob = model.continuous_var(name="ob")
    nt = model.continuous_var(name="nt")
    cp = model.continuous_var(name="cp")
    return ob, phi, c, nt, cp


def distance(h1, h2):
    distance = math.sqrt((h2.x - h1.x) ** 2 + (h2.y - h1.y) ** 2 + (h2.z - h1.z) ** 2)
    return distance


def run(hits, M, alpha, gamma, model_path_out, solution_path_out, figure_path_out):
    model = Model(name="Track")
    layers = sorted(list(hits.keys()))
    K = len(layers) + 1
    layers = [0] + layers + [0]

    # create_variables
    ob, phi, c, nt, cp = create_variables(model, hits)

    # First constraints:
    print("---First constraints---")
    cc_1 = 0
    n_p_1 = len(hits[layers[1]]) + 1
    for i in range(1, n_p_1):
        tmp = 0
        for p_2 in range(2, K - 1):
            n_p_2 = len(hits[layers[p_2]]) + 1
            for j in range(1, n_p_2):
                tmp += phi[1, p_2, i, j]
        cc_1 += 1
        cn = "FC_" + str(cc_1)
        model.add_constraint(tmp <= 1, ctname=cn)

    print("---Second constraints---")
    cc_2 = 1
    tmp = 0
    n_p_1 = len(hits[layers[1]]) + 1
    for i in range(1, n_p_1):
        n_p_2 = len(hits[layers[K - 1]]) + 1
        for j in range(1, n_p_2):
            tmp += phi[1, K - 1, i, j]
    cc_2 += 1
    cn = "SC_" + str(cc_2)
    model.add_constraint(tmp == 0, ctname=cn)

    print("---Thirth constraints---")
    cc_3 = 0
    for p_1 in range(2, K - 1):
        n_p_1 = len(hits[layers[p_1]]) + 1
        for i in range(1, n_p_1):
            tmp = 0
            for p_2 in range(p_1 + 1, K):
                n_p_2 = len(hits[layers[p_2]]) + 1
                for j in range(1, n_p_2):
                    tmp += phi[p_1, p_2, i, j]
            cc_3 += 1
            cn = "TC_" + str(cc_3)
            model.add_constraint(tmp <= 1, ctname=cn)

    print("---Fourth constraints---")
    cc_4 = 0
    n_p_2 = len(hits[layers[K - 1]]) + 1
    for j in range(1, n_p_2):
        tmp = 0
        for p_1 in range(2, K - 1):
            n_p_1 = len(hits[layers[p_1]]) + 1
            for i in range(1, n_p_1):
                tmp += phi[p_1, K - 1, i, j]
        cc_4 += 1
        cn = "FoC_" + str(cc_4)
        model.add_constraint(tmp <= 1, ctname=cn)

    print("---Fifth constraints---")
    cc_5 = 0
    for p_2 in range(2, K - 1):
        n_p_2 = len(hits[layers[p_2]]) + 1
        for j in range(1, n_p_2):
            tmp = 0
            for p_1 in range(1, p_2):
                n_p_1 = len(hits[layers[p_1]]) + 1
                for i in range(1, n_p_1):
                    tmp += phi[p_1, p_2, i, j]
            cc_5 += 1
            cn = "FiC_" + str(cc_5)
            model.add_constraint(tmp <= 1, ctname=cn)

    print("---Sixth constraints---")
    cc_6 = 0
    for p_1 in range(2, K - 1):
        n_p_1 = len(hits[layers[p_1]]) + 1
        for i in range(1, n_p_1):
            t_1 = 0
            for p_2 in range(1, p_1):
                n_p_2 = len(hits[layers[p_2]]) + 1
                for j in range(1, n_p_2):
                    t_1 += phi[p_2, p_1, j, i]
            t_2 = 0
            for p_2 in range(p_1 + 1, K):
                n_p_2 = len(hits[layers[p_2]]) + 1
                for j in range(1, n_p_2):
                    t_2 += phi[p_1, p_2, i, j]
            cc_6 += 1
            cn = "SiC_" + str(cc_6)
            model.add_constraint(t_1 == t_2, ctname=cn)

    print("---Seventh constraint---")
    tmp = 0
    n_p_1 = len(hits[layers[1]]) + 1
    for i in range(1, n_p_1):
        for p_2 in range(2, K - 1):
            n_p_2 = len(hits[layers[p_2]]) + 1
            for j in range(1, n_p_2):
                tmp += phi[1, p_2, i, j]
    cn = "SeC_" + str(1)
    model.add_constraint(tmp == nt, ctname=cn)

    print("---Eighth constraint---")
    cn = "EiC_" + str(1)
    model.add_constraint(nt <= len(hits[layers[1]]), ctname=cn)

    print("---Nineth constraints---")
    cc_9 = 0
    betas = []
    cost_lb = 0
    for p_1 in range(1, K - 2):
        n_p_1 = len(hits[layers[p_1]]) + 1
        beta_lower = 1000000
        for i in range(1, n_p_1):
            for p_2 in range(p_1 + 1, K - 1):
                n_p_2 = len(hits[layers[p_2]]) + 1
                for j in range(1, n_p_2):
                    for p_3 in range(p_2 + 1, K):
                        n_p_3 = len(hits[layers[p_3]]) + 1
                        for k in range(1, n_p_3):
                            h_i = hits[layers[p_1]][i - 1]
                            h_j = hits[layers[p_2]][j - 1]
                            h_k = hits[layers[p_3]][k - 1]
                            seg_1 = Segment(h_j, h_i)
                            seg_2 = Segment(h_j, h_k)

                            angle = Angle(seg_1=seg_1, seg_2=seg_2).angle
                            dist = (distance(h_i, h_j) + distance(h_j, h_k))
                            beta = angle * dist
                            betas.append(beta)

                            if beta < beta_lower:
                                beta_lower = beta

                            cc_9 += 1
                            cn = "NiC_" + str(cc_9)
                            model.add_constraint(
                                beta <= c[p_1, i] + (2 - phi[p_1, p_2, i, j] - phi[p_2, p_3, j, k]) * M, ctname=cn)
        cost_lb += beta_lower

    beta_max = max(betas)
    beta_min = min(betas)
    beta_average = sum(betas) / len(betas)
    beta_upper = (beta_average + beta_min) / 2
    print("Max beta x dist = ", beta_max)
    print("Min beta x dist = ", beta_min)
    print("Mean beta x dist = ", beta_average)
    print("Limit beta x dist = ", beta_average)

    cp_tmp = 0
    for p_1 in range(1, K - 2):
        n_p_1 = len(hits[layers[p_1]]) + 1
        for i in range(1, n_p_1):
            cp_tmp += c[p_1, i]
            model.add_constraint(c[p_1, i] <= beta_upper)

    total_hits = sum(len(hits[layers[p]]) for p in range(1, K))
    print("total_hits:", total_hits)
    objective = alpha * (len(hits[layers[1]]) - nt) + gamma * cp

    print("Add lower bound")
    model.add_constraint(cp >= cost_lb * nt)
    # model.add_constraint(cp_tmp >= beta_min * nt)

    # print("Add upper bound")
    # model.add_constraint(cp_tmp <= cost_ub * nt)
    # limit_upper_cost = nt * (beta_min + (beta_upper - beta_min) / total_hits)
    # model.add_constraint(cp_tmp <= limit_upper_cost)

    model.add_constraint(cp >= cp_tmp)
    model.add_constraint(ob >= objective)
    model.set_objective('min', ob)

    model.print_information()
    model.export_as_lp(model_path_out)
    model.solve(log_output=True)
    print(model.solve_status)
    model.solution.export(solution_path_out)

    f = open(solution_path_out)
    result = json.load(f)
    f.close()

    result = result['CPLEXSolution']['variables']

    segments = []
    for var in result:
        print(var)
        var_name = var['name']
        var_value = round(float(var['value']))
        phi_p_p_i_j = var_name.split('_')
        if 'nt' in phi_p_p_i_j[0] or 'c' in phi_p_p_i_j[0] or var_value != 1.0 or 's' in phi_p_p_i_j[0] or 'q' in \
                phi_p_p_i_j[0] or 'ob' in \
                phi_p_p_i_j[0]:
            continue

        p_1 = int(phi_p_p_i_j[1])
        p_2 = int(phi_p_p_i_j[2])
        i = int(phi_p_p_i_j[3])
        j = int(phi_p_p_i_j[4])
        h_1 = hits[layers[p_1]][i - 1]
        h_2 = hits[layers[p_2]][j - 1]
        segments.append([h_1, h_2])
    display(hits, segments, figure_path_out)


def display(hits, segments, out=""):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    layers = list(hits.keys())
    xs = []
    ys = []
    zs = []

    for p in layers:
        h_p = hits[p]
        for h in h_p:
            xs.append(h.x)
            ys.append(h.y)
            zs.append(h.z)
    ax.scatter(xs, ys, zs, marker='o', color='red')

    for segment in segments:
        h1 = segment[0]
        h2 = segment[1]
        ax.plot(xs=[h1.x, h2.x], ys=[h1.y, h2.y], zs=[h1.z, h2.z], color='blue')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    plt.savefig(out)
    plt.show()


if __name__ == '__main__':
    src_path = '../data_selected'
    data_path = src_path + '/15hits/unknown_track/hits.csv'
    model_path_out = "results/15hits/unknown_track/model_docplex_LB_new_dist_test.lp"
    solution_path_out = "results/15hits/unknown_track/solution_new_dist_test.json"
    figure_path_out = "results/15hits/unknown_track/result_new_dist_test.PNG"

    M = 100000
    alpha = 10000
    gamma = 1
    hits = read_hits(data_path)
    result = run(hits, M, alpha, gamma, model_path_out, solution_path_out, figure_path_out)
