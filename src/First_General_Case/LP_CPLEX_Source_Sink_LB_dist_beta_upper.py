from data import *
from read_data import *
from docplex.mp.model import Model
import json
import matplotlib.pyplot as plt


def create_variables(model, hits):
    layers = sorted(list(hits.keys()))
    K = len(layers) - 1

    v = []
    for p_1 in range(0, K):
        n_p_1 = len(hits[layers[p_1]]) + 1
        for p_2 in range(p_1 + 1, K + 1):
            n_p_2 = len(hits[layers[p_2]]) + 1
            for i in range(1, n_p_1):
                for j in range(1, n_p_2):
                    v.append((p_1, p_2, i, j))

    phi = model.binary_var_dict(v, name="phi")

    v = []
    # for p_1 in range(0, K - 1):
    #     n_p_1 = len(hits[layers[p_1]]) + 1
    #     for i in range(1, n_p_1):
    #         for p_2 in range(p_1 + 1, K):
    #             n_p_2 = len(hits[layers[p_2]]) + 1
    #             for j in range(1, n_p_2):
    #                 for p_3 in range(p_2 + 1, K + 1):
    #                     n_p_3 = len(hits[layers[p_3]]) + 1
    #                     for k in range(1, n_p_3):
    #                         v.append((p_1, p_2, p_3, i, j, k))
    for p_1 in range(0, K - 1):
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


def run(hits, M, alpha, gamma, NT, model_path_out, solution_path_out, figure_path_out):
    model = Model(name="Track")
    layers = sorted(list(hits.keys()))
    K = len(layers) - 1

    # create_variables
    ob, phi, c, nt, cp = create_variables(model, hits)
    nt = 6

    tracks = [[[1, 1], [2, 1], [3, 1], [4, 1], [5, 1], [6, 1], [7, 1]],
              [[1, 2], [2, 2], [3, 2], [4, 2], [5, 2], [6, 2], [7, 2]],
              [[1, 3], [2, 4], [3, 4], [4, 3], [5, 4], [6, 4], [7, 4]],
              [[1, 4], [2, 3], [3, 3], [4, 4], [5, 3], [6, 3], [7, 3]],
              [[1, 5], [2, 5], [3, 5], [4, 5], [5, 5], [6, 5], [7, 5]],
              [[1, 6], [2, 6], [3, 6], [4, 6], [5, 6], [6, 6], [7, 6]]]

    first_layer = [t[0] for t in tracks]
    for t in first_layer:
        p_2 = t[0]
        j = t[1]
        model.add_constraint(phi[0, p_2, 1, j] == 1)

    last_layer = [t[-1] for t in tracks]
    for t in first_layer:
        p_1 = t[0]
        i = t[1]
        model.add_constraint(phi[p_1, 8, i, 1] == 1)

    for track in tracks:
        for id in range(len(track)-1):
            p_1 = track[id][0]
            i = track[id][1]
            p_2 = track[id+1][0]
            j = track[id+1][1]
            model.add_constraint(phi[p_1, p_2, i, j] == 1)


    # n_p_2 = len(hits[layers[1]]) + 1
    # for j in range(3, 4):
    #     model.add_constraint(phi[0, 1, 1, j] == 1)
    #     # print('hit =', hits[layers[1]][1].hit_id)
    #     print('hit =', hits[layers[1]][j].hit_id)
    #     print(0, 1, 1, j)
    #
    # for p_1 in range(1, K-1):
    #     n_p_1 = len(hits[layers[p_1]]) + 1
    #     n_p_2 = len(hits[layers[p_1 + 1]]) + 1
    #     for i in range(1, n_p_1):
    #         for j in range(3, 4):
    #             if i == j:
    #                 model.add_constraint(phi[p_1, p_1 + 1, i, j] == 1)
    #                 print('hit =', hits[layers[p_1]][i].hit_id)
    #                 print('hit =', hits[layers[p_1 +1]][j].hit_id)
    #                 print(p_1, p_1 + 1, i, j)
    #
    # n_p_1 = len(hits[layers[7]]) + 1
    # for i in range(3, 4):
    #     model.add_constraint(phi[7, 8, i, 1] == 1)
    #     print('hit =', hits[layers[7]][i].hit_id)
    #     # print('hit =', hits[layers[8]][1].hit_id)
    #     print(7, 8, i, 1)
    # phi[1, 2, 3, 3] = 1
    # phi[2, 3, 3, 3] = 1
    # phi[3, 4, 3, 3] = 1
    # phi[4, 5, 3, 3] = 1
    # phi[1, 2, 4, 4] = 1
    # phi[2, 3, 4, 4] = 1
    # phi[3, 4, 4, 4] = 1
    # phi[4, 5, 4, 4] = 1
    # for k, v in phi.items():
    #     print(k, v)
    print("---First constraints---")
    tmp = 0
    cc_1 = 0
    for p_2 in range(1, K):
        n_p_2 = len(hits[layers[p_2]]) + 1
        for j in range(1, n_p_2):
            tmp += phi[0, p_2, 1, j]
    cc_1 += 1
    cn = "FC_" + str(cc_1)
    model.add_constraint(tmp == nt, ctname=cn)

    # Second constraints:
    print("---Second constraints---")
    tmp = 0
    cc_2 = 0
    for p_1 in range(1, K):
        n_p_1 = len(hits[layers[p_1]]) + 1
        for i in range(1, n_p_1):
            tmp += phi[p_1, K, i, 1]
    cc_2 += 1
    cn = "SC_" + str(cc_2)
    model.add_constraint(tmp == nt, ctname=cn)

    # Third constraints:
    print("---Third constraints---")
    print("---Fourth constraints---")
    cc_3 = 0
    cc_4 = 0
    for p_1 in range(1, K):
        n_p_1 = len(hits[layers[p_1]]) + 1
        for i in range(1, n_p_1):
            t_1 = 0
            for p_2 in range(0, p_1):
                n_p_2 = len(hits[layers[p_2]]) + 1
                for j in range(1, n_p_2):
                    t_1 += phi[p_2, p_1, j, i]

            t_2 = 0
            for p_2 in range(p_1 + 1, K + 1):
                n_p_2 = len(hits[layers[p_2]]) + 1
                for j in range(1, n_p_2):
                    t_2 += phi[p_1, p_2, i, j]
            cc_3 += 1
            cn = "TC_" + str(cc_3)
            model.add_constraint(t_1 == t_2, ctname=cn)

            cc_4 += 1
            cn = "FoC_" + str(cc_4)
            model.add_constraint(t_1 <= 1, ctname=cn)

    print("---Fifth constraints---")
    betas = []
    cc_5 = 0
    min_cost = 0
    for p_1 in range(0, K - 1):
        n_p_1 = len(hits[layers[p_1]]) + 1
        min_beta = 1000000000
        for i in range(1, n_p_1):
            for p_2 in range(p_1 + 1, K):
                n_p_2 = len(hits[layers[p_2]]) + 1
                for j in range(1, n_p_2):
                    for p_3 in range(p_2 + 1, K + 1):
                        n_p_3 = len(hits[layers[p_3]]) + 1
                        for k in range(1, n_p_3):
                            h_i = hits[layers[p_1]][i - 1]
                            h_j = hits[layers[p_2]][j - 1]
                            h_k = hits[layers[p_3]][k - 1]
                            seg_1 = Segment(h_j, h_i)
                            seg_2 = Segment(h_j, h_k)
                            angle = Angle(seg_1=seg_1, seg_2=seg_2).angle
                            dist = distance(h_i, h_j) // 100 + distance(h_j, h_k) // 100
                            beta = angle * dist
                            betas.append(beta)

                            if beta < min_beta:
                                min_beta = beta
                            cc_5 += 1
                            cn = "FiC_" + str(cc_5)
                            model.add_constraint(
                                beta <= c[p_1, i] + ((2 - phi[p_1, p_2, i, j] - phi[p_2, p_3, j, k]) * M),
                                ctname=cn)
        min_cost += min_beta

    beta_upper = sum(betas) / len(betas)
    print("beta upper =", beta_upper)

    cp_tmp = 0
    # for p_1 in range(0, K - 1):
    #     n_p_1 = len(hits[layers[p_1]]) + 1
    #     for i in range(1, n_p_1):
    #         for p_2 in range(p_1 + 1, K):
    #             n_p_2 = len(hits[layers[p_2]]) + 1
    #             for j in range(1, n_p_2):
    #                 for p_3 in range(p_2 + 1, K + 1):
    #                     n_p_3 = len(hits[layers[p_3]]) + 1
    #                     for k in range(1, n_p_3):
    #                         if p_1 != 0 and p_1 != K - 2:
    #                             model.add_constraint(c[p_1, p_2, p_3, i, j, k] <= beta_upper)
    #                         cp_tmp += c[p_1, p_2, p_3, i, j, k]
    for p_1 in range(0, K - 1):
        n_p_1 = len(hits[layers[p_1]]) + 1
        for i in range(1, n_p_1):
            # if p_1 != 0 and p_1 != K - 2:
            #     model.add_constraint(c[p_1, i] <= beta_upper)
            cp_tmp += c[p_1, i]
    total_hits = sum(len(hits[layers[p]]) for p in range(1, K))
    print("total_hits:", total_hits)

    # objective = alpha * (total_hits - nt) + gamma * cp
    objective = cp

    # model.add_constraint(cp >= min_cost * nt)
    model.add_constraint(cp >= cp_tmp)

    model.add_constraint(ob >= objective)
    model.set_objective('min', ob)

    model.print_information()
    model.export_as_lp(model_path_out)
    model.solve(log_output=True)
    print(model.solve_status)
    model.solution.export(solution_path_out)
    segments = []
    with open(solution_path_out) as f:
        result = json.load(f)['CPLEXSolution']['variables']
        for var in result:
            print(var)
            phi_p_p_i_j = var['name'].split('_')
            if "phi" in phi_p_p_i_j[0]:
                p_1 = int(phi_p_p_i_j[1])
                p_2 = int(phi_p_p_i_j[2])
                i = int(phi_p_p_i_j[3])
                j = int(phi_p_p_i_j[4])
                h_1 = hits[layers[p_1]][i - 1]
                h_2 = hits[layers[p_2]][j - 1]
                segments.append([p_1, p_2, h_1, h_2])
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
        p_1 = segment[0]
        p_2 = segment[1]
        h1 = segment[2]
        h2 = segment[3]
        if p_1 == 0:
            ax.plot(xs=[h1.x, h2.x], ys=[h1.y, h2.y], zs=[h1.z, h2.z], color='green')
        elif p_2 == 8:
            ax.plot(xs=[h1.x, h2.x], ys=[h1.y, h2.y], zs=[h1.z, h2.z], color='black')
        else:
            ax.plot(xs=[h1.x, h2.x], ys=[h1.y, h2.y], zs=[h1.z, h2.z], color='blue')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    plt.savefig(out)
    plt.show()


def create_source_sink(hits):
    layers = sorted(list(hits.keys()))
    first_layer = hits[layers[0]]
    last_layer = hits[layers[-1]]

    source = Hit(
        hit_id=10,
        x=sum([h.x for h in first_layer]) / len(first_layer),
        y=sum([h.y for h in first_layer]) / len(first_layer),
        z=sum([h.z for h in first_layer]) / len(first_layer) - 200,
        volume_id=0,
        layer_id=0,
        module_id=0
    )

    sink = Hit(
        hit_id=10,
        x=sum([h.x for h in last_layer]) / len(last_layer),
        y=sum([h.y for h in last_layer]) / len(last_layer),
        z=sum([h.z for h in last_layer]) / len(last_layer) + 200,
        volume_id=0,
        layer_id=0,
        module_id=0
    )

    return source, sink


if __name__ == '__main__':
    src_path = '../data_selected'
    data_path = src_path + '/6hits/unknown_track/hits.csv'
    hits = read_hits(data_path)
    source, sink = create_source_sink(hits)
    hits[0] = [source]
    hits[16] = [sink]

    model_path_out = "results/6hits/unknown_track/model_docplex_LB_dist_ss_test.lp"
    solution_path_out = "results/6hits/unknown_track/solution_LB_dist_ss_test.json"
    figure_path_out = "results/6hits/unknown_track/result_LB_dist_ss_test.PNG"

    M = 10000
    alpha = 10000
    gamma = 1
    NT = 6
    run(hits, M, alpha, gamma, NT, model_path_out, solution_path_out, figure_path_out)
