from data import *
from read_data import *
from docplex.mp.model import Model
import json
import matplotlib.pyplot as plt


def create_variables(model, hits):
    layers = sorted(list(hits.keys()))
    K = len(layers) - 1
    print(K)

    v = []
    for p_1 in range(0, K):
        for p_2 in range(1, K + 1):
            n_p_1 = len(hits[layers[p_1]]) + 1
            n_p_2 = len(hits[layers[p_2]]) + 1
            for i in range(1, n_p_1):
                for j in range(1, n_p_2):
                    v.append((p_1, p_2, i, j))
    phi = model.binary_var_dict(v, name="phi")

    v = []
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


def run(hits, M, alpha, gamma, model_path_out, solution_path_out, figure_path_out, NT):
    model = Model(name="Track")
    layers = sorted(list(hits.keys()))
    K = len(layers) - 1
    # create_variables
    ob, phi, c, nt, cp = create_variables(model, hits)
    nt = NT
    # first constraints:
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
    cc_3 = 0
    for p_1 in range(1, K):
        n_p_1 = len(hits[layers[p_1]]) + 1
        tmp_1 = 0
        tmp_2 = 0
        for i in range(1, n_p_1):
            t_1 = 0
            for p_2 in range(0, p_1):
                n_p_2 = len(hits[layers[p_2]]) + 1
                for j in range(1, n_p_2):
                    tmp_1 += phi[p_2, p_1, j, i]
                    t_1 += phi[p_2, p_1, j, i]
            t_2 = 0
            for p_2 in range(p_1 + 1, K + 1):
                n_p_2 = len(hits[layers[p_2]]) + 1
                for j in range(1, n_p_2):
                    tmp_2 += phi[p_1, p_2, i, j]
                    t_2 += phi[p_1, p_2, i, j]
            cc_3 += 1
            cn = "TC_" + str(cc_3)
            model.add_constraint(t_1 == t_2, ctname=cn)

    print("---Fourth constraints---")
    cc_4 = 0
    for p_1 in range(1, K):
        n_p_1 = len(hits[layers[p_1]]) + 1
        for i in range(1, n_p_1):
            tmp = 0
            for p_2 in range(0, p_1):
                n_p_2 = len(hits[layers[p_2]]) + 1
                for j in range(1, n_p_2):
                    tmp += phi[p_2, p_1, j, i]
            cc_4 += 1
            cn = "TC_" + str(cc_4)
            model.add_constraint(tmp <= 1, ctname=cn)

    print("---Fifth constraints---")
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
                            beta = Angle(seg_1=seg_1, seg_2=seg_2).angle * (
                                    distance(h_i, h_j) // 100 + distance(h_j, h_k) // 100)
                            if beta < min_beta:
                                min_beta = beta
                            cc_5 += 1
                            cn = "FiC_" + str(cc_5)
                            if p_1 == 0:
                                c[p_1, i] = 0
                            else:

                                model.add_constraint(
                                beta <= c[p_1, i] + (2 - phi[p_1, p_2, i, j] - phi[p_2, p_3, j, k]) * M, ctname=cn)
        min_cost += min_beta
    min_cost = min_cost * nt
    cp_tmp = 0
    for p_1 in range(0, K - 1):
        n_p_1 = len(hits[layers[p_1]]) + 1
        for i in range(1, n_p_1):
            cp_tmp += c[p_1, i]

    model.add_constraint(cp >= min_cost)
    model.add_constraint(cp >= cp_tmp)
    total_hits = sum(len(hits[layers[p]]) for p in range(0, K + 1))
    print("total_hits:", total_hits)

    objective = alpha * (total_hits - nt) + gamma * cp

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
        # print("var:", var_name, "-- value:", var_value)
        if 'c' in phi_p_p_i_j[0] or var_value != 1.0 or 's' in phi_p_p_i_j[0] or 'q' in phi_p_p_i_j[0] or 'ob' in \
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
        if h1 in hits[0] or h2 in hits[0]:
            ax.plot(xs=[h1.x, h2.x], ys=[h1.y, h2.y], zs=[h1.z, h2.z], color='green')
        elif h1 in hits[16] or h2 in hits[16]:
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

    hits_volume = read_hits(data_path)
    hits = hits_volume[9]
    source, sink = create_source_sink(hits)
    hits[0] = [source]
    hits[16] = [sink]
    layers = sorted(list(hits.keys()))

    model_path_out = "results/6hits/unknown_track/model_docplex_LB_dist_ss_test.lp"
    solution_path_out = "results/6hits/unknown_track/solution_LB_dist_ss_test.json"
    figure_path_out = "results/6hits/unknown_track/result_LB_dist_ss_test.PNG"

    M = 10000
    alpha = 100000
    gamma = 1
    NT = 6
    run(hits, M, alpha, gamma, model_path_out, solution_path_out, figure_path_out, NT)
