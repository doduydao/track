from data import *
from read_data import *
# import pulp
import cplex
from docplex.mp.model import Model
import json
import random
import matplotlib.pyplot as plt
from docplex.mp.model_reader import ModelReader


def create_variables(model, hits):
    layers = sorted(list(hits.keys()))
    K = len(layers) - 2

    v = []
    for p_1 in range(0, K + 1):
        for p_2 in range(1, K + 2):
            n_p_1 = len(hits[layers[p_1]]) + 1
            n_p_2 = len(hits[layers[p_2]]) + 1
            for i in range(1, n_p_1):
                for j in range(1, n_p_2):
                    v.append((p_1, p_2, i, j))

    phi = model.binary_var_dict(v, name="phi")

    v = []
    for p_1 in range(0, K):
        n_p_1 = len(hits[layers[p_1]]) + 1
        for i in range(1, n_p_1):
            v.append((p_1, i))
    c = model.continuous_var_dict(v, name="c", lb=0)

    ob = model.continuous_var(name="ob")

    return ob, phi, c


def distance(h1, h2):
    distance = math.sqrt((h2.x - h1.x) ** 2 + (h2.y - h1.y) ** 2 + (h2.z - h1.z) ** 2)
    return distance


def run(hits, M, model_path_out, solution_path_out):
    model = Model(name="Track")
    layers = sorted(list(hits.keys()))
    print("layers:", layers)
    K = len(layers) - 2
    print("K=", K)
    # create_variables
    ob, phi, c = create_variables(model, hits)

    nts = [len(hits[l]) for l in layers[1:-1]]
    min_nt = min(nts)
    max_nt = max(nts)

    # first constraints:
    print("---First constraints---")
    tmp = 0
    for p_2 in range(1, K + 1):
        n_p_2 = len(hits[layers[p_2]]) + 1
        for j in range(1, n_p_2):
            tmp += phi[0, p_2, 1, j]
    count_constraint = 1
    c1 = "FC_" + str(count_constraint)
    # model.add_constraint(tmp >= min_nt, ctname=c1 + "1")
    model.add_constraint(tmp == max_nt, ctname=c1 + "2")

    # Second constraints:
    print("---Second constraints---")
    tmp = 0
    for p_1 in range(1, K + 1):
        n_p_1 = len(hits[layers[p_1]]) + 1
        for i in range(1, n_p_1):
            tmp += phi[p_1, K + 1, i, 1]
    count_constraint = 1
    c2 = "SC_" + str(count_constraint)
    # model.add_constraint(tmp >= min_nt, ctname=c2 + "1")
    model.add_constraint(tmp == max_nt, ctname=c2 + "2")

    # Third constraints:
    print("---Third constraints---")
    count_constraint = 0
    for p_1 in range(1, K + 1):
        n_p_1 = len(hits[layers[p_1]]) + 1
        for i in range(1, n_p_1):
            tmp = 0
            for p_2 in range(0, p_1):
                n_p_2 = len(hits[layers[p_2]]) + 1
                for j in range(1, n_p_2):
                    tmp += phi[p_2, p_1, j, i]
            count_constraint += 1
            c3 = "TC_" + str(count_constraint)
            model.add_constraint(tmp <= 1, ctname=c3)

    print("---Fourth, Fifth and Sixth constraints---")
    count_constraint_4 = 0
    count_constraint_5 = 0
    count_constraint_6 = 0
    for p_1 in range(1, K + 1):
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
            for p_2 in range(p_1 + 1, K + 2):
                n_p_2 = len(hits[layers[p_2]]) + 1
                for j in range(1, n_p_2):
                    tmp_2 += phi[p_1, p_2, i, j]
                    t_2 += phi[p_1, p_2, i, j]
            count_constraint_4 += 1
            c4 = "FoC_" + str(count_constraint_4)
            model.add_constraint(t_1 == t_2, ctname=c4)

        count_constraint_5 += 1
        c5 = "FiC_" + str(count_constraint_5)
        model.add_constraint(tmp_1 >= min_nt, ctname=c5)
        count_constraint_5 += 1
        c5 = "FiC_" + str(count_constraint_5)
        model.add_constraint(tmp_1 <= n_p_1 - 1, ctname=c5)
        count_constraint_6 += 1
        c6 = "SiC_" + str(count_constraint_6)
        model.add_constraint(tmp_2 >= min_nt, ctname=c6)
        count_constraint_6 += 1
        c6 = "SiC_" + str(count_constraint_6)
        model.add_constraint(tmp_2 <= n_p_1 - 1, ctname=c6)

    print("--- Seventh and Eighth constraints---")
    count_constraint_7 = 0
    count_constraint_8 = 0
    min_cost = 0
    for p_1 in range(0, K):
        n_p_1 = len(hits[layers[p_1]]) + 1
        for i in range(1, n_p_1):
            min_beta = 1000000000
            max_beta = -1
            for p_2 in range(p_1 + 1, K + 1):
                n_p_2 = len(hits[layers[p_2]]) + 1
                for j in range(1, n_p_2):
                    for p_3 in range(p_2 + 1, K + 2):
                        n_p_3 = len(hits[layers[p_3]]) + 1
                        for k in range(1, n_p_3):
                            h_i = hits[layers[p_1]][i - 1]
                            h_j = hits[layers[p_2]][j - 1]
                            h_k = hits[layers[p_3]][k - 1]
                            seg_1 = Segment(h_j, h_i)
                            seg_2 = Segment(h_j, h_k)
                            beta = Angle(seg_1=seg_1, seg_2=seg_2).angle * (distance(h_i, h_j) + distance(h_j, h_k))
                            if beta < min_beta:
                                min_beta = beta
                            if beta > max_beta:
                                max_beta = beta
                            count_constraint_7 += 1
                            c7 = "SeC_" + str(count_constraint_7)
                            model.add_constraint(
                                beta <= c[p_1, i] + (2 - phi[p_1, p_2, i, j] - phi[p_2, p_3, j, k]) * M, ctname=c7)
            count_constraint_8 += 1
            c8 = "EiC_" + str(count_constraint_8)
            model.add_constraint(c[p_1, i] >= min_beta, ctname=c8)
            count_constraint_8 += 1
            c8 = "EiC_" + str(count_constraint_8)
            model.add_constraint(c[p_1, i] <= max_beta, ctname=c8)
            min_cost += min_beta

    objective = 0
    for p_1 in range(0, K):
        n_p_1 = len(hits[layers[p_1]]) + 1
        for i in range(1, n_p_1):
            objective += c[p_1, i]

    model.add_constraint(objective >= min_cost)
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

    return result['CPLEXSolution']['variables']


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
    hits_path = '../event000001000/volume_id_9/hits-vol_9_FGC_track_min_6_track_with_noise.csv'
    # hits_path = '../event000001000/sel/event000001000-hits-sel-01.csv'
    hits_volume = read_hits(hits_path)
    hits = dict()
    for k, v in hits_volume.items():
        print("Volume id:", k)
        print("No_layers:", len(v))
        hits = v
    source, sink = create_source_sink(hits)
    hits[0] = [source]
    hits[16] = [sink]
    layers = sorted(list(hits.keys()))

    model_path_out = "result_dao_FGC/model_docplex_hits-vol_9_FGC_track_min_6_track_with_noise_2.lp"
    solution_path_out = "result_dao_FGC/solution_hits-vol_9_FGC_track_min_6_track_with_noise_2.json"

    M = 10000
    result = run(hits, M, model_path_out, solution_path_out)
    # with open(solution_path_out, 'r', encoding='utf-8') as f:
    #     result = json.load(f)['CPLEXSolution']['variables']

    segments = []
    for var in result:

        var_name = var['name']
        var_value = round(float(var['value']))
        phi_p_p_i_j = var_name.split('_')
        print("var:", var_name, "-- value:", var_value)
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
    out = "result_dao_FGC/result_hits-vol_9_FGC_track_min_6_track_with_noise_2.PNG"
    display(hits, segments, out)
