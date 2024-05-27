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
    return ob, phi, c


def distance(h1, h2):
    distance = math.sqrt((h2.x - h1.x) ** 2 + (h2.y - h1.y) ** 2 + (h2.z - h1.z) ** 2)
    return distance


def run(hits, M, model_path_out, solution_path_out, figure_path_out):
    model = Model(name="Track")
    layers = sorted(list(hits.keys()))
    K = len(layers) + 1
    layers = [0] + layers + [0]

    # create_variables
    ob, phi, c = create_variables(model, hits)

    nts = [len(hits[l]) for l in layers[1:-1]]
    n_min = min(nts)

    # First constraints:
    print("---First and second constraints---")
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

    cc_2 = 1
    for p_1 in range(2, K - 1):
        n_p_1 = len(hits[layers[p_1]]) + 1
        for i in range(1, n_p_1):
            tmp = 0
            for p_2 in range(p_1 + 1, K):
                n_p_2 = len(hits[layers[p_2]]) + 1
                for j in range(1, n_p_2):
                    tmp += phi[p_1, p_2, i, j]
            cc_2 += 1
            cn = "SC_" + str(cc_2)
            model.add_constraint(tmp <= 1, ctname=cn)

    # First constraints:
    print("---Third and fourth constraints---")
    cc_3 = 0
    n_p_2 = len(hits[layers[K - 1]]) + 1
    for j in range(1, n_p_2):
        tmp = 0
        for p_1 in range(2, K - 1):
            n_p_1 = len(hits[layers[p_1]]) + 1
            for i in range(1, n_p_1):
                tmp += phi[p_1, K - 1, i, j]
        cc_3 += 1
        cn = "TC_" + str(cc_3)
        model.add_constraint(tmp <= 1, ctname=cn)

    cc_4 = 0
    for p_2 in range(2, K - 1):
        n_p_2 = len(hits[layers[p_2]]) + 1
        for j in range(1, n_p_2):
            tmp = 0
            for p_1 in range(1, p_2):
                n_p_1 = len(hits[layers[p_1]]) + 1
                for i in range(1, n_p_1):
                    tmp += phi[p_1, p_2, i, j]
            cc_4 += 1
            cn = "FoC_" + str(cc_4)
            model.add_constraint(tmp <= 1, ctname=cn)

    print("---Fifth constraints---")
    cc_5 = 0
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
            cc_5 += 1
            cn = "FiC_" + str(cc_5)
            model.add_constraint(t_1 == t_2, ctname=cn)

    print("---Sixth constraints---")
    cc_6 = 0
    for p_1 in range(2, K):
        n_p_1 = len(hits[layers[p_1]]) + 1
        tmp = 0
        for i in range(1, n_p_1):
            for p_2 in range(1, p_1):
                n_p_2 = len(hits[layers[p_2]]) + 1
                for j in range(1, n_p_2):
                    tmp += phi[p_2, p_1, j, i]
        cc_6 += 1
        cn = "SiC_" + str(cc_6)
        model.add_constraint(n_min <= tmp, ctname=cn)
        cc_6 += 1
        cn = "SiC_" + str(cc_6)
        model.add_constraint(tmp <= n_p_1, ctname=cn)

    print("---Seventth constraints---")
    cc_7 = 0
    for p_1 in range(1, K - 1):
        n_p_1 = len(hits[layers[p_1]]) + 1
        tmp = 0
        for i in range(1, n_p_1):
            for p_2 in range(p_1 + 1, K):
                n_p_2 = len(hits[layers[p_2]]) + 1
                for j in range(1, n_p_2):
                    tmp += phi[p_1, p_2, i, j]
        cc_7 += 1
        cn = "SeC_" + str(cc_7)
        model.add_constraint(n_min <= tmp, ctname=cn)
        cc_7 += 1
        cn = "SeC_" + str(cc_7)
        model.add_constraint(tmp <= n_p_1, ctname=cn)

    print("--- Eighth constraints---")
    cc_8 = 0

    min_cost = 0
    for p_1 in range(1, K - 2):
        n_p_1 = len(hits[layers[p_1]]) + 1
        min_beta = 1000000000
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
                            beta = Angle(seg_1=seg_1, seg_2=seg_2).angle

                            if beta < min_beta:
                                min_beta = beta

                            cc_8 += 1
                            cn = "EiC_" + str(cc_8)
                            model.add_constraint(
                                beta <= c[p_1, i] + (2 - phi[p_1, p_2, i, j] - phi[p_2, p_3, j, k]) * M, ctname=cn)
        min_cost += min_beta

    objective = 0
    for p_1 in range(1, K - 2):
        n_p_1 = len(hits[layers[p_1]]) + 1
        for i in range(1, n_p_1):
            objective += c[p_1, i]

    min_cost = n_min * min_cost
    model.add_constraint(objective >= min_cost, ctname="LB of objective value")
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
        ax.plot(xs=[h1.x, h2.x], ys=[h1.y, h2.y], zs=[h1.z, h2.z], color='blue')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    plt.savefig(out)
    plt.show()


if __name__ == '__main__':
    src_path = '../data_selected'
    data_path = src_path + '/15hits/unknown_track/hits.csv'

    hits_volume = read_hits(data_path)
    hits = hits_volume[9]

    model_path_out = "results/6hits/unknown_track/model_docplex_LB_new.lp"
    solution_path_out = "results/15hits/unknown_track/solution_new.json"
    figure_path_out = "results/15hits/unknown_track/result_new.PNG"

    M = 10000
    result = run(hits, M, model_path_out, solution_path_out, figure_path_out)
