from data import *
from read_data import *
from docplex.mp.model import Model
import json
import matplotlib.pyplot as plt


def distance(h1, h2):
    distance = math.sqrt((h2.x - h1.x) ** 2 + (h2.y - h1.y) ** 2 + (h2.z - h1.z) ** 2)
    return distance


def define_variables(model, hits):
    layers = list(hits.keys())
    L = len(layers) + 1
    layers = [0] + layers + [0]
    v = []
    for p in range(1, L - 1):
        n_p = len(hits[layers[p]]) + 1
        n_p_1 = len(hits[layers[p + 1]]) + 1
        for i in range(1, n_p):
            for j in range(1, n_p_1):
                v.append((p, i, j))

    f = model.binary_var_dict(v, name="f")

    v = []
    for p_1 in range(1, L - 2):
        n_p_1 = len(hits[layers[p_1]]) + 1
        for i in range(1, n_p_1):
            v.append((p_1, i))

    c = model.continuous_var_dict(v, name="c", lb=0)
    nt = model.continuous_var(name="nt")
    cp = model.continuous_var(name="cp")
    ob = model.continuous_var(name="ob")
    return f, c, nt, cp, ob


def run(hits, NT, M, alpha, gamma, model_path_out, solution_path_out, figure_path_out):
    model = Model(name="Track")
    layers = list(hits.keys())
    L = len(layers) + 1
    layers = [0] + layers + [0]

    f, c, nt, cp, ob = define_variables(model, hits)
    nt = NT
    # Constraints
    print("---First constraints---")
    cc_1 = 0
    for p in range(2, L):
        n_p = len(hits[layers[p]]) + 1
        n_p_1 = len(hits[layers[p - 1]]) + 1
        for i in range(1, n_p):
            tmp = 0
            for j in range(1, n_p_1):
                tmp += f[p - 1, j, i]
            cc_1 += 1
            cn = "FC_" + str(cc_1)
            model.add_constraint(tmp <= 1, ctname=cn)

    print("---Second constraints---")
    cc_2 = 0
    for p in range(1, L - 1):
        n_p = len(hits[layers[p]]) + 1
        n_p_1 = len(hits[layers[p + 1]]) + 1
        for i in range(1, n_p):
            tmp = 0
            for j in range(1, n_p_1):
                tmp += f[p, i, j]
            cc_2 += 1
            cn = "SC_" + str(cc_2)
            model.add_constraint(tmp <= 1, ctname=cn)

    print("---Third constraints---")
    cc_3 = 0

    n_1 = len(hits[layers[1]]) + 1
    n_2 = len(hits[layers[2]]) + 1
    tmp = 0
    for i in range(1, n_1):
        for j in range(1, n_2):
            tmp += f[1, i, j]
    cc_3 += 1
    cn = "TC_" + str(cc_3)
    model.add_constraint(tmp == nt, ctname=cn)

    print("---Fourth constraints---")
    cc_4 = 0
    for p in range(2, L - 1):
        n_p = len(hits[layers[p]]) + 1
        n_p_1 = len(hits[layers[p - 1]]) + 1
        n_p_2 = len(hits[layers[p + 1]]) + 1
        for i in range(1, n_p):
            tmp_1 = 0
            tmp_2 = 0
            for j in range(1, n_p_1):
                tmp_1 += f[p - 1, j, i]

            for j in range(1, n_p_2):
                tmp_2 += f[p, i, j]

            cc_4 += 1
            cn = "FoC_" + str(cc_4)
            model.add_constraint(tmp_1 == tmp_2, ctname=cn)

    print("---Fifth constraints---")
    cc_5 = 0
    cost_LB = 0
    betas = []
    for p in range(1, L - 2):
        n_p = len(hits[layers[p]]) + 1
        n_p_1 = len(hits[layers[p + 1]]) + 1
        n_p_2 = len(hits[layers[p + 2]]) + 1
        min_beta = 1000000000
        for i in range(1, n_p):
            for j in range(1, n_p_1):
                for k in range(1, n_p_2):
                    h_i = hits[layers[p]][i - 1]
                    h_j = hits[layers[p + 1]][j - 1]
                    h_k = hits[layers[p + 2]][k - 1]
                    seg_1 = Segment(h_j, h_i)
                    seg_2 = Segment(h_j, h_k)
                    angle = Angle(seg_1=seg_1, seg_2=seg_2).angle
                    dist = distance(h_i, h_j) // 100 + distance(h_j, h_k) // 100
                    beta = angle * dist
                    betas.append(beta)
                    if min_beta > beta:
                        min_beta = beta
                    cc_5 += 1
                    cn = "FiC_" + str(cc_5)
                    model.add_constraint(
                        beta <= c[p, i] + ((2 - f[p, i, j] - f[p + 1, j, k]) * M),
                        ctname=cn)
        cost_LB += min_beta

    beta_upper = sum(betas) / len(betas)
    print("beta upper =", beta_upper)
    cp_tmp = 0
    for p in range(1, L - 2):
        n_p = len(hits[layers[p]]) + 1
        for i in range(1, n_p):
            model.add_constraint(c[p, i] <= beta_upper)
            cp_tmp += c[p, i]

    objective = alpha * (len(hits[layers[1]]) - nt) + gamma * cp
    model.add_constraint(cp >= cost_LB * nt)
    model.add_constraint(cp >= cp_tmp)

    model.add_constraint(ob >= objective)
    model.set_objective('min', ob)

    model.print_information()
    model.solve(log_output=True)

    model.export_as_lp(model_path_out)
    model.solution.export(solution_path_out)
    f = open(solution_path_out)
    result = json.load(f)
    f.close()

    result = result['CPLEXSolution']['variables']

    segments = []

    for var in result:
        f_p_i_j = var['name'].split('_')
        value = int(round(float(var['value'])))
        if f_p_i_j[0] == 'z' or value != 1:
            continue
        print(var)
        p = int(f_p_i_j[1])
        i = int(f_p_i_j[2])
        j = int(f_p_i_j[3])

        h_1 = hits[layers[p]][i - 1]
        h_2 = hits[layers[p + 1]][j - 1]
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

    model_path_out = "results/15hits/unknown_track/model_docplex_strong_LB.lp"
    solution_path_out = "results/15hits/unknown_track/solution_strong_LB.json"
    figure_path_out = "results/15hits/unknown_track/result_strong_LB.PNG"
    M = 10000
    alpha = 1000
    gamma = 1
    NT = 15
    result = run(hits, NT, M, alpha, gamma, model_path_out, solution_path_out, figure_path_out)
