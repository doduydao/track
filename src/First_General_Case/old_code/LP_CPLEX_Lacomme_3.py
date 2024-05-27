from data import *
from read_data import *
# import pulp
import cplex
from docplex.mp.model import Model
import json
import random
import matplotlib.pyplot as plt


def create_variables(model, hits):
    layers = sorted(list(hits.keys()))
    K = len(layers) - 2

    # create phi_p_p'_i_j variables
    v = []
    for p_1 in range(0, K + 1):
        for p_2 in range(1, K + 2):
            n_p_1 = len(hits[layers[p_1]]) + 1
            n_p_2 = len(hits[layers[p_2]]) + 1
            for i in range(1, n_p_1):
                for j in range(1, n_p_2):
                    v.append((p_1, p_2, i, j))

    phi = model.binary_var_dict(v, name="phi")
    print("No_phi_variables: ", len(phi))
    # print(phi)
    v = []
    for p_2 in range(1, K + 1):
        n_p_2 = len(hits[layers[p_2]]) + 1
        for j in range(1, n_p_2):
            v.append((p_2, j))
    c = model.continuous_var_dict(v, name="c", lb=0)
    print("No_phi_variables: ", len(c))
    c[0, 1] = 0
    v = []
    for p_2 in range(0, K + 1):
        n_p_2 = len(hits[layers[p_2]]) + 1
        for j in range(1, n_p_2):
            v.append((p_2, j))
    s = model.continuous_var_dict(v, name="s", lb=0)
    s[0, 1] = 0
    print("No_phi_variables: ", len(c))

    v = []
    for p in range(1, K + 1):
        v.append(p)
    q = model.continuous_var_dict(v, name="q", lb=0)

    ob = model.continuous_var(name="ob")
    # print("No_phi_variables: ", len(c))

    return phi, c, s, ob, q

def distance(h1, h2):
    distance = math.sqrt((h2.x - h1.x) ** 2 + (h2.y - h1.y) ** 2 + (h2.z - h1.z) ** 2)
    return distance


def run(hits, nt, M, model_path_out, solution_path_out):
    model = Model(name="Track")
    layers = sorted(list(hits.keys()))
    print("layers:", layers)
    K = len(layers) - 2
    print("K=", K)
    # create_variables
    phi, c, s, ob, q = create_variables(model, hits)

    # add constraints
    # first constraints:
    print("---First constraints---")
    count_constraint = 0
    tmp = 0
    for p_2 in range(1, K + 1):
        n_p_2 = len(hits[layers[p_2]]) + 1
        for j in range(1, n_p_2):
            tmp += phi[0, p_2, 1, j]
    c1 = "FC_" + str(count_constraint)
    count_constraint += 1
    model.add_constraint(tmp == nt, ctname=c1)
    print("Number of first constraints:", count_constraint)

    # Second constraints:
    print("---Second constraints---")
    count_constraint = 0
    tmp = 0
    for p_1 in range(1, K + 1):
        n_p_1 = len(hits[layers[p_1]]) + 1
        for i in range(1, n_p_1):
            tmp += phi[p_1, K + 1, i, 1]
    c2 = "SC_" + str(count_constraint)
    count_constraint += 1
    model.add_constraint(tmp == nt, ctname=c2)
    print("Number of second constraints:", count_constraint)

    # Third constraints:
    print("---Third constraints---")
    count_constraint = 0
    for p_1 in range(1, K + 1):
        n_p_1 = len(hits[layers[p_1]]) + 1
        for i in range(1, n_p_1):
            tmp_1 = 0
            for p_2 in range(0, p_1):
                n_p_2 = len(hits[layers[p_2]]) + 1
                for j in range(1, n_p_2):
                    tmp_1 += phi[p_2, p_1, j, i]

            tmp_2 = 0
            for p_2 in range(p_1 + 1, K + 2):
                n_p_2 = len(hits[layers[p_2]]) + 1
                for j in range(1, n_p_2):
                    tmp_2 += phi[p_1, p_2, i, j]
            c3 = "TC_" + str(count_constraint)
            count_constraint += 1
            model.add_constraint(tmp_1 == tmp_2, ctname=c3)
    print("Number of Third constraints:", count_constraint)

    print("---Fourth constraints---")
    count_constraint = 1
    for p_1 in range(1, K + 1):
        n_p_1 = len(hits[layers[p_1]]) + 1
        for i in range(1, n_p_1):
            tmp = 0
            for p_2 in range(0, p_1):
                n_p_2 = len(hits[layers[p_2]]) + 1
                for j in range(1, n_p_2):
                    tmp += phi[p_2, p_1, j, i]
            c4 = "FoC_" + str(count_constraint)
            count_constraint += 1
            model.add_constraint(tmp <= 1, ctname=c4)
    print("Number of Fourth constraints:", count_constraint)

    print("--- Fiveth constraints---")
    count_constraint = 1
    for p_1 in range(0, K):
        n_p_1 = len(hits[layers[p_1]]) + 1
        for i in range(1, n_p_1):
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
                            beta = Angle(seg_1=seg_1, seg_2=seg_2).angle * ((h_k.z//100 - h_i.z//100) ** 2 + (50 - h_k.z//100) )
                            print("beta", beta)
                            c7 = "FiC_" + str(count_constraint)
                            count_constraint += 1
                            model.add_constraint(
                                beta <= (2 - phi[p_1, p_2, i, j] - phi[p_2, p_3, j, k]) * M + c[p_2, j], ctname=c7)
    print("Number of Fiveth constraints:", count_constraint)

    for p_1 in range(0, K):
        n_p_1 = len(hits[layers[p_1]]) + 1
        for i in range(1, n_p_1):
            for p_2 in range(p_1 + 1, K + 1):
                n_p_2 = len(hits[layers[p_2]]) + 1
                for j in range(1, n_p_2):
                    # print("1:", s[p_1, i])
                    # print("2:",c[p_1, j])
                    # print("3:",s[p_2, j])
                    model.add_constraint(s[p_1, i] + c[p_2, j] <= s[p_2, j] + (1 - phi[p_1, p_2, i, j]) * M)
    for p in range(1, K + 1):
        n_p = len(hits[layers[p]]) + 1
        tmp = 0
        for i in range(1, n_p):
            tmp += s[p, i]
        model.add_constraint(q[p] == tmp)

    for p in range(2, K + 1):
        model.add_constraint(q[p - 1] <= q[p])

    objective = 0
    print("K=", K)
    for p_1 in range(1, K + 1):
        n_p_1 = len(hits[layers[p_1]]) + 1
        for i in range(1, n_p_1):
            objective += c[p_1, i]

    sum_s = 0
    for p in range(1, K + 1):
        n_p_1 = len(hits[layers[p]]) + 1
        for i in range(1, n_p_1):
            sum_s += s[p, i]
    # model.add_constraint(sum_s >= 700)

    model.add_constraint(ob >= q[K])
    model.add_constraint(ob >= objective)
    model.set_objective('min', ob)
    # model.add_constraint(ob <= 300)

    # model.add_constraint(objective >= 0.009, ctname="LB")

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
        z=sum([h.z for h in first_layer]) / len(first_layer) - 1,
        volume_id=0,
        layer_id=0,
        module_id=0
    )
    sink = Hit(
        hit_id=10,
        x=sum([h.x for h in last_layer]) / len(last_layer),
        y=sum([h.y for h in last_layer]) / len(last_layer),
        z=sum([h.z for h in last_layer]) / len(last_layer) + 1,
        volume_id=0,
        layer_id=0,
        module_id=0
    )

    return source, sink


if __name__ == '__main__':
    hits_path = '../event000001000/sel/event000001000-hits-sel-01.csv'
    hits_volume = read_hits(hits_path)
    hits = dict()
    for k, v in hits_volume.items():
        print("Volume id:", k)
        print("No_layers:", len(v))
        hits = v

    # for p, hp in hits.items():
    #     z = hp[0].z // 100
    #
    #
    #     print("layer: ", p, "--", )
    #     # print()

    source, sink = create_source_sink(hits)
    hits[0] = [source]
    hits[16] = [sink]
    # print(hits)
    layers = sorted(list(hits.keys()))
    # print(layers)
    model_path_out = "result_f2_Lacomme/model_docplex.lp"
    solution_path_out = "result_f2_Lacomme/solution.json"

    nt = 6
    M = 10000
    result = run(hits, nt, M, model_path_out, solution_path_out)

    # del hits[0]
    # del hits[16]
    segments = []
    for var in result:
        # print(var)
        var_name = var['name']
        var_value = var['value']
        phi_p_p_i_j = var_name.split('_')
        print(var)
        if 'c' in phi_p_p_i_j[0] or var['value'] != '1.0':
            continue

        p_1 = int(phi_p_p_i_j[1])
        p_2 = int(phi_p_p_i_j[2])

        # if p_1 == 0 or p_2 == len(layers)-1:
        #     continue

        i = int(phi_p_p_i_j[3])
        j = int(phi_p_p_i_j[4])
        h_1 = hits[layers[p_1]][i - 1]
        h_2 = hits[layers[p_2]][j - 1]
        segments.append([h_1, h_2])
    out = "result_f2_Lacomme/result.PNG"
    display(hits, segments, out)
