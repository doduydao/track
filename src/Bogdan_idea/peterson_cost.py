from data import *
from docplex.mp.model import Model
import json
import matplotlib.pyplot as plt


def define_variables(model, costs):
    var = set()
    for cost in costs:
        i_j = cost.id[0]
        j_k = cost.id[1]
        var.add(i_j)
        var.add(j_k)
    var = sorted(var, key=lambda x: (x[0], x[1]))
    x = model.binary_var_dict(var, name='x')

    # ob = model.continuous_var(name="ob")
    return x


def run(list_hits, costs, m, M, model_path_out, solution_path_out, figure_path_out):
    # define model
    model = Model(name="Track")

    # create variables
    x = define_variables(model, costs)

    # create objective function
    N = len(list_hits)

    first_part = 0
    sum_segments = 0
    for cost in costs:
        i_j = cost.id[0]
        j_k = cost.id[1]
        cos_beta = cost.cos_beta
        sum_distance = cost.sum_distance
        first_part += (-(cos_beta ** m) / sum_distance) * x[i_j] * x[j_k]
        sum_segments += x[i_j]

    second_part = M * ((sum_segments - N) ** 2)

    i_s = set()
    j_s = set()
    for k in list(x.keys()):
        i_s.add(k[0])
        j_s.add(k[1])

    # third_part = 0
    for i in i_s:
        t_1 = 0
        for k in list(x.keys()):
            if i == k[0]:
                t_1 += x[k]
        model.add_constraint(t_1 == 1)

    # fourth_part = 0
    for j in j_s:
        t_2 = 0
        for k in list(x.keys()):
            if j == k[1]:
                t_2 += x[k]
        model.add_constraint(t_2 == 1)

    fifth_part = 0
    for k in list(x.keys()):
        j = k[1]
        t_1 = 0
        for k_1 in list(x.keys()):
            if k_1[1] == j:
                t_1 += x[k_1]
        t_2 = 0
        for k_2 in list(x.keys()):
            if k_2[0] == j:
                t_2 += x[k_2]
        if str(t_1) != '0' and str(t_2) != '0':
            model.add_constraint(t_1 == t_2)

    ob = first_part + second_part

    model.set_objective("min", ob)

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
        print(var)

        x_i_j = var['name'].split('_')
        if 'x' in x_i_j[0]:
            i = int(x_i_j[1])
            j = int(x_i_j[2])
            h_1 = list_hits[i]
            h_2 = list_hits[j]
            segments.append([h_1, h_2])

    display(list_hits, segments, figure_path_out)


def display(hits, segments, out=""):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    xs = []
    ys = []
    zs = []

    for h in hits:
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


def build_segments(hits_1, hits_2, list_hits):
    segments = []
    for i in hits_1:
        for j in hits_2:
            if i.index is None:
                index_i = list_hits.index(i)
                i.set_index(index_i)
            if j.index is None:
                index_j = list_hits.index(j)
                j.set_index(index_j)
            segments.append(Segment(i, j))
    return segments


def get_costs(list_hits, hits):
    layers = list(hits.keys())
    costs = []
    all_segments = []
    for l in layers[1:-1]:
        hits_l = hits[l]
        first_part = []
        for i in range(max(layers[0], l - 3), l):
            hits_i = hits[i]
            segs_l_i = build_segments(hits_i, hits_l, list_hits)
            first_part.append(segs_l_i)
            all_segments += segs_l_i

        second_part = []
        for i in range(l + 1, min(l + 3, layers[-1]) + 1):
            hits_i = hits[i]
            segs_l_i = build_segments(hits_l, hits_i, list_hits)
            second_part.append(segs_l_i)
            all_segments += segs_l_i

        for f in first_part:
            for s in second_part:
                for seg_f in f:
                    for seg_s in s:
                        if seg_f.hit_2.hit_id == seg_s.hit_1.hit_id:
                            costs.append(Cost(seg_f, seg_s))
    print(len(all_segments))
    return costs


if __name__ == '__main__':
    src_path = '../../src/data_selected'
    folder = '/1hits/'
    data_path = src_path + folder + 'known_track/hits.csv'
    hits_by_layers = read_hits(data_path)[9]
    list_hits = []
    for hs in list(hits_by_layers.values()):
        list_hits += hs

    costs = get_costs(list_hits, hits_by_layers)

    # for cost in costs:
    #     cos_beta = cost.cos_beta
    #     sum_distance = cost.sum_distance
    #     print("cost = ", (-cos_beta ** 7) / sum_distance)

    model_path_out = "results/1hits/known_track/model.lp"
    solution_path_out = "results/1hits/known_track/solution.json"
    figure_path_out = "results/1hits/known_track/result.PNG"

    m = 1
    M = 1
    result = run(list_hits, costs, m, M, model_path_out, solution_path_out, figure_path_out)