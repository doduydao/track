from data import *
from docplex.mp.model import Model
import json
import matplotlib.pyplot as plt
import os


def define_variables(model, costs):
    var = set()
    for cost in costs.keys():
        i_j = cost[0]
        j_k = cost[1]
        var.add(i_j)
        var.add(j_k)
    var = sorted(var, key=lambda x: (x[0], x[1]))
    x = model.binary_var_dict(var, name='x')

    ob = model.continuous_var(name="ob")

    ss = model.continuous_var(name="ss")
    fp = model.continuous_var(name="fp")
    sp = model.continuous_var(name="sp")
    return x, ob, ss, fp, sp


def run(list_hits, costs, m, M, model_path_out, solution_path_out, figure_path_out):
    # define model
    model = Model(name="Track")

    # create variables
    x, ob, ss, fp, sp = define_variables(model, costs)

    # create objective function
    hit_last_layer = hits_by_layers[m]
    N = len(list_hits) - len(hit_last_layer)
    print(N)
    first_part = 0
    sum_segments = 0
    for id, cost in costs.items():
        i_j = id[0]
        j_k = id[1]
        first_part += cost * x[i_j] * x[j_k]
        sum_segments += x[i_j]

    second_part = M * (((sum_segments / N) - 1) ** 2)

    # third_part
    for k in list(x.keys()):
        i = k[0]
        t_1 = 0
        for k_1 in list(x.keys()):
            if i == k_1[0]:
                t_1 += x[k_1]
        model.add_constraint(t_1 == 1)

    # fourth_part
    for k in list(x.keys()):
        j = k[1]
        t_2 = 0
        for k_1 in list(x.keys()):
            if j == k_1[1]:
                t_2 += x[k_1]
        model.add_constraint(t_2 == 1)

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
    model.add_constraint(ss >= sum_segments)
    model.add_constraint(fp <= first_part)
    model.add_constraint(sp <= second_part)

    model.add_constraint(ob <= fp + sp)
    model.set_objective("max", ob)

    model.print_information()
    model.solve(log_output=True)

    model.export_as_lp(model_path_out)
    model.solution.export(solution_path_out)
    f = open(solution_path_out)
    result = json.load(f)
    f.close()

    result = result['CPLEXSolution']['variables']

    display(list_hits, result, figure_path_out)


def display(hits, result, out=""):
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

    for var in result:
        print(var)
        x_i_j = var['name'].split('_')
        if 'x' in x_i_j[0]:
            i = int(x_i_j[1])
            j = int(x_i_j[2])
            h1 = list_hits[i]
            h2 = list_hits[j]
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
    # all_segments = []
    for l in layers[1:-1]:
        hits_l = hits[l]
        first_part = []
        for i in range(max(layers[0], l - 3), l):
            hits_i = hits[i]
            segs_l_i = build_segments(hits_i, hits_l, list_hits)
            first_part.append(segs_l_i)
            # all_segments += segs_l_i

        second_part = []
        for i in range(l + 1, min(l + 3, layers[-1]) + 1):
            hits_i = hits[i]
            segs_l_i = build_segments(hits_l, hits_i, list_hits)
            second_part.append(segs_l_i)
            # all_segments += segs_l_i

        for f in first_part:
            for s in second_part:
                for seg_f in f:
                    for seg_s in s:
                        if seg_f.hit_2.hit_id == seg_s.hit_1.hit_id:
                            cost = Cost(seg_f, seg_s)
                            cos_beta = cost.cos_beta

                            # if cost.id == ((62, 94), (94, 124)):
                            #     print(cost.id,"cos_beta:",cos_beta, -(cos_beta**7)/cost.sum_distance)
                            #
                            # if cost.id == ((67, 94), (94, 121)):
                            #     print(cost.id,"cos_beta:", cos_beta, -(cos_beta**7)/cost.sum_distance)

                            if cos_beta >= math.cos(math.pi / 500):
                                costs.append(cost)
    all_segments = set()
    for cost in costs:
        all_segments.add(cost.seg_1.id)
        all_segments.add(cost.seg_2.id)

    print("number of segments:", len(all_segments))
    return costs


def write_costs(costs, path, m):
    data = dict()
    for cost in costs:
        cos_beta = cost.cos_beta
        sum_distance = cost.sum_distance
        str_key = str(cost.id)
        data[str_key] = (cos_beta ** m) / sum_distance
    with open(path, 'w') as f:
        json.dump(data, f, ensure_ascii=False, indent=4)
    print('wrote!')


def load_costs(path):
    with open(path) as f:
        data = json.load(f)
    costs = dict()
    for k, v in data.items():
        costs[eval(k)] = v

    return costs


def check_path(path):
    if os.path.exists(path) == False:
        os.mkdir(path)
    else:
        print("Folder exist")


if __name__ == '__main__':
    data_selected_path = '../../src/data_selected'
    out_path = '/Users/doduydao/daodd/PycharmProjects/track/src/Bogdan_idea/results'
    folder = '/100hits/known_track/'
    check_path(out_path + folder)
    data_path = data_selected_path + folder + 'hits.csv'
    costs_path_out = out_path + folder + "costs.json"
    model_path_out = out_path + folder + "model.lp"
    solution_path_out = out_path + folder + "solution.json"
    figure_path_out = out_path + folder + "result.PNG"

    hits_by_layers = read_hits(data_path)[9]
    list_hits = []
    for hs in list(hits_by_layers.values()):
        list_hits += hs
    # costs = get_costs(list_hits, hits_by_layers)
    m = 7
    M = 1
    # write_costs(costs, costs_path_out, m)
    costs = load_costs(costs_path_out)
    result = run(list_hits, costs, m, M, model_path_out, solution_path_out, figure_path_out)
