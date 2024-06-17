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


def run(list_hits, costs, m, alpha, beta, model_path_out, solution_path_out, figure_path_out):
    # define model
    model = Model(name="Track")

    # create variables
    x, ob, ss, fp, sp = define_variables(model, costs)

    # create objective function
    hit_last_layer = hits_by_layers[m]
    N = len(list_hits) - len(hit_last_layer)
    print("N =", N)
    first_part = 0
    segments = set()
    for id, cost in costs.items():
        i_j = id[0]
        j_k = id[1]
        first_part += cost * x[i_j] * x[j_k]
        segments.add(i_j)
        segments.add(j_k)

    sum_segments = sum([x[s] for s in segments])
    # second_part
    second_part = (sum_segments - N) ** 2

    # third_part
    third_part = 0
    for k in list(x.keys()):
        i = k[0]
        t_1 = 0
        for k_1 in list(x.keys()):
            if i == k_1[0]:
                t_1 += x[k_1]
        third_part += (1 - t_1) ** 2

    # fourth_part
    fourth_part = 0
    for k in list(x.keys()):
        j = k[1]
        t_2 = 0
        for k_1 in list(x.keys()):
            if j == k_1[1]:
                t_2 += x[k_1]
        fourth_part += (1 - t_2) ** 2

    ob = -first_part + alpha * second_part + beta * (third_part + fourth_part)

    # model.add_constraint(ss <= sum_segments)
    # model.add_constraint(fp <= first_part)
    # model.add_constraint(ob <= first_part)
    model.set_objective("min", ob)

    model.print_information()
    model.solve(log_output=True)

    model.export_as_lp(model_path_out)
    if model.solution == None:
        print("No solution!")
    else:
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

    count_x = 0
    for var in result:
        print(var)
        x_i_j = var['name'].split('_')
        if 'x' in x_i_j[0] and round(float(var['value'])) == 1.0:
            count_x += 1
            i = int(x_i_j[1])
            j = int(x_i_j[2])
            h1 = list_hits[i]
            h2 = list_hits[j]
            ax.plot(xs=[h1.x, h2.x], ys=[h1.y, h2.y], zs=[h1.z, h2.z], color='blue')
    print("count x:", count_x)
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


def get_costs(list_hits, hits, beta_max):
    layers = list(hits.keys())
    print("Layers: ", layers)
    costs = []
    single_segs = []
    for l in layers[0:-1]:
        print("Combine h-layer: ", l)
        hits_l = hits[l]
        for i in range(l + 1, min(l + 3, layers[-1]) + 1):
            print("   with h-layer: ", i)
            hits_i = hits[i]
            segs = build_segments(hits_l, hits_i, list_hits)
            single_segs.append(segs)

    for segs_1 in single_segs:
        for segs_2 in single_segs:
            for seg_1 in segs_1:
                for seg_2 in segs_2:
                    if ((seg_1.hit_2.hit_id == seg_2.hit_1.hit_id) & ((seg_1.d_l + seg_2.d_l) <= 4) & (
                            seg_1.layer < seg_2.layer)):
                        # print("Combine segment ", seg_1.layer, seg_1.d_l, " with segment " , seg_2.layer, seg_2.d_l, " , Hits: ", seg_1.hit_1.hit_id, seg_1.hit_2.hit_id, seg_2.hit_1.hit_id, seg_2.hit_2.hit_id);
                        cost = Cost(seg_1, seg_2)
                        cos_beta = cost.cos_beta
                        if cos_beta >= math.cos(beta_max):
                            costs.append(cost)
    all_segments = set()
    for cost in costs:
        all_segments.add(cost.seg_1)
        all_segments.add(cost.seg_2)
    print("Number of segments:", len(all_segments))
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
    folder = '/6hits/known_track/'
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

    beta_max = math.pi / 100
    m = 7
    alpha = 1
    beta = 1
    costs = get_costs(list_hits, hits_by_layers, beta_max)
    write_costs(costs, costs_path_out, m)
    costs = load_costs(costs_path_out)
    result = run(list_hits, costs, m, alpha, beta, model_path_out, solution_path_out, figure_path_out)
