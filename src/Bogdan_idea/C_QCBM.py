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
    return x


def run(list_hits, costs, m, model_path_out, solution_path_out, figure_path_out):
    # define model
    model = Model(name="Track")

    # create variables
    x = define_variables(model, costs)

    # create objective function
    hit_last_layer = hits_by_layers[m]
    N = len(list_hits) - len(hit_last_layer)
    print("N=", N)
    ob = 0
    segments = set()
    for id, cost in costs.items():
        i_j = id[0]
        j_k = id[1]
        ob += cost * x[i_j] * x[j_k]
        segments.add(i_j)
        segments.add(j_k)

    # second_part
    sum_segments = sum([x[s] for s in segments])
    ctn = "SP" + str(1)
    model.add_constraint(sum_segments == N, ctname=ctn)

    t_1 = dict()
    t_2 = dict()
    for k in x.keys():
        i = k[0]
        j = k[1]
        if i not in t_1:
            t_1[i] = {j}
        else:
            t_1[i].add(j)

        if j not in t_2:
            t_2[j] = {i}
        else:
            t_2[j].add(i)

    # third_part
    ct = 0
    for i, v in t_1.items():
        tmp = 0
        for j in v:
            tmp += x[(i, j)]
        ct += 1
        ctn = "TP" + str(ct)
        model.add_constraint(tmp == 1, ctname=ctn)

    # fourth_part
    ct = 0
    for j, v in t_2.items():
        tmp = 0
        for i in v:
            tmp += x[(i, j)]
        ct += 1
        ctn = "FP" + str(ct)
        model.add_constraint(tmp == 1, ctname=ctn)

    model.set_objective("min", -ob)
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


def cal_expected_value(list_hits):
    track = dict()
    for hit in list_hits:
        k = hit.particle_id / 1000
        if k not in track:
            track[k] = [hit]
        else:
            track[k].append(hit)
    cost = 0
    for t, hs in track.items():
        for i in range(len(hs) - 2):
            h_i = hs[i]
            h_j = hs[i + 1]
            h_k = hs[i + 2]

            seg_1 = Segment(h_j, h_i)
            seg_2 = Segment(h_j, h_k)
            c = Cost(seg_1, seg_2)
            # print(c.cos_beta ** 7 / c.sum_distance)
            cost += c.cos_beta ** 7 / c.sum_distance

    print("expected cost=", cost)


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
    costs = []
    for l in layers[1:-1]:
        hits_j = hits[l]
        first_part = []
        for i in range(max(layers[0], l - 3), l):
            hits_i = hits[i]
            segs_j_i = build_segments(hits_i, hits_j, list_hits)
            first_part.append(segs_j_i)

        second_part = []
        for k in range(l + 1, min(l + 3, layers[-1]) + 1):
            hits_k = hits[k]
            segs_j_k = build_segments(hits_j, hits_k, list_hits)
            second_part.append(segs_j_k)

        for f in first_part:
            for s in second_part:
                for seg_f in f:
                    for seg_s in s:
                        if seg_f.hit_2.hit_id == seg_s.hit_1.hit_id:
                            if seg_f.gaps + seg_s.gaps <= 4:
                                cost = Cost(seg_f, seg_s)
                                cos_beta = cost.cos_beta
                                if cos_beta >= math.cos(beta_max):
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
    model_path_out = out_path + folder + "model_C_QCBM.lp"
    solution_path_out = out_path + folder + "solution_C_QCBM.json"
    figure_path_out = out_path + folder + "result_C_QCBM.PNG"

    hits_by_layers = read_hits(data_path)[9]
    list_hits = []
    for hs in list(hits_by_layers.values()):
        list_hits += hs

    beta_max = math.pi / 100
    print("beta_max:", beta_max)
    m = 7
    costs = get_costs(list_hits, hits_by_layers, beta_max)
    write_costs(costs, costs_path_out, m)
    costs = load_costs(costs_path_out)
    result = run(list_hits, costs, m, model_path_out, solution_path_out, figure_path_out)
    cal_expected_value(list_hits)
