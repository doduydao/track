from data import *
import json
import matplotlib.pyplot as plt
import neal
import os
import time
from pyqubo import Binary

def define_variables(costs):
    var = set()
    for cost in costs.keys():
        i_j = cost[0]
        j_k = cost[1]
        var.add(i_j)
        var.add(j_k)
    var = sorted(var, key=lambda x: (x[0], x[1]))
    x = dict()
    for v in var:
        x[v] = Binary('x_' + str(v[0]) + "_" + str(v[1]))
    return x


def create_objective_function(list_hits, costs, m, alpha, beta):
    x = define_variables(costs)
    # print(x)
    hit_last_layer = hits_by_layers[m]
    N = len(list_hits) - len(hit_last_layer)
    print("N =", N)

    # first part
    first_part = 0
    segments = set()
    for id, cost in costs.items():
        i_j = id[0]
        j_k = id[1]
        first_part += cost * x[i_j] * x[j_k]
        segments.add(i_j)
        segments.add(j_k)

    sum_segments = sum([x[s] for s in segments])
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

    H = -first_part + alpha * second_part + beta * (third_part + fourth_part)
    return H


def display(list_hits, result, out=""):
    segments = list()
    for k, v in result.items():
        if v == 1:
            if type(k) is str:
                k = [int(e) for e in k.split('_')[1:]]
            h_1 = list_hits[k[0]]
            h_2 = list_hits[k[1]]
            segments.append([h_1, h_2])
    print("No_segments:", len(segments))
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    xs = []
    ys = []
    zs = []

    for h in list_hits:
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
    src_path = '../../src/data_selected'
    folder = '/2hits/known_track/'
    out_path = '/Users/doduydao/daodd/PycharmProjects/track/src/Bogdan_idea/results'
    check_path(out_path + folder)
    data_path = src_path + folder + 'hits.csv'
    costs_path_out = out_path + folder + "costs.json"
    figure_path_out = out_path + folder + "result_qubo.PNG"
    hits_by_layers = read_hits(data_path)[9]
    list_hits = []

    for hs in list(hits_by_layers.values()):
        list_hits += hs

    beta_max = math.pi / 2
    m = 7
    costs = get_costs(list_hits, hits_by_layers, beta_max)
    write_costs(costs, costs_path_out, m)
    costs = load_costs(costs_path_out)

    alpha = 1
    beta = 1

    H = create_objective_function(list_hits, costs, m, alpha, beta)
    # print(H)
    model = H.compile()
    qubo, offset = model.to_qubo()
    print("offset:", offset)
    print(len(qubo.keys()))
    # qubo = {(k[1], k[0]): v for k, v in sorted(qubo.items(), key=lambda x:x[0][1])}
    # for k, v in qubo.items():
    #     print(k, v)

    sampler = neal.SimulatedAnnealingSampler()

    start = time.time()
    response = sampler.sample_qubo(qubo)
    print(response)
    ob_value = response.first.energy
    result = response.first.sample
    print("ob_value:", ob_value)
    end = time.time()
    display(list_hits, result, out=figure_path_out)
    print("time:", end - start)

    # obf = cbf(hits_by_layers, list_hits, costs, m, alpha, beta)
    # qubo_dao, offset_dao = tqubo(obf)
    #
    # print("offset_dao:", offset_dao)
    # print(len(qubo_dao.keys()))
    # qubo_tmp = dict()
    # for k, v in qubo_dao.items():
    #     x_i_j = "x_" + str(k[0][0]) + "_" + str(k[0][1])
    #     x_j_k = "x_" + str(k[1][0]) + "_" + str(k[1][1])
    #     qubo_tmp[(x_i_j, x_j_k)] = v
    # for k, v in qubo_tmp.items():
    #     # print(k, v)
    #     if k in qubo:
    #         if round(v,4) != round(qubo[k],4):
    #             print(k, v, qubo[k])
    #     else:
    #         print(k, v)

    # print(qubo_dao)

    # start = time.time()
    # response = sampler.sample_qubo(qubo_dao)
    # # # response = sampler.sample(bqm)
    # print(response)
    # ob_value = response.first.energy
    # result = response.first.sample
    # print("ob_value:", ob_value)
    # end = time.time()
    # display(list_hits, result, out=figure_path_out)
    # print("time:", end - start)


