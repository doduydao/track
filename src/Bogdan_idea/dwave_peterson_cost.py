from data import *
from docplex.mp.model import Model
import json
import matplotlib.pyplot as plt
import neal
from collections import defaultdict
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite
from matplotlib import pyplot as plt
import networkx as nx
from dimod.binary import BinaryQuadraticModel
import re
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
    return x, ob


def to_qubo(ob):
    str_ob = ob.repr_str().replace("+-", "-").replace("-", "+-")
    elements = str_ob.split("+")
    # print("elements:", elements)
    while "" in elements:
        elements.remove("")

    for i in range(len(elements)):
        e = elements[i]
        if "^2" in e:
            match = re.search(r'^[^x]+x', e)
            # print(match)
            if match:
                i_cut = match.end() - 1
                if e[:i_cut].replace('.', '', 1).replace('-', '').isdigit() or e[:i_cut] == "-":
                    elements[i] = e[:i_cut] + e[i_cut:-2] + "*" + e[i_cut:-2]
                else:
                    elements[i] = e[:-2] + "*" + e[:-2]
            else:
                elements[i] = e[:-2] + "*" + e[:-2]

    if 'x' not in elements[-1]:
        offset = float(elements[-1])
    else:
        offset = 0
        elements.append(0)

    qubo = dict()
    for e in elements[:-1]:
        es = e.replace("*", "").split("x")
        if es[0] == "":
            es[0] = "1"
        coeff = float(es[0])
        vars = es[1:]
        if len(vars) > 1:
            key = list()
            for var in vars:
                var = var.split("_")
                i = int(var[1])
                j = int(var[-1])
                key.append((i, j))
            key = tuple(key)

        else:
            var = vars[0].split("_")
            i = int(var[1])
            j = int(var[-1])
            key = ((i, j), (i, j))

        # if offset != 0:
        #     coeff = coeff / offset

        if key not in qubo:
            qubo[key] = coeff
        else:
            qubo[key] += coeff

    return qubo, offset


def create_objective_function(hits_by_layers,list_hits, costs, m, alpha, beta):
    # define model
    model = Model(name="Track")
    model.float_precision = 15
    # create variables
    x, ob = define_variables(model, costs)

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

    ob = -first_part + beta * second_part + alpha * (third_part + fourth_part)
    # print("sum_segments:", sum_segments)
    # print("---" * 10)
    # print("first_part", first_part)
    # print("---" * 10)
    # print("second_part", second_part)
    # print("---" * 10)
    # print("third_part", third_part)
    # print("---" * 10)
    # print("fourth_part", fourth_part)
    # print("---" * 10)

    return ob


def display(list_hits, result, out=""):
    # segments = []
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

    count_segments = 0
    for k, v in result.items():
        if v == 1:
            # print(k, v)
            count_segments += 1
            h1 = list_hits[k[0]]
            h2 = list_hits[k[1]]
            ax.plot(xs=[h1.x, h2.x], ys=[h1.y, h2.y], zs=[h1.z, h2.z], color='blue')
    print("No segments:", count_segments)
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
        hits_l = hits[l]
        first_part = []
        for i in range(max(layers[0], l - 3), l):
            hits_i = hits[i]
            segs_l_i = build_segments(hits_i, hits_l, list_hits)
            first_part.append(segs_l_i)

        second_part = []
        for i in range(l + 1, min(l + 3, layers[-1]) + 1):
            hits_i = hits[i]
            segs_l_i = build_segments(hits_l, hits_i, list_hits)
            second_part.append(segs_l_i)

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
    src_path = '../../src/data_selected'
    folder = '/2hits/known_track/'
    out_path = '/Users/doduydao/daodd/PycharmProjects/track/src/Bogdan_idea/results'
    check_path(out_path + folder)
    data_path = src_path + folder + 'hits.csv'
    costs_path_out = out_path + folder + "costs.json"
    figure_path_out = out_path + folder + "result_dwave.PNG"
    hits_by_layers = read_hits(data_path)[9]
    list_hits = []
    for hs in list(hits_by_layers.values()):
        list_hits += hs

    beta_max = math.pi / 50
    m = 7
    alpha = 1
    beta = 10000

    costs = get_costs(list_hits, hits_by_layers, beta_max)
    write_costs(costs, costs_path_out, m)
    costs = load_costs(costs_path_out)
    ob_funct = create_objective_function(hits_by_layers,list_hits, costs, m, alpha, beta)
    # print("ob_funct:", ob_funct)
    qubo, offset = to_qubo(ob_funct)
    # print(offset)
    # print(qubo)
    # sampler = EmbeddingComposite(DWaveSampler())
    sampler = neal.SimulatedAnnealingSampler()
    response = sampler.sample_qubo(qubo)
    print(response)
    ob_value = response.first.energy
    result = response.first.sample
    # if offset != 0:
    #     ob_value *= offset
    print("ob_value:", ob_value)

    display(list_hits, result, out=figure_path_out)
