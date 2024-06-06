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


def define_variables(model, costs):
    var = set()
    for cost in costs:
        i_j = cost.id[0]
        j_k = cost.id[1]
        var.add(i_j)
        var.add(j_k)
    var = sorted(var, key=lambda x: (x[0], x[1]))
    x = model.binary_var_dict(var, name='x')
    ob = model.continuous_var(name="ob")
    return x, ob


def to_qubo(ob):
    str_ob = ob.repr_str().replace("+-", "-").replace("-", "+-")
    elements = str_ob.split("+")
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

    qubo = dict()
    # print(elements[:-1])
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

        if offset != 0:
            coeff = coeff / offset
        # if key == ((12, 21), (21, 36)):
        #     print("error:",e)
        if key not in qubo:
            qubo[key] = coeff
        else:
            qubo[key] += coeff
    #
    # for k, v in qubo.items():
    #     print(k, v)

    return qubo, offset


def create_objective_function(hits_by_layers, list_hits, costs, m, M, P_1, P_2, P_3):
    # define model
    model = Model(name="Track")
    model.float_precision = 15

    # create variables
    x, ob = define_variables(model, costs)

    # create objective function
    hit_last_layer = hits_by_layers[m]
    N = len(list_hits) - len(hit_last_layer)

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

    third_part = 0
    for k in list(x.keys()):
        i = k[0]
        tmp = set()
        for k_1 in list(x.keys()):
            if k_1[0] == i:
                tmp.add(k_1)
        tmp = list(tmp)
        for m in range(len(tmp) - 1):
            for n in range(m + 1, len(tmp)):
                # print("abc:", tmp[m], tmp[n])
                third_part += x[tmp[m]] * x[tmp[n]]

    fourth_part = 0
    for k in list(x.keys()):
        j = k[1]
        tmp = set()
        for k_1 in list(x.keys()):
            if k_1[1] == j:
                tmp.add(k_1)
        tmp = list(tmp)
        for m in range(len(tmp) - 1):
            for n in range(m + 1, len(tmp)):
                # print("abc:", tmp[m], tmp[n])
                fourth_part += x[tmp[m]] * x[tmp[n]]

    ob = first_part + second_part + P_1 * third_part + P_2 * fourth_part

    return 1 / 2 * ob


def display(list_hits, result, out=""):
    segments = []
    for k, v in result.items():
        if v == 1:
            print(k)
            h_1 = list_hits[k[0]]
            h_2 = list_hits[k[1]]
            segments.append([h_1, h_2])

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
                            print("cos_beta:", cos_beta)
                            if cos_beta >= math.cos(math.pi/20):
                                costs.append(cost)
    all_segments = set()
    for cost in costs:
        all_segments.add(cost.seg_1)
        all_segments.add(cost.seg_2)
    print("number of segments:", len(all_segments))
    return costs


if __name__ == '__main__':
    src_path = '../../src/data_selected'
    folder = '/6hits/'
    data_path = src_path + folder + 'known_track/hits.csv'
    hits_by_layers = read_hits(data_path)[9]

    list_hits = []
    for hs in list(hits_by_layers.values()):
        list_hits += hs
    costs = get_costs(list_hits, hits_by_layers)

    m = len(hits_by_layers.keys())
    M = 1
    P_1 = 1
    P_2 = 1
    P_3 = 1
    ob_funct = create_objective_function(hits_by_layers, list_hits, costs, m, M, P_1, P_2, P_3)
    qubo, offset = to_qubo(ob_funct)

    # for k, v in qubo.items():
    #     print(k, v)

    # sampler = EmbeddingComposite(DWaveSampler())
    sampler = neal.SimulatedAnnealingSampler()
    response = sampler.sample_qubo(qubo)

    ob_value = response.first.energy
    result = response.first.sample

    if offset != 0:
        ob_value *= offset
    print("ob_value:", ob_value)
    figure_out = 'results' + folder + 'known_track/result_dwave.PNG'
    display(list_hits, result, out=figure_out)
