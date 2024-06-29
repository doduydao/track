from data import *
import json
import matplotlib.pyplot as plt
import neal
import os
import time
from pyqubo import Binary
import dimod
from dwave.system import LeapHybridSampler, DWaveSampler, EmbeddingComposite
import dwave.inspector
import multiprocessing
import itertools


def define_variables(costs):
    var = set()
    for ijk in costs.keys():
        i, j, k = ijk[0], ijk[1], ijk[2]
        var.add((i, j))
        var.add((j, k))
    var = sorted(var, key=lambda x: (x[0], x[1]))
    x = dict()
    for v in var:
        x[v] = Binary('x_' + str(v[0]) + "_" + str(v[1]))
    return x


def create_objective_function(costs, hits, hits_by_layers, alpha, gamma):
    x = define_variables(costs)
    NL = len(hits_by_layers.keys())

    first_part = sum(-alpha * w * x[(ijk[0], ijk[1])] * x[(ijk[1], ijk[2])] for ijk, w in costs.items())

    third_part = 0
    fourth_part = 0
    for h in hits:
        k = h.index
        constraint_out = []
        constraint_in = []
        for k_1 in x.keys():
            if ((h.layer_id < (NL - 1)) and (k_1[0] == k)):
                constraint_out.append(x[(k_1[0], k_1[1])])
            if ((h.layer_id > 1) and (k_1[1] == k)):
                constraint_in.append(x[(k_1[0], k_1[1])])
        if (len(constraint_out) > 0):
            third_part += (1 - sum(constraint_out)) ** 2
        if (len(constraint_in) > 0):
            fourth_part += (1 - sum(constraint_out)) ** 2
    H = first_part + gamma * third_part + gamma * fourth_part
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
    print("No. segments by model :", len(segments))
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


def build_segments(hits_1, hits_2):
    segments = []
    for i in hits_1:
        for j in hits_2:
            segments.append(Segment(i, j))
    return segments


def compute_cost(first_part, second_part, cos_beta_max):
    costs = []
    for seg_f in first_part:
        for seg_s in second_part:
            if seg_f.hit_2.hit_id == seg_s.hit_1.hit_id and seg_f.gaps + seg_s.gaps <= 4:
                cost = Cost(seg_f, seg_s)
                if cost.cos_beta >= cos_beta_max:
                    costs.append(cost)
    return costs


def get_costs(hits, beta_max):
    layers = list(hits.keys())
    costs = []

    cos_beta_max = math.cos(beta_max)
    for l in layers[1:-1]:
        first_part = [build_segments(hits[i], hits[l]) for i in range(max(layers[0], l - 3), l)]
        second_part = [build_segments(hits[l], hits[k]) for k in range(l + 1, min(l + 3, layers[-1]) + 1)]
        arg_triples = list(e + (cos_beta_max,) for e in list(itertools.product(first_part, second_part)))
        with multiprocessing.Pool() as pool:
            result = pool.starmap(compute_cost, arg_triples)
            costs += [item for sublist in result for item in sublist]
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


def cal_expected_value(list_hits, m):
    track = dict()
    for hit in list_hits:
        k = hit.particle_id / 1000
        if k not in track:
            track[k] = [hit]
        else:
            track[k].append(hit)
    cost = 0
    segs = set()
    for t, hs in track.items():
        for i in range(len(hs) - 2):
            h_i = hs[i]
            h_j = hs[i + 1]
            h_k = hs[i + 2]
            seg_1 = Segment(h_i, h_j)
            seg_2 = Segment(h_j, h_k)
            segs.add(seg_1.id)
            segs.add(seg_2.id)
            c = Cost(seg_1, seg_2)
            cost -= c.cos_beta ** m / c.sum_distance

    print("Expected cost :", cost)
    print("Expected No. segments :", len(segs))


def simulate_annealing(H, num_reads):
    sampler = neal.SimulatedAnnealingSampler()
    qubo, offset = H.to_qubo()
    print("offset:", offset)

    start = time.time()
    sampleset = sampler.sample_qubo(qubo, num_reads=num_reads)
    end = time.time()
    print("Execution time:", end - start)
    decoded_samples = H.decode_sampleset(sampleset)
    best_sample = min(decoded_samples, key=lambda x: x.energy)
    ob_value = best_sample.energy
    result = best_sample.sample
    return ob_value, result



if __name__ == '__main__':
    src_path = '../../src/data_selected'
    folder = '/2hits/known_track/'

    out_path = 'results'
    check_path(out_path + folder)
    data_path = src_path + folder + 'hits.csv'
    figure_path_out = out_path + folder + "result_D_AQUBM.PNG"
    # read data
    print("Loading data...")
    start = time.time()
    hits_by_layers, hits = read_hits(data_path)
    end = time.time()
    print("Loaded data! Execution time: ", end - start)
    # weight
    A = 50
    beta_max = math.pi / A
    m = 1

    re_calculate = False
    costs_path_out = out_path + folder + "pi_" + str(A) + "costs.json"
    if os.path.exists(costs_path_out) == False:
        re_calculate = True
    if re_calculate:
        # calculate costs
        print("\n----Compute cost----")
        start = time.time()
        costs = get_costs(hits_by_layers, beta_max)
        end = time.time()
        print('Complete!. Execution time: ', end - start, 's')

        print("\n---Write cost out---")
        print("Path: ", costs_path_out)
        write_costs(costs, costs_path_out, m)

    print("---Load cost---")
    costs = load_costs(costs_path_out)

    alpha = 100
    gamma = 1

    H = create_objective_function(costs, hits, hits_by_layers, alpha, gamma)

    cal_expected_value(hits, m)
    print("\n----Simulating----")
    num_reads = 100
    ob_value, result = simulate_annealing(H, num_reads)
    print("Objective value:", ob_value / alpha)
    display(hits, result, out=figure_path_out)
