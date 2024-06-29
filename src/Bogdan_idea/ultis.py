from data import *
import json
import os
import time
import multiprocessing
import itertools
from tqdm import tqdm
import math


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
    cos_beta_max = math.cos(beta_max)
    layers = list(hits.keys())
    L =  layers[1:-1]
    costs = []
    
    with multiprocessing.Pool() as pool:
        for i in tqdm (range (len(L)), 
               desc="Computing…", 
               ascii=False, ncols=75):
            l = L[i]
            first_part = [build_segments(hits[i], hits[l]) for i in range(max(layers[0], l - 3), l)]
            second_part = [build_segments(hits[l], hits[k]) for k in range(l + 1, min(l + 3, layers[-1]) + 1)]
            arg_triples = list(e + (cos_beta_max,) for e in list(itertools.product(first_part, second_part)))
            result = pool.starmap(compute_cost, arg_triples)
            costs += [item for sublist in result for item in sublist]
        time.sleep(0.01)
    
    # all_segments = set()
    # for cost in costs:
    #     all_segments.add(cost.seg_1.id)
    #     all_segments.add(cost.seg_2.id)
    # print("\nNumber of segments:", len(all_segments))

    return costs


def write_costs(costs, path, m):
    data = dict()
    for cost in costs:
        i_j = cost.id[0]
        j_k = cost.id[1]
        i = i_j[0]
        j = i_j[1]
        k = j_k[1]
        cos_beta = cost.cos_beta
        sum_distance = cost.sum_distance
        str_key = str((i, j, k))
        data[str_key] = (cos_beta ** m) / sum_distance
    with open(path, 'w') as f:
        json.dump(data, f, ensure_ascii=False, indent=4)
    print('Wrote cost!')


def load_costs(path):
    with open(path) as f:
        data = json.load(f)
    costs = dict()
    for k, v in data.items():
        costs[eval(k)] = v
    return costs


def check_path(path):
    if os.path.exists(path) == False:
        print("Path is not exist")
        os.mkdir(path)
        print("Created path!")


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


# if __name__ == '__main__':
#     src_path = '../../src/data_selected'
#     folder = '/2hits/known_track/'
#     out_path = '/Users/doduydao/daodd/PycharmProjects/track/src/Bogdan_idea/results'
#     data_path = src_path + folder + 'hits.csv'
#     figure_path_out = out_path + folder + "result_D_AQUBM.PNG"
#     check_path(out_path + folder)

#     # read data
#     print("Loading data...")
#     start = time.time()
#     hits_by_layers, hits = read_hits(data_path)
#     end = time.time()
#     print("Loaded data! Execution time: ", end - start)
#     # weight
#     A = 50
#     beta_max = math.pi / A
#     m = 1

#     re_calculate = False
#     costs_path_out = out_path + folder + "pi_" + str(A) + "costs.json"
#     if os.path.exists(costs_path_out) == False:
#         re_calculate = True
#     if re_calculate:
#         # calculate costs
#         print("\n----Compute cost----")
#         start = time.time()
#         costs = get_costs(hits_by_layers, beta_max)
#         end = time.time()
#         print('Complete!. Execution time: ', end - start, 's')

#         print("\n---Write cost out---")
#         print("Path: ", costs_path_out)
#         write_costs(costs, costs_path_out, m)

#     print("---Load cost---")
#     costs = load_costs(costs_path_out)