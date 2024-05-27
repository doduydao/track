from data import *
from read_data import *
import math
import matplotlib.pyplot as plt


def calculate_all_angles(hits):
    # tính tất cả các góc nhỏ nhất theo layer
    angles = dict()
    layers = list(hits.keys())
    L = len(layers) + 1
    layers = [0] + layers + [0]

    betas = []

    for p in range(1, L - 2):
        min_angle = 100000000
        n_p = len(hits[layers[p]]) + 1
        n_p_1 = len(hits[layers[p + 1]]) + 1
        n_p_2 = len(hits[layers[p + 2]]) + 1

        for i in range(1, n_p):
            for j in range(1, n_p_1):
                for k in range(1, n_p_2):
                    h_i = hits[layers[p]][i - 1]
                    h_j = hits[layers[p + 1]][j - 1]
                    h_k = hits[layers[p + 2]][k - 1]
                    seg_1 = Segment(h_j, h_i)
                    seg_2 = Segment(h_j, h_k)
                    angle = Angle(seg_1=seg_1, seg_2=seg_2).angle * (distance(h_i, h_j) + distance(h_j, h_k))
                    # angle = Angle(seg_1=seg_1, seg_2=seg_2).angle
                    betas.append(angle)
                    if min_angle > angle:
                        min_angle = angle
        angles[layers[p]] = min_angle * len(hits[layers[p]])
    print(len(betas))
    max_beta = max(betas)
    mean_beta = sum(betas) / len(betas)
    min_beta = min(betas)
    print("Max beta:", max_beta)
    print("Min beta:", min_beta)
    print("Mean beta:", mean_beta)
    print("abc:", (max_beta + mean_beta) / 2)
    dict_betas = dict()
    interval = (max(betas) - min(betas)) / 500
    for i in range(math.ceil(math.pi / interval)):
        dict_betas[i * interval] = 0

    for beta in betas:
        if beta <= mean_beta:
            key = math.floor(beta / interval) * interval

            if key not in dict_betas:
                dict_betas[key] = 1
            else:
                dict_betas[key] += 1

    fre_beta_sorted = sorted(dict_betas.items(), key=lambda x: x[0])
    for i in range(len(fre_beta_sorted) - 1):
        interval = str(round(fre_beta_sorted[i][0])) + "-" + str(round(fre_beta_sorted[i + 1][0]))
        print(interval, fre_beta_sorted[i][1])

    # with open(out, 'w') as file:
    #     json.dump(angles, file)
    # return angles


def distance(h1, h2):
    distance = math.sqrt((h2.x - h1.x) ** 2 + (h2.y - h1.y) ** 2 + (h2.z - h1.z) ** 2)
    return distance


def create_source_sink(hits):
    layers = sorted(list(hits.keys()))
    first_layer = hits[layers[0]]
    last_layer = hits[layers[-1]]

    source = Hit(
        hit_id=10,
        x=sum([h.x for h in first_layer]) / len(first_layer),
        y=sum([h.y for h in first_layer]) / len(first_layer),
        z=sum([h.z for h in first_layer]) / len(first_layer) - 200,
        volume_id=0,
        layer_id=0,
        module_id=0
    )

    sink = Hit(
        hit_id=10,
        x=sum([h.x for h in last_layer]) / len(last_layer),
        y=sum([h.y for h in last_layer]) / len(last_layer),
        z=sum([h.z for h in last_layer]) / len(last_layer) + 200,
        volume_id=0,
        layer_id=0,
        module_id=0
    )

    return source, sink




if __name__ == '__main__':
    src_path = 'data_selected'
    data_path = src_path + '/15hits/unknown_track/hits.csv'
    hits_volume = read_hits(data_path)
    hits = hits_volume[9]
    source, sink = create_source_sink(hits)
    hits[0] = [source]
    hits[16] = [sink]
    # out_angles_path = "angles.json"
    angles = calculate_all_angles(hits)
