import json
from read_data import *
import random
import matplotlib.pyplot as plt

def display(hits, segments, out=""):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    layers = list(hits.keys())
    xs = []
    ys = []
    zs = []

    for p in layers:
        h_p = hits[p]
        for h in h_p:
            xs.append(h.x)
            ys.append(h.y)
            zs.append(h.z)
    ax.scatter(xs, ys, zs, marker='o', color='red')

    for segment in segments:
        h1 = segment[0]
        h2 = segment[1]
        if h1 in hits[0] or h2 in hits[0]:
            ax.plot(xs=[h1.x, h2.x], ys=[h1.y, h2.y], zs=[h1.z, h2.z], color='green')
        elif h1 in hits[16] or h2 in hits[16]:
            ax.plot(xs=[h1.x, h2.x], ys=[h1.y, h2.y], zs=[h1.z, h2.z], color='black')
        else:
            ax.plot(xs=[h1.x, h2.x], ys=[h1.y, h2.y], zs=[h1.z, h2.z], color='blue')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    plt.savefig(out)
    plt.show()

def create_source_sink(hits):
    layers = sorted(list(hits.keys()))
    first_layer = hits[layers[0]]
    last_layer = hits[layers[-1]]

    source = Hit(
        hit_id=10,
        x=sum([h.x for h in first_layer]) / len(first_layer),
        y=sum([h.y for h in first_layer]) / len(first_layer),
        z=sum([h.z for h in first_layer]) / len(first_layer) - 1,
        volume_id=0,
        layer_id=0,
        module_id=0
    )
    sink = Hit(
        hit_id=10,
        x=sum([h.x for h in last_layer]) / len(last_layer),
        y=sum([h.y for h in last_layer]) / len(last_layer),
        z=sum([h.z for h in last_layer]) / len(last_layer) + 1,
        volume_id=0,
        layer_id=0,
        module_id=0
    )

    return source, sink



if __name__ == '__main__':
    hits_path = '../event000001000/sel/event000001000-hits-sel-01.csv'
    hits_volume = read_hits(hits_path)
    hits = dict()
    for k, v in hits_volume.items():
        print("Volume id:", k)
        print("No_layers:", len(v))
        hits = v

    # for p, hp in hits.items():
    #     z = hp[0].z // 100
    #
    #
    #     print("layer: ", p, "--", )
    #     # print()

    source, sink = create_source_sink(hits)
    hits[0] = [source]
    hits[16] = [sink]
    # print(hits)
    layers = sorted(list(hits.keys()))
    # print(layers)
    model_path_out = "result_f2_Lacomme/model_docplex.lp"
    solution_path_out = "/Users/doduydao/daodd/PycharmProjects/Quantum_Research/Tracking/src/result_f2_Lacomme/result.json"
    with open(solution_path_out, 'r', encoding='utf-8') as f:
        result = json.load(f)['result']

    nt = 6
    M = 10000
    # result = run(hits, nt, M, model_path_out, solution_path_out)

    # del hits[0]
    # del hits[16]
    segments = []
    for var_name, var_value in result.items():
        # print(var)
        # var_name = var['name']
        # var_value = var['value']
        phi_p_p_i_j = var_name.split('_')
        print(var_name, var_value)
        # print(var)
        if 'c' in phi_p_p_i_j[0] or var_value != 1.0:
            continue

        p_1 = int(phi_p_p_i_j[1])
        p_2 = int(phi_p_p_i_j[2])

        # if p_1 == 0 or p_2 == len(layers)-1:
        #     continue

        i = int(phi_p_p_i_j[3])
        j = int(phi_p_p_i_j[4])
        h_1 = hits[layers[p_1]][i - 1]
        h_2 = hits[layers[p_2]][j - 1]
        segments.append([h_1, h_2])
    out = "result_f2_Lacomme/result.PNG"
    display(hits, segments, out)
