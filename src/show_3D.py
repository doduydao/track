import matplotlib.pyplot as plt
import numpy as np
from read_data import *
import json


def display(hits, solution, out=""):
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

    for t, h in solution.items():
        for i in range(len(h) - 1):
            # ax.plot_lines(xdata=[], ydata=[ ], zdata=[], color='blue')  # Adjust color as desired
            ax.plot(xs=[h[i].x, h[i + 1].x], ys=[h[i].y, h[i + 1].y], zs=[h[i].z, h[i + 1].z], color='blue')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    plt.savefig(out)
    plt.show()



if __name__ == '__main__':
    # hits_path = '/Users/doduydao/daodd/PycharmProjects/Quantum_Research/Tracking/event000001000/event000001000-hits.csv'
    # hits_path = "/Users/doduydao/daodd/PycharmProjects/Quantum_Research/Tracking/event000001000/sel/event000001000-hits-sel-01.csv"
    # hits_path = 'C:/Users/dddo/PycharmProjects/Quantum_Research/Tracking/event000001000/sel/event000001000-hits-sel-01.csv'
    hits_path = '/Users/doduydao/daodd/PycharmProjects/Quantum_Research/Tracking/event000001000/volume_id_9/hits-vol_9_20_track.csv'

    hits_volume = read_hits(hits_path)
    hits = dict()
    for k, v in hits_volume.items():
        # print(k, v)
        print("Volume id:", k)
        print("No_layers:", len(v))
        hits = v

    # layers = list(hits.keys())
    # f = open('result_pulp_Bogdan_full.json')
    # result = json.load(f)
    # f.close()
    #
    # solution = dict()
    # for var_name, var_value in result.items():
    #     x_p_t_i = var_name.split('_')
    #     if x_p_t_i[0] == 'y' or var_value == 0:
    #         continue
    #     p = int(x_p_t_i[1])
    #     t = int(x_p_t_i[2])
    #     i = int(x_p_t_i[3])
    #     print(var_name, var_value)
    #     if t not in solution:
    #         solution[t] = [hits[layers[p - 1]][i - 1]]
    #     else:
    #         solution[t] += [hits[layers[p - 1]][i - 1]]

    # out = "data.PNG"
    solution = dict()
    out = ""
    # out = "/Users/doduydao/daodd/PycharmProjects/Quantum_Research/Tracking/event000001000/sublayer_2/result_Bogdan_2.PNG"
    display(hits, solution, out)
    print("Dao")