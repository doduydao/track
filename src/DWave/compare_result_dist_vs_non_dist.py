import json
from data import *
import matplotlib.pyplot as plt


def distance(h1, h2):
    distance = math.sqrt((h2.x - h1.x) ** 2 + (h2.y - h1.y) ** 2 + (h2.z - h1.z) ** 2)
    return distance


def calculate_cost_dist_to_non_dist(hits, segs):
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


    tracks = []
    state = [False for _ in segs]

    for i in range(len(segs) - 1):
        track = []
        tmp_seg = segs[i]
        if state[i] == False:
            track.append(segs[i])
            state[i] = True

        for j in range(i + 1, len(segs)):
            if state[j] == False and tmp_seg[1] == segs[j][0]:
                track.append(segs[j])
                state[j] = True
                tmp_seg = segs[j]
        if len(track) > 0:
            tracks.append(track)
    cost = 0
    for track in tracks:
        ct = 0
        for i in range(len(track) - 1):
            for j in range(i + 1, len(track)):
                if track[i][1] == track[j][0]:
                    print(track[i], track[j])
                    h_i = track[i][0]
                    h_j = track[i][1]
                    h_k = track[j][1]

                    ax.plot(xs=[h_i.x, h_j.x], ys=[h_i.y, h_j.y], zs=[h_i.z, h_j.z], color='blue')
                    ax.plot(xs=[h_j.x, h_k.x], ys=[h_j.y, h_k.y], zs=[h_j.z, h_k.z], color='blue')
                    seg_1 = Segment(h_j, h_i)
                    seg_2 = Segment(h_j, h_k)

                    angle = Angle(seg_1=seg_1, seg_2=seg_2).angle
                    # dist = distance(h_i, h_j) + distance(h_j, h_k)
                    dist = 1
                    ct += angle * dist
        print("ct:", ct)
        cost += ct

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    # plt.savefig(out)
    # plt.show()
    return cost


def calculate_cost(segs):
    cost = 0
    for i in range(len(segs) - 1):
        for j in range(i + 1, len(segs)):
            if segs[i][1].hit_id == segs[j][0].hit_id:
                print(segs[i], segs[j])
                h_i = segs[i][0]
                h_j = segs[i][1]
                h_k = segs[j][1]
                seg_1 = Segment(h_j, h_i)
                seg_2 = Segment(h_j, h_k)
                angle = Angle(seg_1=seg_1, seg_2=seg_2).angle
                # dist = distance(h_i, h_j) + distance(h_j, h_k)
                dist = 1
                cost += angle * dist
    return cost


def display(hits, segments):
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
        ax.plot(xs=[h1.x, h2.x], ys=[h1.y, h2.y], zs=[h1.z, h2.z], color='blue')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    # plt.savefig(out)
    # plt.show()


if __name__ == "__main__":
    src_path = '../../src/data_selected/'
    folder = '100hits/'
    data_path = src_path + folder + 'known_track/hits.csv'
    solution_path = "result/" + folder + "known_track/solution_LB_dist_lan_1.json"

    hits = read_hits(data_path)[9]

    with open(solution_path, 'r', encoding='utf-8') as f:
        result = json.load(f)

    layers = [0] + list(hits.keys())
    segments = []
    for var, value in result.items():
        f_p_i_j = var.split('_')
        if 'f' in f_p_i_j[0] and value == 1.0:
            # print(var, value)
            p = int(f_p_i_j[1])
            i = int(f_p_i_j[2])
            j = int(f_p_i_j[3])
            h_1 = hits[layers[p]][i - 1]
            h_2 = hits[layers[p + 1]][j - 1]
            segments.append([h_1, h_2])

    # cost = calculate_cost_dist_to_non_dist(hits, segments)
    cost = calculate_cost(segments)
    print("cost:", cost)
    # display(hits, segments)
