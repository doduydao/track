from data import *
from read_data import *
import matplotlib.pyplot as plt


def distance_to_line(src, dest, point):
    vec_src_dest = np.array([dest.x - src.x, dest.y - src.y, dest.z - src.z])
    vec_point_src = np.array([point.x - src.x, point.y - src.y, point.z - src.z])
    vec_project = np.cross(vec_src_dest, vec_point_src)

    return np.linalg.norm(vec_project)


def considering_track(hits):
    tracks = dict()
    layers = sorted(list(hits.keys()))
    conflict_hits = []
    for p in layers:
        hp = hits[p]
        for h in hp:
            if h.track is not None:
                for t in h.track:
                    if t not in tracks:
                        tracks[t] = [h]
                    else:
                        tracks[t] += [h]

                if len(h.track) > 1:
                    conflict_hits.append(h)

    for t, ht in tracks.items():
        tracks[t] = sorted(ht, key=lambda hit: hit.layer_id)

    for h in conflict_hits:
        ts = h.track
        print(h, h.track)
        min_beta = 100
        best_track = 0
        for t in ts:
            print(t, tracks[t])
            ht = tracks[t]
            id_h = ht.index(h)
            h_i = ht[id_h - 1]
            h_j = h
            h_k = ht[id_h + 1]
            seg_1 = Segment(h_j, h_i)
            seg_2 = Segment(h_j, h_k)
            beta = Angle(seg_1=seg_1, seg_2=seg_2).angle
            if beta < min_beta:
                min_beta = beta
                best_track = t

        h.track = best_track
        # print(best_track)
        ts.remove(best_track)

        for t in ts:
            tracks[t] = tracks[t].remove(h)

    return tracks


def run(hits):
    layers = sorted(list(hits.keys()))

    tracks = []
    track_name = 1
    for i in hits[2]:
        track = None
        min_cost_t_i = 1000000000000000
        for j in hits[14]:
            track_i_j = []
            for p in layers[1:-1]:
                min_dist = 1000000000000000
                best_hit = None
                for h in hits[p]:
                    dist = distance_to_line(i, j, h)
                    if min_dist > dist:
                        best_hit = [dist, h]
                        min_dist = dist
                track_i_j.append(best_hit)
            cost_track = sum([h[0] for h in track_i_j]) ** 2
            if min_cost_t_i > cost_track:
                hs = [i] + [h[1] for h in track_i_j] + [j]
                track = [cost_track, hs]
                min_cost_t_i = cost_track
        hs = track[1]
        state = False
        for hit in hs:
            if hit.track is None:
                hit.track = [track_name]
                state = True
            else:
                state = False
            #     hit.track += [track_name]
        if state == True:
            tracks.append(track)
            track_name += 1

    total_cost = sum([t[0] for t in tracks])
    tracks = [t[1] for t in tracks]
    # tracks = considering_track(hits)
    return total_cost, tracks


def display(hits, tracks):
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

    for hits in tracks:
        for i in range(len(hits) - 1):
            h1 = hits[i]
            h2 = hits[i + 1]
            ax.plot(xs=[h1.x, h2.x], ys=[h1.y, h2.y], zs=[h1.z, h2.z], color='blue')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    # plt.savefig(out)
    plt.show()


if __name__ == '__main__':
    hits_path = '../event000001000/volume_id_9/hits-vol_9_FGC_min_20_track_with_noise.csv'
    # hits_path = '../event000001000/volume_id_9/hits-vol_9_FGC_track_min_6_track_with_noise.csv'
    # hits_path = '../event000001000/volume_id_9/hits-vol_9_131_track.csv'
    # hits_path = '../event000001000/volume_id_9/hits-vol_9_FGC_min_20_track_with_noise.csv'
    hits_volume = read_hits(hits_path)
    hits = dict()
    for k, v in hits_volume.items():
        print("Volume id:", k)
        print("No_layers:", len(v))
        hits = v
    for p, hp in hits.items():
        print(p, len(hp))
    import time

    start = time.time()
    total_cost, tracks = run(hits)
    end = time.time()
    print("Execution time:", end - start)
    print("cost:", total_cost)
    display(hits, tracks)
