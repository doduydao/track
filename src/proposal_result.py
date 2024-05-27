from data import *
from read_data import *
import matplotlib.pyplot as plt


def read_truth_data(path):
    df = pd.read_csv(path)
    return df

def display_bogdan(hits, solution, out=""):
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


def load_hits(path):
    df = pd.read_csv(path)
    list_df = [row.tolist() for index, row in df.iterrows()]
    volumes = dict()

    for i in list_df:
        hit = Hit(
            hit_id=i[0],
            x=i[1],
            y=i[2],
            z=i[3],
            volume_id=i[4],
            layer_id=i[5],
            module_id=i[6],
            particle_id=i[7]
        )
        volume_id = int(hit.volume_id)
        if volume_id not in volumes:
            volumes[volume_id] = [hit]
        else:
            volumes[volume_id] += [hit]
    for id, hits in volumes.items():
        layers = dict()
        for hit in hits:
            layer_id = int(hit.layer_id)
            if layer_id not in layers:
                layers[layer_id] = [hit]
            else:
                layers[layer_id] += [hit]
        volumes[id] = layers
    return volumes

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

    for track, hits in segments.items():
        for i in range(len(hits) - 1):
            h1 = hits[i]
            h2 = hits[i + 1]
            ax.plot(xs=[h1.x, h2.x], ys=[h1.y, h2.y], zs=[h1.z, h2.z], color='blue')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    plt.savefig(out)
    plt.show()


if __name__ == '__main__':
    # hits_path = "../event000001000/volume_id_9/hits-vol_9_20_track.csv"
    # hits_path = "../event000001000/volume_id_9/hits-vol_9_FGC_track_min_6_track_with_noise.csv"
    # hits_path = "../event000001000/volume_id_9/event000001000-volume_id_9.csv"
    hits_path = "../event000001000/volume_id_9/hits-vol_9_FGC_min_20_track_with_noise.csv"

    hits = load_hits(hits_path)[9]

    track = dict()
    for l, hp in hits.items():
        for i in range(len(hp)):
            id = hp[i].hit_id
            particle_id = hp[i].particle_id
            if particle_id not in track:
                track[particle_id] = [hp[i]]
            else:
                track[particle_id] += [hp[i]]
    for t, h in track.items():
        h = sorted(h, key=lambda obj: obj.z)
        track[t] = h
        # if len(h) == 7:
        #     track[t] = h
        # else:
        #     track[t] = []

    out = "../event000001000/volume_id_9/expected_result_9_FGC_min_20_track_with_noise.PNG"
    display(hits, track, out)
