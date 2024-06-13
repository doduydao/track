import numpy as np
import math
import pandas as pd


class Hit:
    def __init__(self, hit_id, particle_id=None, x=0, y=0, z=0, volume_id=None, layer_id=None, module_id=None,
                 selected=None):
        self.hit_id = hit_id
        self.x = x
        self.y = y
        self.z = z
        self.particle_id = particle_id
        self.volume_id = volume_id
        self.layer_id = layer_id
        self.module_id = module_id
        self.selected = selected
        self.track = None


class Segment:
    def __init__(self, hit_1, hit_2):
        self.x = hit_2.x - hit_1.x
        self.y = hit_2.y - hit_1.y
        self.z = hit_2.z - hit_1.z


class Angle:
    def __init__(self, seg_1, seg_2):
        self.angle = self.calculate_angle(seg_1, seg_2)

    def calculate_angle(self, seg_1, seg_2):
        v1 = np.array([seg_1.x, seg_1.y, seg_1.z])
        v2 = np.array([seg_2.x, seg_2.y, seg_2.z])

        dot_product = np.dot(v1, v2)
        magnitude1 = np.linalg.norm(v1)
        magnitude2 = np.linalg.norm(v2)
        return math.pi - np.arccos(dot_product / (magnitude1 * magnitude2))


def read_hits(path):
    df = pd.read_csv(path)
    # print(df)
    list_df = [row.tolist() for index, row in df.iterrows()]
    volumes = dict()

    for i in list_df:
        hit = Hit(
            hit_id=i[0],
            x=i[1],
            y=i[2],
            z=i[3],
            volume_id=i[4],
            layer_id=i[5] / 2,
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


if __name__ == '__main__':
    src_path = '../../src/data_selected'
    folder = '/50hits/'

    data_path = src_path + folder + 'known_track/hits.csv'
    model_path_out = "result" + folder + "known_track/model_docplex_CQM_no_dist_no_LB.lp"
    solution_path = "result" + folder + "known_track/solution_dwave_no_dist_no_LB.json"
    out = "result" + folder + "known_track/result_dwave_no_dist_no_LB.PNG"
    hits = read_hits(data_path)[9]

    track = dict()
    for p, hp in hits.items():
        print(p)
        for h in hp:
            p_id = h.particle_id / 10000000000
            if p_id not in track:
                track[p_id] = [h]
            else:
                track[p_id] += [h]

    cost = 0
    count = 1
    for pa, hs in track.items():
        t = []
        ct = 0
        for i in range(len(hs) - 2):
            h_i = hs[i]
            h_j = hs[i + 1]
            h_k = hs[i + 2]
            seg_1 = Segment(h_j, h_i)
            seg_2 = Segment(h_j, h_k)
            angle = Angle(seg_1=seg_1, seg_2=seg_2).angle
            beta = angle
            ct += beta
        # print('Track', count, ':', ct)
        cost += ct
    print('cost of all track is', cost)
