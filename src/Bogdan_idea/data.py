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
        self.index = None

    def set_index(self, index):
        self.index = index


class Segment:
    def __init__(self, hit_1, hit_2):
        self.id = (hit_1.index, hit_2.index)
        self.hit_1 = hit_1
        self.hit_2 = hit_2
        self.d_x = hit_2.x - hit_1.x
        self.d_y = hit_2.y - hit_1.y
        self.d_z = hit_2.z - hit_1.z


class Cost:
    def __init__(self, seg_1, seg_2):
        self.id = (seg_1.id, seg_2.id)
        self.seg_1 = seg_1
        self.seg_2 = seg_2
        self.hit_1 = seg_1.hit_1
        self.hit_2 = seg_1.hit_2
        self.hit_3 = seg_2.hit_2
        self.cos_beta = self.calculate_cos_beta()
        self.sum_distance = self.calculate_distance()


    def calculate_cos_beta(self):
        v1 = np.array([self.seg_1.d_x, self.seg_1.d_y, self.seg_1.d_z])
        v2 = np.array([self.seg_2.d_x, self.seg_2.d_y, self.seg_2.d_z])
        dot_product = np.dot(v1, v2)
        norm_v1 = np.linalg.norm(v1)
        norm_v2 = np.linalg.norm(v2)
        return dot_product / (norm_v1 * norm_v2)

    def calculate_distance(self):
        v1 = np.array([self.seg_1.d_x, self.seg_1.d_y, self.seg_1.d_z])
        v2 = np.array([self.seg_2.d_x, self.seg_2.d_y, self.seg_2.d_z])
        norm_v1 = np.linalg.norm(v1)
        norm_v2 = np.linalg.norm(v2)
        return norm_v1 + norm_v2

def read_hits(path):
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
                  layer_id =i[5]/2,
                  module_id = i[6],
                  particle_id= i[7]
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
    folder = '/6hits/'

    data_path = src_path + folder + 'known_track/hits.csv'
    model_path_out = "result" + folder + "known_track/model_docplex_CQM_no_dist_no_LB.lp"
    solution_path = "result" + folder + "known_track/solution_dwave_no_dist_no_LB.json"
    out = "result" + folder + "known_track/result_dwave_no_dist_no_LB.PNG"
    hits = read_hits(data_path)[9]

    track = dict()
    for p, hp in hits.items():
        print(p)
        for h in hp:
            p_id = h.particle_id/10000000000
            if p_id not in track:
                track[p_id] = [h]
            else:
                track[p_id] += [h]

    expect_value = 0
    for pa, hs in track.items():
        t = []
        ct = 0
        m = len(hs)
        for i in range(len(hs) - 2):
            h_i = hs[i]
            h_j = hs[i + 1]
            h_k = hs[i + 2]
            seg_1 = Segment(h_j, h_i)
            seg_2 = Segment(h_j, h_k)
            cost = Cost(seg_1=seg_1, seg_2=seg_2)
            cos_beta = cost.cos_beta
            ct += (-(cost.cos_beta ** m) / cost.sum_distance)
        print('Track', pa, ':', ct)
        expect_value += ct
    print('Expect value of all track is', expect_value)