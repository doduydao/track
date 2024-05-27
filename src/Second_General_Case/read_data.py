import pandas as pd
from data import Hit
import random

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
                  layer_id =i[5],
                  module_id = i[6]
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

def select_hits_from_sublayer(hits, sublayer_selected, no_hits):
    sublayer_dict = dict()
    for hit in hits:  # Hp : tập các hit có trên cùng 1 layer
        if hit.z not in sublayer_dict:
            sublayer_dict[hit.z] = [hit]
        else:
            sublayer_dict[hit.z] += [hit]
    sublayer_dict = [i[1] for i in sorted(sublayer_dict.items(), key=lambda x: x[0])]
    hits = sublayer_dict[sublayer_selected]

    idx = set()
    if no_hits != -1:
        while len(idx) < no_hits:
            idx.add(random.randint(0, len(hits)-1))

        hits_random = []
        for i in idx:
            hits_random.append(hits[i])
        return hits_random
    else:
        return hits

def select_hits_from_volume(hits, volume_id):
    return hits[volume_id]


def write_out_selected_hits(hits_volume, pathout):

    hit_id = []
    x = []
    y = []
    z = []
    volume_id = []
    layer_id = []
    module_id = []

    for layer, hits in hits_volume.items():
        for hit in hits:
            hit_id.append(hit.hit_id)
            x.append(hit.x)
            y.append(hit.y)
            z.append(hit.z)
            volume_id.append(hit.volume_id)
            layer_id.append(hit.layer_id)
            module_id.append(hit.module_id)
    print(len(hit_id))
    data_dict = {"hit_id": hit_id,
                 "x":x,
                 "y":y,
                 "z":z,
                 "volume_id":volume_id,
                 "layer_id":layer_id,
                 "module_id":module_id}
    df = pd.DataFrame.from_dict(data_dict)
    df.to_csv(pathout, index=False, sep=',')
    print("Done!")


def decode_particle_id(data):
    """Decode particle_id into vertex id, generation, etc.
    """
    components = [
        ('vertex_id',    0xfff0000000000000, 13 * 4),
        ('primary_id',   0x000ffff000000000, 9 * 4),
        ('generation',   0x0000000fff000000, 6 * 4),
        ('secondary_id', 0x0000000000fff000, 3 * 4),
        ('process',      0x0000000000000fff, 0),
    ]
    pid = data['particle_id'].values.astype('u8')
    for name, mask, shift in components:
        data[name] = (pid & mask) >> shift
    return data

# def select_from_truth_hits(truth_path, no_truth_hits, hits_volume):
#     df = pd.read_csv(truth_path)
#     df = decode_particle_id(df)
#     list_df = [row.tolist() for index, row in df.iterrows()]
#     volumes = dict()
#
#     all_hits = []
#     for i in list_df:
#         hit = Hit(
#             hit_id=i[0],
#             particle_id = i[1],
#             x=i[2],
#             y=i[3],
#             z=i[4]
#         )
#         all_hits.append(hit)
#
#     for p, hp in hits_volume.items():
#         for h in hp:
#
#
#
#
#     #     volume_id = int(hit.volume_id)
#     #     if volume_id not in volumes:
#     #         volumes[volume_id] = [hit]
#     #     else:
#     #         volumes[volume_id] += [hit]
#     # for id, hits in volumes.items():
#     #     layers = dict()
#     #     for hit in hits:
#     #         layer_id = int(hit.layer_id)
#     #         if layer_id not in layers:
#     #             layers[layer_id] = [hit]
#     #         else:
#     #             layers[layer_id] += [hit]
#     #     volumes[id] = layers
#
#     return volumes

if __name__ == '__main__':
    hits_path = '/Users/doduydao/daodd/PycharmProjects/Quantum_Research/Tracking/event000001000/event000001000-hits.csv'
    all_hits = read_hits(hits_path)

    volume_id = 9
    no_hits = -1
    sublayer_selected = 2
    hits_volume = select_hits_from_volume(all_hits, volume_id)

    for k, v in hits_volume.items():
        hits_layer = select_hits_from_sublayer(v, sublayer_selected, no_hits)
        hits_volume[k] = hits_layer
    # for k, v in hits_volume.items():
    #     print(k, len(v))


    # pathout = "/Users/doduydao/daodd/PycharmProjects/Quantum_Research/Tracking/event000001000/sublayer_2/event000001000-hits_volume_9.csv"
    # write_out_selected_hits(hits_volume, pathout)

    # truth_path = "/Users/doduydao/daodd/PycharmProjects/Quantum_Research/Tracking/event000001000/event000001000-truth.csv"
    # no_truth_hits = 10
    # select_from_truth_hits(truth_path, no_truth_hits, hits_volume)




