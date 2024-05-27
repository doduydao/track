import numpy as np
import math


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
