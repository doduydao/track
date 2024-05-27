import matplotlib
import networkx as nx
import dimod
from dimod import ConstrainedQuadraticModel, Binary, quicksum
from docplex.mp.model import Model
from dwave.system import LeapHybridCQMSampler
import json

from ..data import *
from ..read_data import *

try:
    import matplotlib.pyplot as plt
except ImportError:
    matplotlib.use("agg")
    import matplotlib.pyplot as plt


def distance(h1, h2):
    distance = math.sqrt((h2.x - h1.x) ** 2 + (h2.y - h1.y) ** 2 + (h2.z - h1.z) ** 2)
    return distance


def build_model(hits, model_path_out):
    model = Model(name="Track")

    layers = list(hits.keys())
    print("layers:", layers)
    no_layer = len(layers)
    no_hits = len(list(hits.values())[0])
    print(no_layer, no_hits)

    f = model.binary_var_dict(
        [(p, i, j) for p in range(1, no_layer) for i in range(1, no_hits + 1) for j in range(1, no_hits + 1)],
        name="f")

    objective = 0
    for p in range(1, no_layer - 1):
        for i in range(1, no_hits + 1):
            for j in range(1, no_hits + 1):
                for k in range(1, no_hits + 1):
                    h_i = hits[layers[p - 1]][i - 1]
                    h_j = hits[layers[p]][j - 1]
                    h_k = hits[layers[p + 1]][k - 1]
                    seg_1 = Segment(h_j, h_i)
                    seg_2 = Segment(h_j, h_k)
                    objective += f[p, i, j] * f[p+1, j, k] * Angle(seg_1=seg_1, seg_2=seg_2).angle

    model.set_objective('min', objective)

    # Constraints
    # first constraints:
    print("---First constraints---")
    count_constraint = 0

    for j in range(1, no_hits + 1):
        for p in range(1, no_layer):
            tmp = 0
            for i in range(1, no_hits + 1):
                tmp += f[p, i, j]
            constraint_name = "FC_" + str(count_constraint)
            count_constraint += 1
            model.add_constraint(tmp == 1, ctname=constraint_name)
    print("Number of first constraints:", count_constraint)
    # print()
    # Second constraints:
    print("---Second constraints---")
    count_constraint = 0
    for i in range(1, no_hits + 1):
        for p in range(1, no_layer):
            tmp = 0
            for j in range(1, no_hits + 1):
                tmp += f[p, i, j]
            constraint_name = "SC_" + str(count_constraint)
            count_constraint += 1
            model.add_constraint(tmp == 1, ctname=constraint_name)
    print("Number of second constraints:", count_constraint)
    # Addition constraints:
    model.print_information()
    model.export_as_lp(model_path_out)
    print("Wrote model")


def run_hybrid_solver(cqm):
    """Solve CQM using hybrid solver."""

    # Initialize the CQM solver
    sampler = LeapHybridCQMSampler()

    # Solve the problem using the CQM solver
    sampleset = sampler.sample_cqm(cqm, label='Track finding')
    feasible_sampleset = sampleset.filter(lambda row: row.is_feasible)

    try:
        sample = feasible_sampleset.first.sample
        energy = feasible_sampleset.first.energy
        run_time = sampleset.info['run_time'] / 1000000
    except:
        print("\nNo feasible solutions found.")
        exit()

    return sample, energy, run_time


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
        ax.plot(xs=[h1.x, h2.x], ys=[h1.y, h2.y], zs=[h1.z, h2.z], color='blue')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    plt.savefig(out)
    plt.show()


# ------- Main program -------
if __name__ == "__main__":

    hits_path = '../event000001000/volume_id_9/hits-vol_9_20_track.csv'
    hits_volume = read_hits(hits_path)
    hits = dict()
    for k, v in hits_volume.items():
        # print(k, v)
        print("Volume id:", k)
        print("No_layers:", len(v))
        hits = v

    # hits_volume
    layers = list(hits.keys())

    model_path_out = "result_f2_dao/model_docplex_CQM.lp"
    # model_path_out = "result_f2_dao/model_docplex.lp"
    build_model(hits, model_path_out)

    with open(model_path_out, 'rb') as f:
        cqm = dimod.lp.load(f)

    # import time

    # start = time.time()
    sample, energy, run_time = run_hybrid_solver(cqm)
    # end = time.time()
    print("Run time:", run_time)
    print("Objective value:", energy)
    solution_path = "result_f2_dao/solution_dwave.json"
    with open(solution_path, 'w', encoding='utf-8') as f:
        json.dump(sample, f, ensure_ascii=False, indent=4)

    with open(solution_path, 'r', encoding='utf-8') as f:
        result = json.load(f)
    # print(result)
    segments = []
    for var, value in result.items():
        f_p_i_j = var.split('_')
        if value == 0:
            continue
        print(var)
        p = int(f_p_i_j[1])
        i = int(f_p_i_j[2])
        j = int(f_p_i_j[3])
        h_1 = hits[layers[p - 1]][i - 1]
        h_2 = hits[layers[p]][j - 1]
        segments.append([h_1, h_2])
    out = "result_f2_dao/result_dwave.PNG"
    display(hits, segments, out)
