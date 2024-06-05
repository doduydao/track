import matplotlib
import dimod
from docplex.mp.model import Model
from dwave.system import LeapHybridCQMSampler
import json

from data import *

try:
    import matplotlib.pyplot as plt
except ImportError:
    matplotlib.use("agg")
    import matplotlib.pyplot as plt


def distance(h1, h2):
    distance = math.sqrt((h2.x - h1.x) ** 2 + (h2.y - h1.y) ** 2 + (h2.z - h1.z) ** 2)
    return distance


def assign_exact(model, f, hits):
    track = dict()
    for p, hp in hits.items():
        for h in hp:
            p_id = h.particle_id / 10000000000
            if p_id not in track:
                track[p_id] = [hp.index(h)]
            else:
                track[p_id] += [hp.index(h)]

    for t, hs in track.items():
        for p in range(len(hs) - 1):
            i = hs[p] + 1
            j = hs[p + 1] + 1
            model.add_constraint(f[p + 1, i, j] == 1)


def build_model(hits, model_path_out):
    model = Model(name="Track")
    layers = list(hits.keys())
    layers = [0] + layers
    L = len(layers)
    no_hits = len(list(hits.values())[1]) + 1
    print(L, no_hits)

    f = model.binary_var_dict(
        [(p, i, j) for p in range(1, L - 1) for i in range(1, no_hits) for j in range(1, no_hits)],
        name="f")

    assign_exact(model, f, hits)

    objective = 0
    for p in range(1, L - 2):
        for i in range(1, no_hits):
            for j in range(1, no_hits):
                for k in range(1, no_hits):
                    h_i = hits[layers[p]][i - 1]
                    h_j = hits[layers[p + 1]][j - 1]
                    h_k = hits[layers[p + 2]][k - 1]
                    seg_1 = Segment(h_j, h_i)
                    seg_2 = Segment(h_j, h_k)
                    angle = Angle(seg_1=seg_1, seg_2=seg_2).angle

                    beta = angle

                    objective += f[p, i, j] * f[p + 1, j, k] * beta
    model.set_objective('min', objective)

    # first constraints:
    print("---First constraints---")
    ct_1 = 0
    for p in range(1, L - 1):
        for j in range(1, no_hits):
            tmp = 0
            for i in range(1, no_hits):
                tmp += f[p, i, j]
            ct_1 += 1
            ctname = "FC_" + str(ct_1)
            model.add_constraint(tmp == 1, ctname=ctname)

    # Second constraints:
    print("---Second constraints---")
    ct_2 = 0
    for p in range(1, L - 1):
        for i in range(1, no_hits):
            tmp = 0
            for j in range(1, no_hits):
                tmp += f[p, i, j]
            ct_2 += 1
            ctname = "SC_" + str(ct_2)
            model.add_constraint(tmp == 1, ctname=ctname)

    model.print_information()
    model.export_as_lp(model_path_out)
    # model.solve(log_output=True)
    # print(model.solution)

    print("Wrote model")


def run_hybrid_solver(cqm):
    """Solve CQM using hybrid solver."""
    print("Solving CQM ...")
    # Initialize the CQM solver
    sampler = LeapHybridCQMSampler()
    # print(sampler.properties["minimum_time_limit_s"])
    # sampler.properties["minimum_time_limit_s"] = 10
    print(sampler.min_time_limit(cqm))
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
    src_path = '/Users/doduydao/daodd/PycharmProjects/track/src/data_selected'
    folder = '/50hits/'

    data_path = src_path + folder + 'known_track/hits.csv'
    model_path_out = "result" + folder + "known_track/CQM_no_dist_no_LB_exact.lp"
    solution_path = "result" + folder + "known_track/solution_no_dist_no_LB_exact.json"
    out = "result" + folder + "known_track/result_dwave_no_dist_no_LB_exact.PNG"
    hits = read_hits(data_path)[9]

    build_model(hits, model_path_out)

    with open(model_path_out, 'rb') as f:
        cqm = dimod.lp.load(f)

    sample, energy, run_time = run_hybrid_solver(cqm)

    print("Run time:", run_time)
    print("Objective value:", energy)

    with open(solution_path, 'w', encoding='utf-8') as f:
        json.dump(sample, f, ensure_ascii=False, indent=4)

    with open(solution_path, 'r', encoding='utf-8') as f:
        result = json.load(f)

    layers = [0] + list(hits.keys())
    segments = []
    for var, value in result.items():
        f_p_i_j = var.split('_')
        if 'f' in f_p_i_j[0] and value == 1.0:
            print(var, value)
            p = int(f_p_i_j[1])
            i = int(f_p_i_j[2])
            j = int(f_p_i_j[3])
            h_1 = hits[layers[p]][i - 1]
            h_2 = hits[layers[p + 1]][j - 1]
            segments.append([h_1, h_2])
    display(hits, segments, out)
