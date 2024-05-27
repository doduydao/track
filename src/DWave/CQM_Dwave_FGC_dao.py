import matplotlib
import networkx as nx
import dimod
from dimod import ConstrainedQuadraticModel, Binary, quicksum
from dwave.system import LeapHybridCQMSampler
import json

from data import *
from read_data import *

try:
    import matplotlib.pyplot as plt
except ImportError:
    matplotlib.use("agg")
    import matplotlib.pyplot as plt


def run_hybrid_solver(cqm):
    """Solve CQM using hybrid solver."""

    # Initialize the CQM solver
    sampler = LeapHybridCQMSampler()

    # Solve the problem using the CQM solver
    sampleset = sampler.sample_cqm(cqm, label='Track finding FGC')
    feasible_sampleset = sampleset.filter(lambda row: row.is_feasible)
    try:
        sample = feasible_sampleset.first.sample
    except:
        print("\nNo feasible solutions found.")
        exit()
    return sample

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

    hits_path = '/Users/doduydao/daodd/PycharmProjects/Quantum_Research/Tracking/src/data_selected/15hits/unknown_track/hits.csv'
    hits_volume = read_hits(hits_path)
    hits = hits_volume[9]
    # hits_volume

    model_file = "/Users/doduydao/daodd/PycharmProjects/Quantum_Research/Tracking/src/First_General_Case/results/15hits/unknown_track/model_docplex_LB_new_dist_test.lp"

    with open(model_file, 'rb') as f:
        cqm = dimod.lp.load(f)

    sample = run_hybrid_solver(cqm)

    solution_path = "result/solution_dwave.json"
    with open(solution_path, 'w', encoding='utf-8') as f:
        json.dump(sample, f, ensure_ascii=False, indent=4)


    with open(solution_path, 'r', encoding='utf-8') as f:
        result = json.load(f)
    # print(result)
    segments = []
    layers = list(hits.keys())
    layers = [0] + layers + [0]
    for var, value in result.items():
        value = round(value)
        phi_p_p_i_j = var.split('_')
        if 'nt' in phi_p_p_i_j[0] or 'c' in phi_p_p_i_j[0] or value != 1.0 or 's' in phi_p_p_i_j[0] or 'q' in \
                phi_p_p_i_j[0] or 'ob' in \
                phi_p_p_i_j[0]:
            continue

        p_1 = int(phi_p_p_i_j[1])
        p_2 = int(phi_p_p_i_j[2])
        i = int(phi_p_p_i_j[3])
        j = int(phi_p_p_i_j[4])
        h_1 = hits[layers[p_1]][i - 1]
        h_2 = hits[layers[p_2]][j - 1]
        segments.append([h_1, h_2])

    out = "result/result_dwave.PNG"
    display(hits, segments, out)
