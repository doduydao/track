from data import *
from read_data import *
from docplex.mp.model import Model
import json
import matplotlib.pyplot as plt


def distance(h1, h2):
    distance = math.sqrt((h2.x - h1.x) ** 2 + (h2.y - h1.y) ** 2 + (h2.z - h1.z) ** 2)
    return distance


def run(hits, model_path_out, solution_path_out, figure_path_out):
    model = Model(name="Track")
    layers = list(hits.keys())
    print("layers:", layers)
    no_layer = len(layers)
    no_hits = len(list(hits.values())[0])
    print(no_layer, no_hits)

    f = model.binary_var_dict(
        [(p, i, j) for p in range(1, no_layer) for i in range(1, no_hits + 1) for j in range(1, no_hits + 1)],
        name="f")

    z = model.continuous_var_dict(
        [(p, i, j, k) for p in range(1, no_layer - 1) for i in range(1, no_hits + 1) for j in range(1, no_hits + 1)
         for
         k in
         range(1, no_hits + 1)], name="z", lb=0, ub=1)

    objective = 0
    LB = 0

    for p in range(1, no_layer - 1):
        beta_lower = 100000
        for i in range(1, no_hits + 1):
            for j in range(1, no_hits + 1):
                for k in range(1, no_hits + 1):
                    h_i = hits[layers[p - 1]][i - 1]
                    h_j = hits[layers[p]][j - 1]
                    h_k = hits[layers[p + 1]][k - 1]
                    seg_1 = Segment(h_j, h_i)
                    seg_2 = Segment(h_j, h_k)
                    angle = Angle(seg_1=seg_1, seg_2=seg_2).angle
                    dist = 1
                    beta = angle * dist
                    if beta < beta_lower:
                        beta_lower = beta
                    objective += z[p, i, j, k] * beta
        LB += beta_lower
    model.set_objective('min', objective)

    model.add_constraint(objective >= LB * no_hits, ctname="LB of objective value")

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
    print("---Addition constraints---")
    count_constraint = 0
    for p in range(1, no_layer - 1):
        for i in range(1, no_hits + 1):
            for j in range(1, no_hits + 1):
                for k in range(1, no_hits + 1):
                    c1 = f[p, i, j] + f[p + 1, j, k] - z[p, i, j, k] <= 1
                    c2 = z[p, i, j, k] <= f[p, i, j]
                    c3 = z[p, i, j, k] <= f[p + 1, j, k]
                    c4 = z[p, i, j, k] >= 0
                    constraint_1_name = "AC_" + str(count_constraint) + "_1"
                    constraint_2_name = "AC_" + str(count_constraint) + "_2"
                    constraint_3_name = "AC_" + str(count_constraint) + "_3"
                    constraint_4_name = "AC_" + str(count_constraint) + "_4"
                    count_constraint += 4

                    model.add_constraint(c1, ctname=constraint_1_name)
                    model.add_constraint(c2, ctname=constraint_2_name)
                    model.add_constraint(c3, ctname=constraint_3_name)
                    model.add_constraint(c4, ctname=constraint_4_name)
    print("Number of addition constraints:", count_constraint)

    model.print_information()
    model.solve(log_output=True)

    model.export_as_lp(model_path_out)
    model.solution.export(solution_path_out)
    f = open(solution_path_out)
    result = json.load(f)
    f.close()

    result = result['CPLEXSolution']['variables']

    segments = []

    for var in result:
        print(var)
        f_p_i_j = var['name'].split('_')
        if 'f' in f_p_i_j[0]:

            p = int(f_p_i_j[1])
            i = int(f_p_i_j[2])
            j = int(f_p_i_j[3])

            h_1 = hits[layers[p - 1]][i - 1]
            h_2 = hits[layers[p]][j - 1]
            segments.append([h_1, h_2])

    display(hits, segments, figure_path_out)


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


if __name__ == '__main__':
    src_path = '/Users/doduydao/daodd/PycharmProjects/Quantum_Research/Tracking/src/data_selected'
    data_path = src_path + '/6hits/known_track/hits.csv'
    hits = read_hits(data_path)[9]

    model_path_out = "results/6hits/known_track/model_docplex_LB_no_dist.lp"
    solution_path_out = "results/6hits/known_track/solution_LB_no_dist.json"
    figure_path_out = "results/6hits/known_track/result_LB_no_dist.PNG"

    result = run(hits, model_path_out, solution_path_out, figure_path_out)
