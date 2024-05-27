from data import *
from read_data import *
# import pulp
import cplex
from docplex.mp.model import Model
import json
from show_3D import *
import random


def run(hits, model_path_out, solution_path_out):
    model = Model(name="Track")

    layers = list(hits.keys())
    print(layers)
    no_layer = len(layers)
    no_hits = len(list(hits.values())[0])
    no_track = no_hits
    print(no_layer, no_hits)

    x = model.binary_var_dict(
        [(p, t, i) for t in range(1, no_track + 1) for p in range(1, no_layer + 1) for i in range(1, no_hits + 1)],
        name="x")

    # y = model.binary_var_dict(
    #     [(p, t, i, j, k) for t in range(1, no_track + 1) for p in range(1, no_layer - 1) for i in range(1, no_hits + 1)
    #      for j in range(1, no_hits + 1) for k in range(1, no_hits + 1)], name="y")

    y = model.continuous_var_dict(
        [(p, t, i, j, k) for t in range(1, no_track + 1) for p in range(1, no_layer - 1) for i in range(1, no_hits + 1)
         for j in range(1, no_hits + 1) for k in range(1, no_hits + 1)], name="y", ub=1, lb=0)

    objective = 0
    for t in range(1, no_track + 1):
        for i in range(1, no_hits + 1):
            for j in range(1, no_hits + 1):
                for k in range(1, no_hits + 1):
                    for p in range(1, no_layer - 1):
                        h_i = hits[layers[p - 1]][i - 1]
                        h_j = hits[layers[p]][j - 1]
                        h_k = hits[layers[p + 1]][k - 1]
                        seg_1 = Segment(h_j, h_i)
                        seg_2 = Segment(h_j, h_k)
                        # print(Angle(seg_1=seg_1, seg_2=seg_2).angle)
                        objective += y[p, t, i, j, k] * Angle(seg_1=seg_1, seg_2=seg_2).angle
                    # break

    model.set_objective('min', objective)
    # Constraints
    # first constraints:
    print("---First constraints---")
    count_constraint = 0
    for i in range(1, no_hits + 1):
        for p in range(1, no_layer + 1):
            tmp = 0
            for t in range(1, no_track + 1):
                tmp += x[p, t, i]
            constraint_name = "FC_" + str(count_constraint)
            count_constraint += 1
            model.add_constraint(tmp == 1, ctname=constraint_name)
    print("Number of first constraints:", count_constraint)
    # print()
    # Second constraints:
    print("---Second constraints---")
    count_constraint = 0
    for t in range(1, no_track + 1):
        for p in range(1, no_layer + 1):
            tmp = 0
            for i in range(1, no_hits + 1):
                tmp += x[p, t, i]
            constraint_name = "SC_" + str(count_constraint)
            count_constraint += 1
            model.add_constraint(tmp == 1, ctname=constraint_name)
    print("Number of second constraints:", count_constraint)
    # Addition constraints:
    print("---Addition constraints---")
    count_constraint = 0
    for t in range(1, no_track + 1):
        for p in range(1, no_layer - 1):
            for i in range(1, no_hits + 1):
                for j in range(1, no_hits + 1):
                    for k in range(1, no_hits + 1):
                        c1 = x[p, t, i] + x[p + 1, t, j] + x[p + 2, t, k] - y[p, t, i, j, k] <= 2
                        c2 = y[p, t, i, j, k] <= x[p, t, i]
                        c3 = y[p, t, i, j, k] <= x[p + 1, t, j]
                        c4 = y[p, t, i, j, k] <= x[p + 2, t, k]
                        c5 = y[p, t, i, j, k] >= 0
                        constraint_1_name = "AC_" + str(count_constraint) + "_1"
                        constraint_2_name = "AC_" + str(count_constraint) + "_2"
                        constraint_3_name = "AC_" + str(count_constraint) + "_3"
                        constraint_4_name = "AC_" + str(count_constraint) + "_4"
                        constraint_5_name = "AC_" + str(count_constraint) + "_5"
                        count_constraint += 5

                        model.add_constraint(c1, ctname=constraint_1_name)
                        model.add_constraint(c2, ctname=constraint_2_name)
                        model.add_constraint(c3, ctname=constraint_3_name)
                        model.add_constraint(c4, ctname=constraint_4_name)
                        model.add_constraint(c5, ctname=constraint_5_name)
    print("Number of addition constraints:", count_constraint)

    model.print_information()
    model.solve(log_output=True)

    model.export_as_lp(model_path_out)
    model.solution.export(solution_path_out)
    f = open(solution_path_out)
    result = json.load(f)
    f.close()

    return result['CPLEXSolution']['variables']


def pick_random_hits(no_track, hits):
    layers = list(hits.keys())
    new_hits = dict()
    for p, hp in hits.items:
        idx = []
        while len(idx) < no_track:
            id = random.randint(0, len(hp))
            if id not in idx:
                idx.append(id)
        new_hp = [hp[i] for i in idx]
        new_hits[p] = new_hp
    return new_hits


if __name__ == '__main__':
    # hits_path = 'C:\\Users\dddo\PycharmProjects\Quantum_Research\Tracking\event000001000\event000001000-hits.csv'
    # hits_path = 'C:\\Users\dddo\PycharmProjects\Quantum_Research\Tracking\event000001000\sel\event000001000-hits-sel-01.csv'
    hits_path = '/Users/doduydao/daodd/PycharmProjects/Quantum_Research/Tracking/event000001000/sel/event000001000-hits-sel-01.csv'
    # hits_path = 'C:\\Users\dddo\PycharmProjects\Quantum_Research\Tracking\event000001000\sublayer_2\event000001000-hits_random.csv'
    hits_volume = read_hits(hits_path)
    hits = dict()
    for k, v in hits_volume.items():
        # print(k, v)
        print("Volume id:", k)
        print("No_layers:", len(v))
        hits = v

    # hits_volume
    layers = list(hits.keys())

    model_path_out = "/Users/doduydao/daodd/PycharmProjects/Quantum_Research/Tracking/src/result_f1/model_docplex.lp"
    solution_path_out = "/Users/doduydao/daodd/PycharmProjects/Quantum_Research/Tracking/src/result_f1/solution.json"
    result = run(hits, model_path_out, solution_path_out)
    solution = dict()
    for var in result:
        print(var)
        x_p_t_i = var['name'].split('_')
        if x_p_t_i[0] == 'y':
            continue
        p = int(x_p_t_i[1])
        t = int(x_p_t_i[2])
        i = int(x_p_t_i[3])
        # print(var['name'])

        if t not in solution:
            solution[t] = [hits[layers[p - 1]][i - 1]]
        else:
            solution[t] += [hits[layers[p - 1]][i - 1]]

    out = "/Users/doduydao/daodd/PycharmProjects/Quantum_Research/Tracking/src/result_f1/result.PNG"
    display(hits, solution, out)
