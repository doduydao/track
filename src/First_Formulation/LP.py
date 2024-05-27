from data import *
from read_data import *
import pulp
import cplex
import json
from show_3D import *
def run(hits, out_result):
    # define model
    model = pulp.LpProblem('Track_Finding', pulp.LpMinimize)

    # define variable
    layers = list(hits.keys())
    print(layers)
    no_layer = len(layers)
    no_hits = len(list(hits.values())[0])
    no_track = no_hits
    print(no_layer, no_hits)

    x = dict()
    for t in range(1, no_track + 1):
        for p in range(1, no_layer + 1):
            for i in range(1, no_hits + 1):
                x_p_t_i = pulp.LpVariable(name=f'x_{p}_{t}_{i}', cat='Binary')
                key = str(p) + str(t) + str(i)
                x[key] = x_p_t_i
    # print(x)
    y = dict()
    for t in range(1, no_track + 1):
        for p in range(1, no_layer + 1):
            for i in range(1, no_hits + 1):
                for j in range(1, no_hits + 1):
                    for k in range(1, no_hits + 1):
                        y_p_t_i_j_k = pulp.LpVariable(name=f'y_{p}_{t}_{i}_{j}_{k}', cat='Binary')
                        key = str(p) + str(t) + str(i) + str(j) + str(k)
                        y[key] = y_p_t_i_j_k
    # print(y)
    # Objective function
    objective = None
    for t in range(1, no_track + 1):
        for i in range(1, no_hits + 1):
            for j in range(1, no_hits + 1):
                for k in range(1, no_hits + 1):
                    for p in range(0, no_layer - 2):
                        h_i = hits[layers[p]][i - 1]
                        h_j = hits[layers[p + 1]][j - 1]
                        h_k = hits[layers[p + 2]][k - 1]
                        seg_1 = Segment(h_i, h_j)
                        seg_2 = Segment(h_j, h_k)
                        key = str(p+1) + str(t) + str(i) + str(j) + str(k)
                        objective += y[key] * Angle(seg_1=seg_1, seg_2=seg_2).angle
    model.setObjective(objective)

    # Constraints

    # first constraints:
    # print("First constraints:")
    count_constraint = 1
    for i in range(1, no_hits + 1):
        for p in range(1, no_layer + 1):
            tmp = 0
            for t in range(1, no_track + 1):
                key = str(p) + str(t) + str(i)
                tmp += x[key]
            constraint = tmp == 1
            # print(constraint)
            constraint_name = "FC_" + str(count_constraint)
            count_constraint += 1
            model.addConstraint(constraint, name=constraint_name)

    # print()
    # Second constraints:
    # print("Second constraints:")
    count_constraint = 1
    for t in range(1, no_track + 1):
        for p in range(1, no_layer + 1):
            tmp = 0
            for i in range(1, no_hits + 1):
                key = str(p) + str(t) + str(i)
                tmp += x[key]
            constraint = tmp == 1
            # print(constraint)

            constraint_name = "SC_" + str(count_constraint)
            count_constraint += 1
            model.addConstraint(constraint, name=constraint_name)
    # print()

    # Addition constraints:
    # print(y_p_t_i_j_k)

    print("Addition constraints")
    count_constraint = 1
    for t in range(1, no_track + 1):
        for p in range(1, no_layer-1):
            for i in range(1, no_hits + 1):
                for j in range(1, no_hits + 1):
                    for k in range(1, no_hits + 1):
                        key_1 = str(p) + str(t) + str(i)
                        key_2 = str(p+1) + str(t) + str(j)
                        key_3 = str(p+2) + str(t) + str(k)
                        key_y = str(p) + str(t) + str(i) + str(j) + str(k)

                        constraint_1 = y[key_y] <= x[key_1]
                        constraint_2 = y[key_y] <= x[key_2]
                        constraint_3 = y[key_y] <= x[key_3]
                        constraint_4 = x[key_1] + x[key_2] + x[key_3] - y[key_y] <= 2

                        constraint_1_name = "AC_"+ str(count_constraint)+"_1"
                        constraint_2_name = "AC_" + str(count_constraint)+"_2"
                        constraint_3_name = "AC_" + str(count_constraint)+"_3"
                        constraint_4_name = "AC_" + str(count_constraint)+"_4"
                        count_constraint += 1
                        model.addConstraint(constraint_1, name=constraint_1_name)
                        model.addConstraint(constraint_2, name=constraint_2_name)
                        model.addConstraint(constraint_3, name=constraint_3_name)
                        model.addConstraint(constraint_4, name=constraint_4_name)



    print(model, file=open('model_dao.txt', 'w'))
    # Solving
    # solver = pulp.getSolver('CPLEX_CMD')
    solver = pulp.CPLEX_PY()
    model.solve(solver)
    print(model.sol_status)

    # model.writeLP("just_to_be_sure.txt")

    vars = dict()

    for var in model.variables():
        vars[var.name] = var.varValue
        print(var.name, var.varValue)

    with open(out_result, "w") as outfile:
        json.dump(vars, outfile)

    return vars

if __name__ == '__main__':
    # hits_path = '/Users/doduydao/daodd/PycharmProjects/Quantum_Research/Tracking/event000001000/event000001000-hits.csv'
    # hits_path = '/Users/doduydao/daodd/PycharmProjects/Quantum_Research/Tracking/event000001000/sel/event000001000-hits-sel-01.csv'
    hits_path = '/Users/doduydao/daodd/PycharmProjects/Quantum_Research/Tracking/event000001000/sel/event000001000-hits-sel-vol_9_10_track.csv'

    hits_volume = read_hits(hits_path)
    hits = dict()
    for k, v in hits_volume.items():
        print("Volume id:", k)
        print("No_layers:", len(v))
        hits = v

    layers = list(hits.keys())

    # no_track = 2
    # no_layer = 5
    # for l in layers[:no_layer]:
    #     hs = hits[l]
    #     hits_test[l] = hs
        # for h in hs:
        #     print(l, h.id)
    out_result = 'result_pulp_dao.json'
    result = run(hits, out_result)

    # f = open(out_result)
    # result = json.load(f)
    # f.close()

    solution = dict()
    for var_name, var_value in result.items():
        x_p_t_i = var_name.split('_')
        if x_p_t_i[0] == 'y' or var_value == 0:
            continue
        p = int(x_p_t_i[1])
        t = int(x_p_t_i[2])
        i = int(x_p_t_i[3])
        print(var_name, var_value)
        if t not in solution:
            solution[t] = [hits[layers[p - 1]][i - 1]]
        else:
            solution[t] += [hits[layers[p - 1]][i - 1]]

    out = "C:\\Users\dddo\PycharmProjects\Quantum_Research\Tracking\event000001000\sublayer_2\\result_Bogdan.PNG"
    # out = "/Users/doduydao/daodd/PycharmProjects/Quantum_Research/Tracking/event000001000/sublayer_2/result_bogdan_2.PNG"
    display(hits, solution, out)
