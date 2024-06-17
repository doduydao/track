import unittest
import dwave_peterson_cost_qubo as peterson_qubo
import math

src_path = '../../src/data_selected'
folder = '/2hits/known_track/'
out_path = '/Users/doduydao/daodd/PycharmProjects/track/src/Bogdan_idea/results'
peterson_qubo.check_path(out_path + folder)
data_path = src_path + folder + 'hits.csv'
costs_path_out = out_path + folder + "costs.json"
hits_by_layers = peterson_qubo.read_hits(data_path)[9]
list_hits = []

for hs in list(hits_by_layers.values()):
    list_hits += hs

beta_max = math.pi
m = 7
alpha = 1
beta = 1
costs = peterson_qubo.get_costs(list_hits, hits_by_layers, beta_max)
peterson_qubo.write_costs(costs, costs_path_out, m)
costs = peterson_qubo.load_costs(costs_path_out)


class PetersonTest(unittest.TestCase):

    def test_define_variables(self):
        x_1 = peterson_qubo.define_variables(costs)

        gold = []

        s_j, e_j = 2, 7
        gold += [(i, j) for i in [0, 1] for j in range(s_j, e_j + 1)]

        s_j, e_j = 4, 9
        gold += [(i, j) for i in [2, 3] for j in range(s_j, e_j + 1)]

        s_j, e_j = 6, 11
        gold += [(i, j) for i in [4, 5] for j in range(s_j, e_j + 1)]

        s_j, e_j = 8, 13
        gold += [(i, j) for i in [6, 7] for j in range(s_j, e_j + 1)]

        s_j, e_j = 10, 13
        gold += [(i, j) for i in [8, 9] for j in range(s_j, e_j + 1)]

        s_j, e_j = 12, 13
        gold += [(i, j) for i in [10, 11] for j in range(s_j, e_j + 1)]

        self.assertEqual(len(x_1.keys()), len(gold))
        self.assertEqual(sorted(list(x_1.keys())), sorted(gold))  # add assertion here

    def test_create_objective_function(self):
        H = peterson_qubo.create_objective_function(list_hits, costs, m, alpha, beta)
        gold = 0

        x = peterson_qubo.define_variables(costs)
        hit_last_layer = hits_by_layers[m]
        N = len(list_hits) - len(hit_last_layer)
        print("N =", N)

        segments = set()
        first_part = 0
        for id, cost in costs.items():
            i_j = id[0]
            j_k = id[1]
            first_part += cost * x[i_j] * x[j_k]
            segments.add(i_j)
            segments.add(j_k)

        sum_segments = sum([x[s] for s in segments])
        second_part = (sum_segments - N) ** 2

        part_3 = ((x[(0, 2)] + x[(1, 2)] - 1) ** 2 +
                  (x[(0, 3)] + x[(1, 3)] - 1) ** 2 +
                  (x[(0, 3)] + x[(1, 3)] - 1) ** 2 +
                  (x[(0, 4)] + x[(1, 4)] + x[(2, 4)] + x[(3, 4)] - 1) ** 2 +
                  (x[(0, 5)] + x[(1, 5)] + x[(2, 5)] + x[(3, 5)] - 1) ** 2 +
                  (x[(0, 6)] + x[(1, 6)] + x[(2, 6)] + x[(3, 6)] + x[(4, 6)] + x[(5, 6)] - 1) ** 2 +
                  (x[(0, 7)] + x[(1, 7)] + x[(2, 7)] + x[(3, 7)] + x[(4, 7)] + x[(5, 7)] - 1) ** 2


                  )

        self.assertEqual(H, gold)  # add assertion here


if __name__ == '__main__':
    unittest.main()

    # ptest = PetersonTest(costs)
    # ptest.test_define_variables()
