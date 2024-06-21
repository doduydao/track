import json
# import pandas as pd
import dimod
def load_dwave_json_result(path):
    with open(path, 'r') as file:
        data = json.load(file)
    # print(data.keys())
    sampleset = dimod.SampleSet.from_serializable(data)
    # print(sampleset.first)

    feasible_sampleset = sampleset.filter(lambda row: row.is_feasible)
    print(feasible_sampleset)
    sample = feasible_sampleset.first.sample
    energy = feasible_sampleset.first.energy
    run_time = sampleset.info['run_time'] / 1000000

    print('energy:',energy)
    print('run_time:', run_time)
    return sample, energy, run_time

if __name__ == '__main__':
    path = '/Users/doduydao/daodd/PycharmProjects/track/src/Bogdan_idea/results/75hits/known_track/answer_01cb8a57-4788-4044-bb26-cc4c153a63e1.json'
    load_dwave_json_result(path)
