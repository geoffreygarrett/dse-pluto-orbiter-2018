"""
Created by Alejandro Daniel Noel
"""

import json
import os
import numpy as np

from trajectory_tool.monte_carlo_analysis.obj_def import Trajectory

folder = 'earth-jupiter-pluto'

base_path = os.path.dirname(os.path.realpath(__file__))
results_dir = os.path.join(base_path, folder, 'case_data')
boosted_dir = os.path.join(base_path, folder, 'boosted_data')
results_files = [os.path.join(results_dir, file) for file in os.listdir(results_dir) if 'cases' in str(file)]
boosted_files = [os.path.join(boosted_dir, file) for file in os.listdir(boosted_dir) if 'cases' in str(file)]
itinerary = folder.split('-')

cases = []
for results_file in results_files:
    cases += [Trajectory.from_dict(case_dict) for case_dict in json.load(open(results_file, 'r'))]

boosted_cases = []
for boosted_file in boosted_files:
    boosted_cases += [Trajectory.from_dict(case_dict) for case_dict in json.load(open(boosted_file, 'r'))]


def find_best_case(_cases):
    _best_case = None
    _lowest_dv = 1e10
    for _case in _cases:
        if _case.delta_v < _lowest_dv:
            _best_case = _case
            _lowest_dv = _case.delta_v
    return _best_case


print("Best from {} original cases:".format(len(cases)))
best_case = find_best_case(cases)
for key, value in best_case.get_dict().items():
    print(key, value)
print("delta V =", best_case.delta_v)

print("\nBest from {} boosted cases:".format(len(boosted_cases)))
best_case = find_best_case(boosted_cases)
for key, value in best_case.get_dict().items():
    print(key, value)
print("delta V =", best_case.delta_v)

