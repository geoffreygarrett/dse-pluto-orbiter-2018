"""
Created by Alejandro Daniel Noel

The purpose of this script is to find from the results if other combinations
of launch dates and flight times between assists could improve the trajectory.
It does so by finding nearby dates for the computed gravity assist events and
checking of a gravity assist is feasible for each pair of incoming and outgoing
velocities.
"""

import json
import os

from trajectory_tool.monte_carlo_analysis.obj_def import Trajectory

folder = 'earth-jupiter-pluto'

base_path = os.path.dirname(os.path.realpath(__file__))
results_dir = os.path.join(base_path, folder, 'case_data')
boosted_dir = os.path.join(base_path, folder, 'boosted_data')

# CREATES DIRECTORY IF IT DOESN NOT EXIST
directory = boosted_dir
if not os.path.exists(directory):
    os.makedirs(directory)

results_files = [os.path.join(results_dir, file) for file in os.listdir(results_dir) if 'cases' in str(file)]
itinerary = folder.split('-')

cases = []
for results_file in results_files:
    cases += [Trajectory.from_dict(case_dict) for case_dict in json.load(open(results_file, 'r'))]

days_for_close_check = 2

boosted_cases = {}
labels_next_save = []

index_bar = 0
count_done = 0
for idx, assist_planet in zip(range(1, len(itinerary) - 1), itinerary[1:-1]):
    for case in cases:
        count_done += 1
        if len(labels_next_save) >= 10000:
            print("Saving 10000 boosted cases ({}% done)".format(100 * float(count_done) / ((len(itinerary) - 2) * len(cases))))
            json.dump([boosted_cases[k].get_dict() for k in labels_next_save],
                      open(os.path.join(boosted_dir, "boosted_cases_{}.json".format(index_bar)), 'w'),
                      indent=4)
            labels_next_save = []
            index_bar += 1
        assist_date = case.planet_dates[idx]
        for temp_case in cases:
            if temp_case.case == case.case:
                continue
            temp_date = temp_case.planet_dates[idx]
            if abs((assist_date - temp_date).days) <= days_for_close_check:
                key = str(int(case.case)) + '-' + str(int(temp_case.case)) + assist_planet
                _key = str(int(temp_case.case)) + '-' + str(int(case.case)) + assist_planet
                if key in boosted_cases.keys() or _key in boosted_cases.keys():
                    continue
                boosted_cases[key] = Trajectory(case=len(boosted_cases),
                                                itinerary=itinerary,
                                                planet_dates=case.planet_dates[:idx + 1] + temp_case.planet_dates[idx + 1:],
                                                velocities=case.velocities[:idx + 1] + temp_case.velocities[idx + 1:],
                                                planet_velocities=case.planet_velocities[:idx + 1] + temp_case.planet_velocities[idx + 1:],
                                                delta_v_dep_arr=[case.delta_v_dep_arr[0], temp_case.delta_v_dep_arr[-1]],
                                                delta_v_assists=case.delta_v_assists[:idx + 1] + temp_case.delta_v_assists[idx + 1:])
                labels_next_save.append(key)
                # print("({}) Joining at {} cases {} and {} (dates {} and {})".format(len(boosted_cases),
                #                                                                     assist_planet,
                #                                                                     case.case,
                #                                                                     temp_case.case,
                #                                                                     case.planet_dates[idx].strftime('%Y-%m-%d'),
                #                                                                     temp_case.planet_dates[idx].strftime('%Y-%m-%d')))

if len(labels_next_save):
    json.dump([boosted_cases[k].get_dict() for k in labels_next_save],
              open(os.path.join(boosted_dir, "boosted_cases_{}.json".format(index_bar)), 'w'),
              indent=4)
