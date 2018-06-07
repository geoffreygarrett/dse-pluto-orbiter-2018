"""
Created by Alejandro Daniel Noel
"""

import datetime
import json
import os

from trajectory_tool.TestSolutionsMatt import sols2send_1, sols2send_2, sols2send_3, sols2send_4
from trajectory_tool.core import TrajectoryTool
from trajectory_tool.monte_carlo_analysis.obj_def import Trajectory


def save_legs(_init_case, _case_num, _trajectories, results_dir):
    json.dump(_trajectories, open(os.path.join(results_dir, "cases_{}-{}.json".format(_init_case, _case_num)), "w"), indent=4)


def do_runner(itinerary_folder):
    base_path = os.path.dirname(os.path.realpath(__file__))

    folder = itinerary_folder
    descr_dir = os.path.join(base_path, folder, 'case_descr')
    results_dir = os.path.join(base_path, folder, 'case_data')

    # CREATES DIRECTORY IF IT DOESN NOT EXIST
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    case_descr_files = [os.path.join(descr_dir, file) for file in os.listdir(descr_dir) if 'cases' in str(file)]
    print(case_descr_files)
    case_descr_files.sort()
    itinerary = folder.split('-')

    # Load all cases to run form the descr files
    cases = []
    next_case = 0
    for descr_file in case_descr_files.copy():
        with open(descr_file, 'r') as f:
            for line in f.readlines():
                if line[0].isdigit():
                    case_num = eval(line.split(':')[0].strip())
                    if case_num == next_case:
                        cases.append({'id': case_num,
                                      'launch_date': datetime.datetime.strptime(line.split(':')[1].split(',', maxsplit=1)[0].strip(), '%Y-%m-%d'),
                                      'durations': eval(line.split(':')[1].split(',', maxsplit=1)[1].strip())})
                        next_case += 1

    # Compute legs
    tt = TrajectoryTool()
    batch_size = 500
    batch_init_case = cases[0]['id']
    trajectories = []
    for i, case in enumerate(cases):
        # print(i)
        result = None
        try:
            result = tt.process_itinerary(case, itinerary, _mode='delta_v')
            print('worked')
        except ValueError as e:
            continue
        trajectories.append(Trajectory(case=case['id'], itinerary=itinerary,
                                       planet_dates=[p['d'].datetime for p in result],
                                       velocities=[result[0]['v']['d'].value,
                                                   *[p['v'][k].value for p in result[1:-1] for k in ['a', 'd']],
                                                   result[-1]['v']['a'].value],
                                       planet_velocities=[p['v']['p'].value for p in result],
                                       delta_v_dep_arr=[result[0]['dv'], result[-1]['dv']],
                                       delta_v_assists=[p['dv'] for p in result[1:-1]]))

        if (i + 1) % batch_size == 0:
            # Save legs until now
            save_legs(batch_init_case, case['id'], [t.get_dict() for t in trajectories], results_dir)
            batch_init_case = case['id'] + 1
            trajectories = []

    # final save
    if len(trajectories):
        save_legs(batch_init_case, cases[-1]['id'], [t.get_dict() for t in trajectories], results_dir)


# print(sols2send)

def connect(itinerary_list):
    itinerary_string = ''
    for body in itinerary_list:
        if body == 'pluto':
            itinerary_string += body
        else:
            itinerary_string += body + '-'
    return itinerary_string


# print(connect(sols2send_1[0]))

# EDIT THE X FOR A NEW SOLUTION RUN
sols_string = []
for sol in sols2send_3:
    sol_string = connect(sol)
    sols_string.append(sol_string)

print(sols_string)
for sol in sols_string:
    print(sol)
    do_runner(sol)
