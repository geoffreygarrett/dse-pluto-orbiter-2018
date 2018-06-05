"""
Created by Alejandro Daniel Noel
"""

import datetime
import os
import json

from trajectory_tool.core import TrajectoryTool
from trajectory_tool.monte_carlo_analysis.obj_def import Trajectory

base_path = os.path.dirname(os.path.realpath(__file__))

folder = 'earth-jupiter-pluto'
descr_dir = os.path.join(base_path, folder, 'case_descr')
results_dir = os.path.join(base_path, folder, 'case_data')
case_descr_files = [os.path.join(descr_dir, file) for file in os.listdir(descr_dir) if 'cases' in str(file)]
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


def save_legs(_init_case, _case_num, _trajectories):
    json.dump(_trajectories, open(os.path.join(results_dir, "cases_{}-{}.json".format(_init_case, _case_num)), "w"), indent=4)


# Compute legs
tt = TrajectoryTool()
batch_size = 500
batch_init_case = cases[0]['id']
trajectories = []
for i, case in enumerate(cases):
    print(i)
    # Run case
    result = tt.process_itinerary(case, itinerary, _mode='fast')
    trajectories.append(Trajectory(case=case['id'], itinerary=itinerary,
                                   planet_dates=[p['d'].datetime for p in result],
                                   velocities=[result[0]['v']['d'].value,
                                               *[p['v'][k].value for p in result[1:-1] for k in ['a', 'd']],
                                               result[-1]['v']['a'].value],
                                   planet_velocities=[p['v']['p'].value for p in result]))

    if (i+1) % batch_size == 0:
        # Save legs until now
        save_legs(batch_init_case, case['id'], [t.get_dict() for t in trajectories])
        batch_init_case = case['id'] + 1
        trajectories = []

# final save
if len(trajectories):
    save_legs(batch_init_case, cases[-1]['id'], [t.get_dict() for t in trajectories])
