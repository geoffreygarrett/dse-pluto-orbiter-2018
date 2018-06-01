"""
Created by Alejandro Daniel Noel
"""

import os
import datetime
base_path = os.path.dirname(os.path.realpath(__file__))


folder = 'earth-venus-jupiter-pluto'
descr_dir = os.path.join(base_path, folder, 'case_descr')
case_descr_files = [os.path.join(descr_dir, file) for file in os.listdir(descr_dir) if 'cases' in str(file)]
case_descr_files.sort()
itinerary = eval(open(case_descr_files[0], 'r').readline().replace('Itinerary: ', '').strip())
next_case = 1 + eval(open(os.path.join(descr_dir, 'log.txt'), 'r').readline().replace('last case done: ', '').strip())


cases = []
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

a = 0
for case in cases:
    # Run case
    ...
    open(os.path.join(descr_dir, 'log.txt'), 'w').write("")




