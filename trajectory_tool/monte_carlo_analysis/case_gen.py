"""
Created by Alejandro Daniel Noel
"""

import datetime

import numpy as np
import os

from trajectory_tool.ale_solar_system.solar_system import SolarSystem
from trajectory_tool.ale_solar_system.time_utils import datetime_to_jd

from trajectory_tool.TestSolutionsMatt import sols2send

def unit_vector(vector):
    return vector / np.linalg.norm(vector)


def angle_between(body1, body2, epoch, ss, dt=datetime.timedelta(days=0), mind_direction=True):
    r1 = ss.coordinates_of(body1, time=datetime_to_jd(epoch))
    r2 = ss.coordinates_of(body2, time=datetime_to_jd(epoch + dt))
    v1_u = unit_vector(r1)
    v2_u = unit_vector(r2)
    angle = np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
    if mind_direction:
        r1 = ss.coordinates_of(body1, time=datetime_to_jd(epoch + datetime.timedelta(days=1)))
        r2 = ss.coordinates_of(body2, time=datetime_to_jd(epoch + dt + datetime.timedelta(days=1)))
        v1_u = unit_vector(r1)
        v2_u = unit_vector(r2)
        angle2 = np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
        return angle if angle2 < angle else 2 * np.pi - angle
    else:
        return angle


def pseudo_distance(body1, body2, initial_epoch, ss, dt=datetime.timedelta(days=0)):
    r1 = ss.coordinates_of(body1, time=datetime_to_jd(initial_epoch))
    r2 = ss.coordinates_of(body2, time=datetime_to_jd(initial_epoch + dt))
    r1_norm = np.linalg.norm(r1)
    r2_norm = np.linalg.norm(r2)
    rad_dist = abs(r2_norm - r1_norm)
    circ_dist = angle_between(body1, body2, initial_epoch, ss, dt) * (r1_norm + r2_norm) * 0.5
    return rad_dist + circ_dist

def save_cases(year, cases_str, folder, itinerary):
    base_path = os.path.dirname(os.path.realpath(__file__))
    filename = "cases_{}.txt".format(year)
    filepath = os.path.join(base_path, folder, 'case_descr', filename)

    # CREATES DIRECTORY IF IT DOESN NOT EXIST
    directory = os.path.join(base_path, folder, 'case_descr')
    if not os.path.exists(directory):
        os.makedirs(directory)
    print("Saving {}".format(filepath))
    print(filepath)
    with open(filepath, 'w') as f:
        f.write("Itinerary: {}\n".format(str(itinerary)))
        f.write("case #: launch date, tof of each leg\n")
        f.write(cases_str)





def do_case_gen(itinerary):
    launch_window = 2022, 2030
    days_of_month = [1, 10, 20]
    #itinerary = ['earth', 'jupiter', 'pluto']
    #CREATES FOLDER NAME
    folder = ''
    for body in itinerary:
        if body == 'PLUTO_BARYCENTER' or body == 'pluto':
            folder += body
        else:
            folder += body + '-'

    flight_times = np.arange(10.0, 25.0, 0.25)
    #folder = "earth-jupiter-pluto2"


    launch_dates = [datetime.date(year, month, day)
                    for year in range(launch_window[0], launch_window[1])
                    for month in range(1, 12 + 1)
                    for day in days_of_month]

    ss = SolarSystem()


    # def save_cases(year, cases_str):
    #     base_path = os.path.dirname(os.path.realpath(__file__))
    #     filename = "cases_{}.txt".format(year)
    #     filepath = os.path.join(base_path, folder, 'case_descr', filename)
    #
    #     #CREATES DIRECTORY IF IT DOESN NOT EXIST
    #     directory = os.path.join(base_path, folder, 'case_descr')
    #     if not os.path.exists(directory):
    #         os.makedirs(directory)
    #     print("Saving {}".format(filepath))
    #     print(filepath)
    #     with open(filepath, 'w') as f:
    #         f.write("Itinerary: {}\n".format(str(itinerary)))
    #         f.write("case #: launch date, tof of each leg\n")
    #         f.write(cases_str)


    cases = ""

    # Check for previous partial case generations
    case_num = 0
    # descr_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), folder, 'case_descr')
    # case_descr_files = [os.path.join(descr_dir, file) for file in os.listdir(descr_dir) if 'cases' in str(file)]
    # for descr_file in case_descr_files.copy():
    #     with open(descr_file, 'r') as f:
    #         for line in f.readlines():
    #                 case_num = max(case_num, eval(line.split(':')[0].strip()))

    for launch_date in launch_dates:
        if launch_date.month == 1 and launch_date.day == 1:
            save_cases(launch_date.year, cases, folder, itinerary)
            cases = ""
        for tof in flight_times:
            # Compute pseudo travel times between itinerary way-points
            pseudo_times = np.zeros(len(itinerary)-1)
            for _ in range(8):  # In case pseudo_distance() depends on planets positions, iterate a bit for some convergence
                pseudo_distances = [0] * (len(itinerary) - 1)
                for i, (first, second) in enumerate(zip(itinerary, itinerary[1:])):
                    pseudo_distances[i] = pseudo_distance(first, second,
                                                          launch_date + datetime.timedelta(days=365*sum(pseudo_times[:i])), ss,
                                                          datetime.timedelta(days=365 * pseudo_times[i]))
                pseudo_times = np.array(pseudo_distances) * tof / sum(pseudo_distances)
            cases += "{}: {}, {}\n".format(str(case_num).ljust(6), launch_date.strftime("%Y-%m-%d"), str(list(pseudo_times)))
            case_num += 1
        print(launch_date)
        if launch_date.month == 12 and launch_date.day == days_of_month[-1]:
            save_cases(launch_date.year, cases, folder, itinerary)


for itinerary in sols2send:
    do_case_gen(itinerary)