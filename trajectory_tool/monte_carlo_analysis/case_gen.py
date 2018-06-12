"""
Created by Alejandro Daniel Noel
"""

import datetime

import numpy as np
import os
import copy

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


def fake_flight_time(body1, body2, initial_epoch, ss, dt=0):
    r1 = ss.coordinates_of(body1, time=datetime_to_jd(initial_epoch))
    r2 = ss.coordinates_of(body2, time=datetime_to_jd(initial_epoch + datetime.timedelta(days=365.25 * dt)))
    r1_norm = np.linalg.norm(r1)
    r2_norm = np.linalg.norm(r2)
    rad_dist = abs(r2_norm - r1_norm)
    r_ave = (r1_norm + r2_norm) * 0.5
    v = np.sqrt(1.327124e20 / (r1_norm * 1000.0)) / 1000.0 + 6
    for _ in range(8):
        circ_dist = angle_between(body1, body2, initial_epoch, ss, datetime.timedelta(days=365.25 * dt)) * r_ave
        dt = (rad_dist + circ_dist) / v / 3600 / 24
    return dt


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


inner_planets_periods = {'mercury': 0.2411, 'venus': 0.69863, 'earth': 1.0, 'mars': 1.8822}


def do_case_gen(itinerary):
    launch_window = 2022, 2030
    days_of_month = [1, 10, 20]
    steps_in_inner_planets_rev = 10

    folder = '-'.join(itinerary)

    launch_dates = [datetime.date(year, month, day)
                    for year in range(launch_window[0], launch_window[1])
                    for month in range(1, 12 + 1)
                    for day in days_of_month]

    ss = SolarSystem()

    cases = ""
    case_num = 0

    flight_times = [np.zeros(len(itinerary) - 1)]
    for launch_date in launch_dates:
        for i, (planet1, planet2) in enumerate(zip(itinerary[:-1], itinerary[1:])):
            if planet2 in inner_planets_periods.keys():
                for time_delta in np.linspace(0.0, inner_planets_periods[planet2], steps_in_inner_planets_rev):
                    flight_times.append(copy.deepcopy(flight_times[-1]))
                    flight_times[-1][i:] = 0
                    departure_date = launch_date + datetime.timedelta(days=365.25 * sum(flight_times[-1][:i]) + time_delta)
                    flight_times[-1][i] = fake_flight_time(planet1, planet2, departure_date, ss)
            else:
                flight_times.append(copy.deepcopy(flight_times[-1]))
                flight_times[-1][i:] = 0
                departure_date = launch_date + datetime.timedelta(days=365.25 * sum(flight_times[-1][:i]))
                flight_times[-1][i] = fake_flight_time(planet1, planet2, departure_date, ss)

        for flight_time in flight_times:
            if sum(flight_time) < 1.0:
                continue
            cases += "{}: {}, {}\n".format(str(case_num).ljust(6), launch_date.strftime("%Y-%m-%d"), str(list(flight_time)))
            case_num += 1
            flight_times = [np.zeros(len(itinerary) - 1)]

        if launch_date.month == 12 and launch_date.day == days_of_month[-1]:
            save_cases(launch_date.year, cases, folder, itinerary)
            cases = ""


do_case_gen(['earth', 'jupiter', 'pluto'])
# for itinerary in sols2send:
#     do_case_gen(itinerary)