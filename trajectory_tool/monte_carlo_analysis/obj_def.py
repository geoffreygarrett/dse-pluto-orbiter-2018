"""
Created by Alejandro Daniel Noel
"""

import datetime
import numpy as np


class Trajectory:
    def __init__(self, case, itinerary, planet_dates, velocities, planet_velocities):
        self.case = case
        self.itinerary = itinerary
        self.planet_dates = planet_dates
        self.velocities = velocities
        self.planet_velocities = planet_velocities
        self.delta_v_dep_arr = []
        self.delta_v_assists = []

    def get_dict(self):
        return {'case': self.case,
                'itinerary': self.itinerary,
                'planet_dates': [date.strftime('%Y-%m-%d') for date in self.planet_dates],
                'velocities': [list(v) for v in self.velocities],
                'planet_velocities': [list(pv) for pv in self.planet_velocities]}

    @classmethod
    def from_dict(cls, dictionary):
        return Trajectory(case=dictionary['case'],
                          itinerary=dictionary['itinerary'],
                          planet_dates=[datetime.datetime.strptime(d, '%Y-%m-%d') for d in dictionary['planet_dates']],
                          velocities=[np.array(v) for v in dictionary['velocities']],
                          planet_velocities=[np.array(v) for v in dictionary['planet_velocities']])

    @property
    def delta_v(self):
        if len(self.delta_v_dep_arr) == 0 or len(self.delta_v_assists) == 0:
            return None
        else:
            return np.linalg.norm(sum(self.delta_v_dep_arr) + sum(self.delta_v_assists))

