from random import randint, choice, uniform
from copy import deepcopy
from datetime import datetime
from trajectory_tool.core import *

import scipy.stats as ss
import numpy as np

LAST_LEG_DURATION_BIAS = 1.5
START_EPOCH = datetime.datetime(2025, 1, 1, 0, 0, 0, 0)

planet_mapping = {'1': Mercury,
                  '2': Venus,
                  '3': Earth,
                  '4': Mars,
                  '5': Jupiter,
                  '6': Saturn,
                  '7': Uranus,
                  '8': Neptune,
                  '9': Pluto}


# Here is the definition of the fitness function.
def fitness_function(chromosone_singleton, chromosone, tt):
    assert chromosone_singleton.unary_schema == list('123456789')
    assert chromosone_singleton.total_schema == '0000 00000 00000 0000'

    _raw_itinerary, _body_list = chromosone_singleton.mapper(chromosone)
    try:
        results = tt.process_itinerary(_raw_itinerary, _body_list, _mode='fast')
        delta_v = sum(results[i]['dv'] for i in range(len(results)))
    except ValueError:
        delta_v = 100
    return 15.0 - delta_v


class Chromosone(object):

    @staticmethod
    def mapper(chromosone):
        sequences = chromosone.split(' ')

        assists = list(filter(lambda a: a != '00000', sequences[1:-1]))
        print(assists)

        _id = chromosone
        _launch_day = START_EPOCH + datetime.timedelta(days=int(sequences[0]))
        _durations = []
        _body_list = []

        for assist in assists:
            planet = planet_mapping[assist[0]].__str__().split(' ')[0].lower()
            duration = float(assist[1:]) / 365
            _durations.append(duration)
            _body_list.append(planet)

        _durations = _durations + [float(sequences[-1]) / 365]

        _raw_itinerary = {'id': _id,
                          'launch_date': _launch_day,
                          'durations': _durations}

        _body_list = ['earth'] + _body_list + ['pluto']

        return _raw_itinerary, _body_list

    def __init__(self, _unary_schema, _total_schema):
        self._unary_schema = _unary_schema
        self._total_schema = _total_schema
        self._tt = TrajectoryTool()

    def fitness(self, _chromosone, _fitness_function, limit=None):
        return fitness_function(self, _chromosone, self._tt)

    @property
    def unary_schema(self):
        return self._unary_schema

    @property
    def total_schema(self):
        return self._total_schema

    def crossover(self, chromosone1, chromosone2):
        pass

    # @property
    def random_chromosone(self, duration):
        duration_days = duration * 365

        chromosone = deepcopy(self._total_schema)
        # (0000) 00000 00000 0000
        chromosone = chromosone.replace('0000', str(randint(1, 9999)).zfill(4), 1)

        # 0000 (0)0000 (0)0000 0000
        _rng_qty_assists = randint(1, len(self._total_schema.split(' ')) - 2)
        _rng_planets = []

        for i in range(_rng_qty_assists):
            _temp = deepcopy(self._unary_schema)
            _temp = list(set(_temp) - set(['9']))  # Pluto removed as gravity assist.
            if len(_rng_planets) > 0:
                _temp = list(set(_temp) - set(_rng_planets[-1]))  # No multi-revolutionary gravity assists.
            if i is 0:
                _temp = list(set(_temp) - set(['3']))  # Remove Earth as possibility for first assist.
            _rng_planets.append(choice(_temp))
        sequence_planets = [' {}0000'.format(planet) for planet in _rng_planets[:]]

        # 0000 0(0000) 0(0000) 0000
        _rng_durations_of_legs_temp = [uniform(1.0, 100.0) for _ in range(_rng_qty_assists)] + \
                                      [uniform(LAST_LEG_DURATION_BIAS * 1.0, LAST_LEG_DURATION_BIAS * 100.0)]
        _rng_duration_of_legs_days = np.array(_rng_durations_of_legs_temp) / np.sum(
            _rng_durations_of_legs_temp) * duration_days
        sequence_planets = [sequence_planets[i][1:].replace('0000', str(int(_rng_duration_of_legs_days[i])).zfill(4))
                            for i in
                            range(len(sequence_planets))]
        sequence_planets = [' ' + seq for seq in sequence_planets]
        chromosone = chromosone.replace(' 00000' * len(_rng_planets), ''.join(sequence_planets), 1)

        # 0000 00000 00000 (0000)
        temp = chromosone.split(' ')
        temp[-1] = temp[-1].replace('0000', str(int(_rng_duration_of_legs_days[-1])).zfill(4), 1)
        chromosone = ' '.join(temp)

        return chromosone

    def mutation(self, chromosone):
        pass


class Population(object):

    @staticmethod
    def random_durations(size):
        # Sample a total trajectory time from a normal distribution.
        range = np.arange(-4, 4, 0.1)
        rangeU, rangeL = range + 0.5, range - 0.5
        prob = ss.norm.cdf(rangeU, scale=2) - ss.norm.cdf(rangeL, scale=2)
        prob = prob / prob.sum()
        nums = np.random.choice(range, size=size, p=prob)
        nums = 20 + nums
        return nums

    def __init__(self, _chromosone, _population_size):
        # Chromosone singleton and schema format.
        self._random_durations = self.random_durations(_population_size)
        self._random_durations_population = [choice(self._random_durations) for _ in range(_population_size)]
        self._population_size = _population_size
        self._chromosone = _chromosone
        self._unary_schema = _chromosone.unary_schema
        self._total_schema = _chromosone.total_schema
        self._genesis_generation = [self._chromosone.random_chromosone(duration)
                                    for duration in self._random_durations_population]
        self._current_generation = self._genesis_generation

    @property
    def genesis_generation(self):
        return self._genesis_generation

    @property
    def current_generation_fitness(self):
        return [self._chromosone.fitness(_chromosone, fitness_function) for _chromosone in self._current_generation]

    # @property
    # def current_generation(self):
    #     current_generation()
    #     while True:

    @property
    def chromosone(self):
        return self._chromosone


# class EvolutionaryAlgorithim(object):
#     def __init__(self, _chromosone, _fitness_function):
#         # Chromosone singleton and schema format.
#         self._chromosone = _chromosone
#         self._unary_schema = _chromosone.unary_schema
#         self._total_schema = _chromosone.total_schema
#
#         # Fitness function.
#         self._fitness_function = _fitness_function


if __name__ == '__main__':
    tt = TrajectoryTool()

    # Chromosone singleton setup.
    _unary_schema = list('123456789')
    _total_schema = '0000 00000 00000 0000'
    Chromosone = Chromosone(_unary_schema=_unary_schema,
                            _total_schema=_total_schema)

    # Population singleton setup.
    # _population_size = 100000
    # Population = Population(Chromosone, _population_size=_population_size)

    raw, bodyl = Chromosone.mapper('6262 62990 00000 5367')
    # print(raw)
    # print(bodyl)
    tt.process_itinerary(raw, bodyl, _mode='plot3D')

    # code = np.array(Population.genesis_generation)
    # results = np.array(Population.current_generation_fitness)
    #
    # print(results)
    #
    # print(results[results > 0])
    # print(code[results > 0])
