from random import randint, choice, uniform, random
from datetime import datetime
from trajectory_tool.core import *
from trajectory_tool.genetic_algorithim_analysis import DIR_GA
import scipy.stats as ss
import os
import numpy as np
from functools import lru_cache
import pandas as pd
import difflib
import _thread

LAST_LEG_DURATION_BIAS = 1.5
START_EPOCH = datetime.datetime(2025, 1, 1, 0, 0, 0, 0)
MUTATION_RATE = 0.1
POPULATION_SIZE = 100
SIMILARITY_FILTER = 0.9

FIRST_LEG_LIMIT_UPPER = 5250
LAST_LEG_LIMIT_UPPER = 5860

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
def fitness_function(chromosome_singleton, chromosome, tt):
    assert chromosome_singleton.unary_schema == list('123456789')
    assert chromosome_singleton.total_schema == '0000 00000 00000 0000'

    _raw_itinerary, _body_list = chromosome_singleton.mapper(chromosome)
    _total_dur = sum(_raw_itinerary['durations'])
    try:
        results = tt.process_itinerary(_raw_itinerary, _body_list, _mode='delta_v')
        delta_v_legs = [results[i]['dv'] for i in range(len(results))]
        delta_v = sum(delta_v_legs)

        # Propulsion leg penalty total
        propulsion_total_penalty = 8

        # First leg penalty
        if delta_v_legs[0] > FIRST_LEG_LIMIT_UPPER/1000.:
            first_leg_penalty_factor = 3
        else:
            first_leg_penalty_factor = 0

        # Last leg penalty
        if delta_v_legs[-1] > LAST_LEG_LIMIT_UPPER/1000.:
            last_leg_penalty_factor = 1
        else:
            last_leg_penalty_factor = 0

        if _total_dur >= 24:
            duration_penalty = 1000

        else:
            duration_penalty = 0

    except (ValueError, RuntimeError):  # Gravity assist/ Lambert solution not possible.
        delta_v = 200
        duration_penalty = 0
        last_leg_penalty_factor = 0
        first_leg_penalty_factor = 0
        propulsion_total_penalty = 0

    return (15.0 - delta_v - duration_penalty)
            # - propulsion_total_penalty * (last_leg_penalty_factor + first_leg_penalty_factor))     # Propulsion req.


class Chromosome(object):

    @staticmethod
    def mutation(_chromosome):
        sequences = _chromosome.split(' ')
        seqs = []
        for seq in sequences:
            seq = list(seq)
            if len(seq) is 5:   # Gravity assist
                if seq[0] is '0':
                    seqs.append(''.join(seq))
                    continue

                ##############################################

                ##############################################

                for i, gene in enumerate(seq):
                    if random() <= MUTATION_RATE:
                        seq[i] = choice(_unary_schema)
                    else:
                        pass
                seqs.append(''.join(seq))

            elif len(seq) is 4:    # Arrival/departure
                for i, gene in enumerate(seq):
                    if random() <= MUTATION_RATE:
                        seq[i] = choice(_unary_schema)
                    else:
                        pass
                seqs.append(''.join(seq))

        return ' '.join(seqs)

    @staticmethod
    def similarity(chromosome1, chromosome2):
        difflib.SequenceMatcher(chromosome1, chromosome2).ratio()

    @staticmethod
    def mapper(chromosome):
        sequences = chromosome.split(' ')

        assists = list(filter(lambda a: a != '00000', sequences[1:-1]))

        _id = chromosome
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

    # @lru_cache(maxsize=500)
    def fitness(self, _chromosome, _fitness_function, limit=None):
        return fitness_function(self, _chromosome, self._tt)

    @property
    def unary_schema(self):
        return self._unary_schema

    @property
    def total_schema(self):
        return self._total_schema

    def crossover(self, chromosome1, chromosome2):
        pass   # TODO: Finish

    def random_chromosome(self, duration, assists='any'):
        duration_days = duration * 365

        chromosome = deepcopy(self._total_schema)
        # (0000) 00000 00000 0000
        chromosome = chromosome.replace('0000', str(randint(1, 9999)).zfill(4), 1)

        # 0000 (0)0000 (0)0000 0000
        if assists is 'any':
            _rng_qty_assists = randint(1, len(self._total_schema.split(' ')) - 2)

        elif assists is '1':
            _rng_qty_assists = 1

        elif assists is '2':
            _rng_qty_assists = 2

        elif assists is '3':
            _rng_qty_assists = 3

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
        chromosome = chromosome.replace(' 00000' * len(_rng_planets), ''.join(sequence_planets), 1)

        # 0000 00000 00000 (0000)
        temp = chromosome.split(' ')
        temp[-1] = temp[-1].replace('0000', str(int(_rng_duration_of_legs_days[-1])).zfill(4), 1)
        chromosome = ' '.join(temp)

        return chromosome


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


    @staticmethod
    def remove_similar(generation_df, percentage_similarity):
        pass
        # return [x['']]

    def _filter(self, generation_df, number, similarity_threshold=None):
        # for row in generation_df:
        #
        # # filtered0 = generation_df.apply(self._chromosome.similarity(generation_df[]))
        filtered1 = generation_df.nlargest(number, columns='Fitness').reset_index(drop=True)
        return filtered1

    def __init__(self, _chromosome, _population_size, assists='any'):
        # chromosome singleton and schema format.
        self._random_durations = self.random_durations(_population_size)
        self._random_durations_population = [choice(self._random_durations) for _ in range(_population_size)]
        self._population_size = _population_size
        self._chromosome = _chromosome
        self._unary_schema = _chromosome.unary_schema
        self._total_schema = _chromosome.total_schema
        self._genesis_generation = [self._chromosome.random_chromosome(duration, assists)
                                    for duration in self._random_durations_population]

        # self._current_generation = self._genesis_generation

        self._current_generation = pd.DataFrame(columns=['Fitness', 'Chromosome'])
        self._current_generation['Chromosome'] = self._genesis_generation
        self._current_generation['Fitness'] = self.current_generation_fitness
        self._current_generation = self._current_generation.nlargest(_population_size, columns='Fitness').reset_index(drop=True)

        # self.current_generation2.reindex_axis(sorted(self.current_generation2['Fitness'].index), axis=1)
        self.generation_stats = pd.DataFrame(columns=['Generation', 'Best', 'Fitness'])
        self._generations = 0

    def elite_class(self, qty):
        return self._current_generation.nlargest(qty, columns='Fitness').reset_index(drop=True)

    def mutate(self):
        mutated_chromosomes = [self._chromosome.mutation(cromo) for cromo in
                               self.current_generation['Chromosome'].tolist()]

        mutated_fitness = self.groups_fitness(mutated_chromosomes)

        mutated_df = pd.DataFrame()
        mutated_df['Fitness'] = mutated_fitness
        mutated_df['Chromosome'] = mutated_chromosomes
        collected_df = self._current_generation.append(mutated_df).drop_duplicates()

        self._current_generation = self._filter(collected_df, _population_size)
        return self._current_generation

    def filter(self):
        return self.current_generation.nlargest(_population_size, columns='Fitness').reset_index(drop=True)

    def groups_fitness(self, list):
        gen_fit = [self._chromosome.fitness(_chromosome, fitness_function) for _chromosome in
                   list]
        return gen_fit

    @property
    def fittest(self):
        return self.elite_class(1).values[0]

    @property
    def genesis_generation(self):
        return self._genesis_generation

    @property
    def current_generation_fitness(self):
        gen_fit = [self._chromosome.fitness(_chromosome, fitness_function) for _chromosome in
                   self._current_generation['Chromosome'].tolist()]
        return gen_fit

    @property
    def current_generation(self):
        return self._current_generation

    @property
    def chromosome_singleton(self):
        return self._chromosome


class EvolutionaryAlgorithim(object):
    @staticmethod
    def evolve(population):
        return EvolutionaryAlgorithim._mutate_population(EvolutionaryAlgorithim._crossover_population(population))

    @staticmethod
    def _crossover_population(population):
        return population

    @staticmethod
    def save_generation(population, identity=None):
        if identity:
            newpath = os.path.join(DIR_GA,'generations_{0}_test'.format(str(id).zfill(4)))
            if not os.path.exists(newpath):
                os.makedirs(newpath)
            else:
                pass
            num = 0
            while True:
                if not os.path.exists(os.path.join(DIR_GA, newpath, 'gen_{}'.format(num))):
                    population.current_generation.to_csv('gen_{}'.format(num))
                else:
                    num+=1
        else:
            newpath = os.path.join(DIR_GA,'generations_{}'.format(str(id).zfill(4)))
            if not os.path.exists(newpath):
                os.makedirs(newpath)
            else:
                pass
                num=0
                while True:
                    if not os.path.exists(os.path.join(DIR_GA, newpath, 'gen_{}'.format(num))):
                        population.current_generation.to_csv('gen_{}'.format(num))
                    else:
                        num+=1


    @staticmethod
    def _mutate_population(population):
        return population.mutate()

    @staticmethod
    def _filter_population(population):
        return population.filter()


if __name__ == '__main__':
    tt = TrajectoryTool()
    to_do = ['evolve', 'plot', 'stats']
    TO_DO = 2
    INSPECT = '6665 30667 61425 6667'

    # chromosome singleton setup.
    _unary_schema = list('123456789')
    _total_schema = '0000 00000 00000 0000'
    Chromosome = Chromosome(_unary_schema=_unary_schema,
                            _total_schema=_total_schema)

    # Population singleton setup.
    _population_size = POPULATION_SIZE
    Population = Population(Chromosome, _population_size=_population_size, assists='any')

    if to_do[TO_DO] is 'plot':
        raw, bodyl = Chromosome.mapper(INSPECT)
        results = tt.process_itinerary(raw, bodyl, _mode='plot3D')

    if to_do[TO_DO] is 'stats':
        raw, bodyl = Chromosome.mapper(INSPECT)
        results = tt.process_itinerary(raw, bodyl, _mode='delta_v')
        print([results[k]['dv'] for k in range(len(results))])
        print(sum([results[k]['dv'] for k in range(len(results))]))

    if to_do is 'evolve':
        count = 0
        while Population.fittest[0] < 5:
            # if count % 5 == 0:
                # print(Population.current_generation)
            result, top = Population.fittest[0], Population.fittest[1]
            if result >= 0:
                EvolutionaryAlgorithim.save_generation(Population)

            print('Fitness: {:0.2f}'.format(result).ljust(20), 'Chromosome: {}'.format(top))
            EvolutionaryAlgorithim.evolve(Population)
            count+=1
