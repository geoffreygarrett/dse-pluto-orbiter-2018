from random import randint, choice, uniform, random
from datetime import datetime
from trajectory_tool.core import *
import scipy.stats as ss
import numpy as np
from functools import lru_cache
import pandas as pd

LAST_LEG_DURATION_BIAS = 1.5
START_EPOCH = datetime.datetime(2025, 1, 1, 0, 0, 0, 0)
MUTATION_RATE = 0.2
POPULATION_SIZE = 100

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
        results = tt.process_itinerary(_raw_itinerary, _body_list, _mode='fast')
        delta_v = sum(results[i]['dv'] for i in range(len(results)))

        if _total_dur >= 24:
            duration_penalty = 1000

        else:
            duration_penalty = 0

    except (ValueError, RuntimeError):
        delta_v = 200
        duration_penalty = 0

    return 15.0 - delta_v - duration_penalty


class Chromosome(object):

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

    @lru_cache(maxsize=500)
    def fitness(self, _chromosome, _fitness_function, limit=None):
        return fitness_function(self, _chromosome, self._tt)

    @property
    def unary_schema(self):
        return self._unary_schema

    @property
    def total_schema(self):
        return self._total_schema

    def crossover(self, chromosome1, chromosome2):
        pass

    def mutation(self, _chromosome):
        # sequences = _chromosome.splut(' ')
        # assist_seq = sequences[1:-1]
        sequences = _chromosome.split(' ')
        seqs = []
        for seq in sequences:
            seq = list(seq)
            if len(seq) is 5:   # Gravity assist
                if seq[0] is '0':
                    seqs.append(''.join(seq))
                    continue
                #####

                ######
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


    def random_chromosome(self, duration):
        duration_days = duration * 365

        chromosome = deepcopy(self._total_schema)
        # (0000) 00000 00000 0000
        chromosome = chromosome.replace('0000', str(randint(1, 9999)).zfill(4), 1)

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
    def _filter(generation_df, number):
        return generation_df.nlargest(number, columns='Fitness').reset_index(drop=True)

    def __init__(self, _chromosome, _population_size):
        # chromosome singleton and schema format.
        self._random_durations = self.random_durations(_population_size)
        self._random_durations_population = [choice(self._random_durations) for _ in range(_population_size)]
        self._population_size = _population_size
        self._chromosome = _chromosome
        self._unary_schema = _chromosome.unary_schema
        self._total_schema = _chromosome.total_schema
        self._genesis_generation = [self._chromosome.random_chromosome(duration)
                                    for duration in self._random_durations_population]

        # self._current_generation = self._genesis_generation

        self._current_generation = pd.DataFrame(columns=['Fitness', 'Chromosome'])
        self._current_generation['Chromosome'] = self._genesis_generation
        self._current_generation['Fitness'] = self.current_generation_fitness
        self._current_generation = self._current_generation.nlargest(_population_size, columns='Fitness').reset_index(drop=True)

        # self.current_generation2.reindex_axis(sorted(self.current_generation2['Fitness'].index), axis=1)
        self.generation_stats = pd.DataFrame(columns=['Generation', 'Best', 'Fitness'])
        self._generations = 0
        print(self._current_generation)

    def elite_class(self, qty):
        return self._current_generation.nlargest(qty, columns='Fitness').reset_index(drop=True)

    def mutate(self):
        mutated_chromosomes = [self._chromosome.mutation(cromo) for cromo in
                               self.current_generation['Chromosome'].tolist()]

        mutated_fitness = self.groups_fitness(mutated_chromosomes)

        mutated_df = pd.DataFrame()
        mutated_df['Fitness'] = mutated_fitness
        mutated_df['Chromosome'] = mutated_chromosomes

        self._current_generation = self._filter(self._current_generation.append(mutated_df), _population_size)
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
    def _mutate_population(population):
        return population.mutate()

    @staticmethod
    def _filter_population(population):
        return population.filter()


if __name__ == '__main__':
    tt = TrajectoryTool()

    # chromosome singleton setup.
    _unary_schema = list('123456789')
    _total_schema = '0000 00000 00000 0000'
    Chromosome = Chromosome(_unary_schema=_unary_schema,
                            _total_schema=_total_schema)

    # Population singleton setup.
    _population_size = POPULATION_SIZE
    Population = Population(Chromosome, _population_size=_population_size)

    raw, bodyl = Chromosome.mapper('1449 50662 00000 7996')
    results = tt.process_itinerary(raw, bodyl, _mode='full')
    print([results[i]['dv'] for i in range(len(results))])

    print([results[i]['v']['p'] for i in range(len(results))])
    print([results[i]['v']['a'] for i in range(len(results))])
    #
    # while Population.fittest[0] < 4:
    #     result, top = Population.fittest[0], Population.fittest[1]
    #     print('Fitness: {:0.2f}'.format(result).ljust(20), 'Chromosome: {}'.format(top))
    #     EvolutionaryAlgorithim.evolve(Population)
