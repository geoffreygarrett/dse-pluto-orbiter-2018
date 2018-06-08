from random import randint, choice, uniform, random
from datetime import datetime
from trajectory_tool.core import *
import scipy.stats as ss
import numpy as np
from functools import lru_cache

LAST_LEG_DURATION_BIAS = 1.5
START_EPOCH = datetime.datetime(2025, 1, 1, 0, 0, 0, 0)
MUTATION_RATE = 0.1
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
    try:
        results = tt.process_itinerary(_raw_itinerary, _body_list, _mode='fast')
        delta_v = sum(results[i]['dv'] for i in range(len(results)))
    except (ValueError, RuntimeError):
        delta_v = 200
    return 15.0 - delta_v


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
        self._generations = []
        self._current_generation = self._genesis_generation
        self._generations.append([self._current_generation])
        self._historical_fitness = []
        self._historical_fitness.append(self.current_generation_fitness)

    # def _print_population(self):
    #     print('-'*40)
    def elite_class(self, qty):
        test = np.array(self._historical_fitness[-1])
        temp = np.argpartition(-test, qty)
        result_args = temp[:qty]
        temp = np.partition(-test, qty)
        result = -temp[:qty]
        top = [self._current_generation[i] for i in result_args]
        return list(zip(result, top))

    def mutate(self):
        _next_generation = [self._chromosome.mutation(cromo) for cromo in self._current_generation]
        self._historical_fitness.append(self.current_generation_fitness)
        self._current_generation = _next_generation
        self._historical_fitness.append(self.current_generation_fitness)

        _contending_generation = self._current_generation + _next_generation
        _contending_scores = self._historical_fitness[-2] + self._historical_fitness[-1]
        _contending_order = np.array(_contending_scores).argsort()

        self._historical_fitness.append(list(np.array(_contending_scores)[_contending_order >= len(_contending_scores)
                                                                          - POPULATION_SIZE]))

        self._current_generation = list(np.array(_contending_generation)[_contending_order >= len(_contending_scores)
                                                                         - POPULATION_SIZE])
        return self._current_generation


    @property
    def fittest(self):
        return self.elite_class(1)

    @property
    def genesis_generation(self):
        return self._genesis_generation

    @property
    def current_generation_fitness(self):
        gen_fit = [self._chromosome.fitness(_chromosome, fitness_function) for _chromosome in self._current_generation]
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

    # print(Population.elite_class(10))
    # for result, top in Population.elite_class(10):
    #     print('{:0.2f}'.ljust(10).format(result), top)

    raw, bodyl = Chromosome.mapper('6786 57343 00000 9593')
    # # print(raw)
    # # print(bodyl)
    tt.process_itinerary(raw, bodyl, _mode='plot3D')

    # code = np.array(Population.genesis_generation)
    # results = np.array(Population.current_generation_fitness)
    # #
    # # print(results)
    # #
    # print(results[results > -20])
    # print(code[results > -20])

    # for i in range(1000):
    #     print(Chromosome.mutation('5679 58178 00000 9761'))

    # while Population.fittest[0][0] < 2:
    #     result, top = Population.fittest[0]
    #     print('Fitness: {:0.2f}'.format(result).ljust(20), 'Chromosome: {}'.format(top))
    #     EvolutionaryAlgorithim.evolve(Population)
