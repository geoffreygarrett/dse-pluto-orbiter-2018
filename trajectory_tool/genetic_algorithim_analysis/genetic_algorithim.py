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
SIMILARITY_FILTER = 1.0
# ELITE_QUANTITY = 1.0
CROSSOVER_THRESHOLD = 0.30

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
        results = tt.process_itinerary(_raw_itinerary, _body_list, _mode='delta_v', verbose=False)
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
    def ensure_ordered(_chromosone):
        sequence = _chromosone.split(' ')
        start = sequence[0]
        end = sequence[-1]
        re_seq = []
        cnt = 0
        for seq in sequence[1:-1]:
            if seq[0] == '0':
                cnt += 1
                continue
            else:
                re_seq.append(seq)
        pad = ['00000'] * cnt
        return ' '.join([start]+re_seq+pad+[end])

    @staticmethod
    def similarity(chromosome1, chromosome2):
        return difflib.SequenceMatcher(a=chromosome1, b=chromosome2).ratio()

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
        sequence1 = chromosome1.split(' ')
        sequence2 = chromosome2.split(' ')

        array = [sequence1, sequence2]
        transpose = [1,0]
        idx_child1 = [randint(0,1) for _ in range(len(sequence1))]
        idx_child2 = [transpose[i] for i in idx_child1]

        child1 = [array[idx_child1[i]][i] for i, j in enumerate(idx_child1)]
        child2 = [array[idx_child2[i]][i] for i, j in enumerate(idx_child2)]

        return self.ensure_ordered(' '.join(child1)), self.ensure_ordered(' '.join(child2))

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

    def _similar_filtration(self, generation_df, percentage_similarity):
        idx_unique = 0
        chromosome_list = generation_df['Chromosome'].tolist()
        drop_idx = []
        count=0
        all_drop_idx = []
        for idx, chromo in enumerate(chromosome_list):
            unique = chromosome_list[idx_unique]
            if count is 0:
                similarity = 0
            else:
                similarity = self._chromosome.similarity(unique, chromo)
            count+=1
            if similarity >= percentage_similarity:
                # print('Unique: {}, Similarity: {:0.2f}, Dropping: {}'.format(unique, similarity, generation_df['Chromosome'][idx]))
                self.add([self._chromosome.random_chromosome(self.random_durations(1))])
                drop_idx.append(idx)
            elif len(drop_idx)>=1:
                idx_unique = idx
                drop_idx = []
                all_drop_idx.append(drop_idx)
        if len(drop_idx) > 1:
            generation_df.drop(generation_df.index[[drop_idx[0], drop_idx[-1]]], inplace=True)
        return generation_df.reset_index(drop=True)

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

    def add(self, list_chromosomes):
        list_chromosomes_fitness = self.groups_fitness(list_chromosomes)
        df_chromosomes = pd.DataFrame()
        df_chromosomes['Fitness'] = list_chromosomes_fitness
        df_chromosomes['Chromosome'] = list_chromosomes
        self._current_generation.append(df_chromosomes)

    def refine(self):
        self._current_generation = self._current_generation.drop_duplicates()
        self._current_generation = self.elite_class(self._population_size)
        self._current_generation = self._similar_filtration(self._current_generation, SIMILARITY_FILTER)
        return self._current_generation

    def crossover(self):
        crossover_threshold = CROSSOVER_THRESHOLD
        chosen_breeders = [self._current_generation['Chromosome'].tolist()[i]
                           for i in range(int(crossover_threshold*self._population_size))]
        random_mates = [choice(self._current_generation['Chromosome']) for _ in range(len(chosen_breeders))]
        pairs = zip(chosen_breeders, random_mates)
        children_chromosomes_zip = [self._chromosome.crossover(ch1, ch2) for ch1, ch2 in pairs]
        children_chromosomes = list([*zip(*children_chromosomes_zip)])[0]
        return children_chromosomes

    def mutate(self):
        mutated_chromosomes = [self._chromosome.mutation(cromo) for cromo in
                               self.current_generation['Chromosome'].tolist()]

        mutated_fitness = self.groups_fitness(mutated_chromosomes)

        mutated_df = pd.DataFrame()
        mutated_df['Fitness'] = mutated_fitness
        mutated_df['Chromosome'] = mutated_chromosomes
        collected_df = self._current_generation.append(mutated_df).drop_duplicates()

        # self._current_generation = self._filter(collected_df, _population_size)
        self._current_generation = collected_df
        return self._current_generation

    def filter(self):
        return self.current_generation.nlargest(_population_size, columns='Fitness').reset_index(drop=True)

    def groups_fitness(self, list):
        gen_fit = [self._chromosome.fitness(_chromosome, fitness_function) for _chromosome in list]
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
        population._generations += 1
        return EvolutionaryAlgorithim._mutate_population(EvolutionaryAlgorithim._crossover_population(population))

    @staticmethod
    def _crossover_population(population):
        children = population.crossover()
        population.add(children)
        return population

    @staticmethod
    def save_generation(population, identity=None):
        if not identity:
            newpath = os.path.join(DIR_GA,'generations_{}_test'.format(str('0').zfill(4)))
            if not os.path.exists(newpath):
                os.makedirs(newpath)
            else:
                pass
            num = 0
            while True:
                if not os.path.exists(os.path.join(DIR_GA, newpath, 'gen_{}'.format(num))):
                    population.current_generation.to_csv(path_or_buf=os.path.join(DIR_GA, newpath, 'gen_{}'.format(num)))
                    break
                else:
                    num+=1
        else:
            newpath = os.path.join(DIR_GA,'generations_{}'.format(str(identity).zfill(4)))
            if not os.path.exists(newpath):
                os.makedirs(newpath)
            else:
                pass
                num=0
                while True:
                    if not os.path.exists(os.path.join(DIR_GA, newpath, 'gen_{}'.format(num))):
                        population.current_generation.to_csv(path_or_buf=os.path.join(DIR_GA, newpath, 'gen_{}'.format(num)))
                        break
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
    to_do = ['evolve', 'plot', 'stats', 'other']
    TO_DO = 1
    # 1.25 2642 50549 30658 7364
    # 0.37 2647 50519 21248 6931
    # - 1.37 8241 50411 61476 6857
    # 2.46 1847 50549 00000 8199
    # 0.69 2257 50698 20954 6998
    INSPECT = '2257 50698 20954 6998'

    # chromosome singleton setup.
    _unary_schema = list('123456789')
    _total_schema = '0000 00000 00000 0000'
    Chromosome = Chromosome(_unary_schema=_unary_schema,
                            _total_schema=_total_schema)

    # Population singleton setup.
    _population_size = POPULATION_SIZE
    Population = Population(Chromosome, _population_size=_population_size, assists='1')

    if to_do[TO_DO] is 'plot':
        raw, bodyl = Chromosome.mapper(INSPECT)
        results = tt.process_itinerary(raw, bodyl, _mode='plot3D')

    if to_do[TO_DO] is 'stats':
        raw, bodyl = Chromosome.mapper(INSPECT)
        results = tt.process_itinerary(raw, bodyl, _mode='delta_v')
        print([results[k]['dv'] for k in range(len(results))])
        print(sum([results[k]['dv'] for k in range(len(results))]))

    if to_do[TO_DO] is 'evolve':
        count = 0
        while Population.fittest[0] < 5:
            try:
                result, top = Population.fittest[0], Population.fittest[1]
                if result >= 2.4:
                    EvolutionaryAlgorithim.save_generation(Population, 't1')

                print('Gen: {}'.format(Population._generations).ljust(20),'Fitness: {:0.2f}'.format(result).ljust(20),
                      'Chromosome: {}'.format(top))
                EvolutionaryAlgorithim.evolve(Population)
                Population.refine()
                count += 1

            except (KeyboardInterrupt, SystemExit):
                print('bye!')
                EvolutionaryAlgorithim.save_generation(Population)
                raise
            # except:
            # report error and proceed
            # if count % 5 == 0:
            #     print(Population.current_generation)

    if to_do[TO_DO] is 'other':
        z1, z2 = Chromosome.crossover('1449 50662 00000 7999', '0000 00000 00000 0000')
        print(Chromosome.similarity(z1,z2))
        print(difflib.SequenceMatcher(a=z1,b=z2).ratio())
        # print(z1)
        # print(z2)
        # X1 = [pd.DataFrame.from_csv(os.path.join(DIR_GA,'generations_0000_test','gen_{}'.format(i)))['Fitness'].tolist()[0] for i in range(124)]
        # # X2 = [pd.DataFrame.from_csv(os.path.join(DIR_GA,'generations_0000_test','gen_{}'.format(i)))['Chromosome'].tolist()[0].split(' ')[0] for i in range(124)]
        # X2 = np.arange(0,124)
        # x = pd.DataFrame.from_csv(os.path.join(DIR_GA,'generations_0000_test','gen_70'))
        # X = x['Chromosome'].tolist()
        # X = [int(i.split(' ')[0]) for i in X]
        # # plt.plot(X, 15-np.array(x['Fitness'].tolist()))
        # plt.plot(X2,X1)
        # plt.show()

