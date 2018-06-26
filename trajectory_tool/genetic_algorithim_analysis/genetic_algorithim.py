from random import randint, choice, uniform, random
from datetime import datetime
from trajectory_tool.core2 import *
from trajectory_tool.genetic_algorithim_analysis import DIR_GA
from trajectory_tool.mnag_mission_analysis.interplanetary_trajectory import InterplanetaryTrajectory
import scipy.stats as ss
import matplotlib
import os
import numpy as np
from functools import lru_cache
import pandas as pd
import difflib
import _thread

LAST_LEG_DURATION_BIAS = 1.5
START_EPOCH = datetime.datetime(2025, 1, 1, 0, 0, 0, 0)
SIMILARITY_FILTER = 0.75
CROSSOVER_THRESHOLD = 0.4

FIRST_LEG_LIMIT_UPPER = 6000
LAST_LEG_LIMIT_UPPER = 6000

planet_mapping = {'1': Mercury,
                  '2': Venus,
                  '3': Earth,
                  '4': Mars,
                  '5': Jupiter,
                  '6': Saturn,
                  '7': Uranus,
                  '8': Neptune,
                  '9': Pluto}

CROSSOVER_METHOD = 'old'

# pd.set_option('display.height', 1000)
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)


# Here is the definition of the fitness function.
def fitness_function(chromosome_singleton, chromosome, ft=None):
    assert chromosome_singleton.unary_schema == list('123456789')
    assert chromosome_singleton.total_schema == '0000 00000 00000 0000'

    _raw_itinerary = chromosome_singleton.mapper(chromosome)
    _total_dur = sum(_raw_itinerary['durations'])
    try:
        tt = InterplanetaryTrajectory()
        tt.process_itinerary(_raw_itinerary)
        tt.basic_analysis()
        delta_v = tt.temp_dv

        if list(chromosome.split(' ')[1])[0] is '3':
            earth_penalty = 1000

        else:
            earth_penalty = 0

        if _total_dur >= 24:
            duration_penalty = 1000

        else:
            duration_penalty = 0

    except (ValueError, RuntimeError):  # Gravity assist/ Lambert solution not possible.
        delta_v = 200
        duration_penalty = 0
        earth_penalty = 0
    return 15.0 - delta_v - duration_penalty - earth_penalty


class Chromosome(object):

    @staticmethod
    def mutation(_chromosome, mutation_rate):
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
                    if random() <= mutation_rate:
                        seq[i] = choice(_unary_schema)
                    else:
                        pass
                seqs.append(''.join(seq))

            elif len(seq) is 4:    # Arrival/departure
                for i, gene in enumerate(seq):
                    if random() <= mutation_rate:
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
        _raw_itinerary['planetary_nodes'] = _body_list
        return _raw_itinerary

    @staticmethod
    def mapper_old(chromosome):
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
        # _raw_itinerary['planetary_nodes'] = _body_list
        return _raw_itinerary, _body_list

    def __init__(self, _unary_schema, _total_schema, tt):
        self._unary_schema = _unary_schema
        self._total_schema = _total_schema
        self._tt = tt

    # @lru_cache(maxsize=500)
    def fitness(self, _chromosome, _fitness_function):
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
        transpose = [1, 0]
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
            count += 1
            if similarity >= percentage_similarity:
                self.add([self._chromosome.random_chromosome(self.random_durations(1))])
                drop_idx.append(idx)
            elif len(drop_idx) >= 1:
                idx_unique = idx
                drop_idx = []
                all_drop_idx.append(drop_idx)
        if len(drop_idx) > 1:
            generation_df.drop(generation_df.index[[drop_idx[0], drop_idx[-1]]], inplace=True)
        return generation_df.reset_index(drop=True)

    def __init__(self, _chromosome, _population_size, _mutation_rate, assists='any'):
        # chromosome singleton and schema format.
        self._random_durations = self.random_durations(_population_size)
        self._random_durations_population = [choice(self._random_durations) for _ in range(_population_size)]
        self._population_size = _population_size
        self._chromosome = _chromosome
        self._unary_schema = _chromosome.unary_schema
        self._total_schema = _chromosome.total_schema
        self._mutation_rate = _mutation_rate
        self._genesis_generation = [self._chromosome.random_chromosome(duration, assists)
                                    for duration in self._random_durations_population]

        # self._current_generation = self._genesis_generation

        self._current_generation = pd.DataFrame(columns=['Fitness', 'Chromosome', 'Fitness_ratio', 'Fitness_cum'])
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

    def replace_with(self, list_chromosomes):
        list_chromosomes_fitness = self.groups_fitness(list_chromosomes)
        df_chromosomes = pd.DataFrame()
        df_chromosomes['Chromosome'] = list_chromosomes
        df_chromosomes['Fitness'] = list_chromosomes_fitness
        self._current_generation = df_chromosomes

    @staticmethod
    def fitness_ratio(fitness_list):
        offset_at_zero = np.array(fitness_list) - np.min(fitness_list)
        generation_fitness_list = (1/np.max(offset_at_zero))*offset_at_zero
        fitness_ratio_list = generation_fitness_list/np.sum(generation_fitness_list)*100
        return fitness_ratio_list

    def refine(self):
        if CROSSOVER_METHOD is 'new':
            self._current_generation = self.elite_class(self._population_size)
            self._current_generation['Fitness_ratio'] = self.fitness_ratio(self._current_generation['Fitness'])
            self._current_generation['Fitness_cum'] = self._current_generation['Fitness_ratio'].cumsum()
        else:
            self._current_generation = self._current_generation.drop_duplicates()
            self._current_generation = self.elite_class(self._population_size)
            self._current_generation = self._similar_filtration(self._current_generation, SIMILARITY_FILTER)
        return self._current_generation

    @staticmethod
    def random_ratio(df):
        return df[df['Fitness_cum'] >= uniform(0, 100)]['Chromosome'].tolist()[0]

    def crossover(self):
        crossover_threshold = CROSSOVER_THRESHOLD
        fitness_ratio = self.fitness_ratio(self._current_generation['Fitness'])
        self._current_generation['Fitness_ratio'] = pd.Series(np.array(fitness_ratio),
                                                              index=self._current_generation.index)
        self._current_generation['Fitness_cum'] = self._current_generation['Fitness_ratio'].cumsum()
        chosen_breeders = [self.random_ratio(self._current_generation) for _ in range(int(0.5*self._population_size))]
        random_mates = [self.random_ratio(self._current_generation) for _ in range(int(0.5*self._population_size))]
        if CROSSOVER_METHOD is 'new':
            fitness_ratio = self.fitness_ratio(self._current_generation['Fitness'])
            self._current_generation['Fitness_ratio'] = pd.Series(np.array(fitness_ratio),
                                                                  index=self._current_generation.index)
            self._current_generation['Fitness_cum'] = self._current_generation['Fitness_ratio'].cumsum()
        else:
            chosen_breeders = [self._current_generation['Chromosome'].tolist()[i]
                               for i in range(int(crossover_threshold*self._population_size))]
            random_mates = [choice(self._current_generation['Chromosome']) for _ in range(len(chosen_breeders))]
        pairs = zip(chosen_breeders, random_mates)
        children_chromosomes_zip = [self._chromosome.crossover(ch1, ch2) for ch1, ch2 in pairs]
        children_chromosomes = list([*zip(*children_chromosomes_zip)])[0]
        return children_chromosomes

    def mutate(self):
        mutated_chromosomes = [self._chromosome.mutation(cromo, self._mutation_rate) for cromo in
                               self.current_generation['Chromosome'].tolist()]
        mutated_fitness = self.groups_fitness(mutated_chromosomes)
        mutated_df = pd.DataFrame()
        mutated_df['Fitness'] = mutated_fitness
        mutated_df['Chromosome'] = mutated_chromosomes
        collected_df = self._current_generation.append(mutated_df).drop_duplicates()
        self._current_generation = collected_df
        return self._current_generation

    def filter(self):
        return self.current_generation.nlargest(_population_size, columns='Fitness').reset_index(drop=True)

    def groups_fitness(self, list):
        gen_fit = [self._chromosome.fitness(_chromosome, fitness_function) for _chromosome in list]
        return gen_fit

    @property
    def fittest(self):
        return self._current_generation.iloc[[0]]

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
        if CROSSOVER_METHOD is 'new':
            children = population.crossover()
            population.replace_with(children)
        else:
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
    tt2 = InterplanetaryTrajectory()
    to_do = ['evolve', 'plot', 'stats', 'other']
    TO_DO = 0

    # 1.25 2642 50549 30658 7364
    # 0.37 2647 50519 21248 6931
    # 0.98 2247 50699 20963 6999
    # - 1.37 8241 50411 61476 6857
    # 2.46 1847 50549 00000 8199
    # 0.69 2257 50698 20954 6998
    INSPECT = '1057 50756 00000 7998'
    INSPECT = '1449 50662 00000 7999'

    # chromosome singleton setup.
    _unary_schema = list('123456789')
    _total_schema = '0000 00000 00000 0000'
    Chromosome = Chromosome(_unary_schema=_unary_schema,
                            _total_schema=_total_schema,
                            tt = tt2)

    if to_do[TO_DO] is 'plot':
        raw, bodyl = Chromosome.mapper(INSPECT)
        results = tt.stationary_process_itinerary(raw, bodyl, mode='plot3D')

    if to_do[TO_DO] is 'stats':
        raw, bodyl = Chromosome.mapper(INSPECT)
        results = tt.stationary_process_itinerary(raw, bodyl, mode='vector_evaluation')
        # print(sum(results[i]['dv'] for i in range(len(results))))
        vc = np.sqrt(Earth.k/(6378*u.km + 180*u.km))
        v_inf2 = ((np.linalg.norm(results[0]['v']['d']-results[0]['v']['p'])) * u.km / u.s)
        print(v_inf2)

        ans = vc*(np.sqrt(2+(v_inf2/vc)**2)-1)
        print('here fab: ', ans.to(u.km/u.s))

    if to_do[TO_DO] is 'evolve':
        # Population singleton setup.
        # MUTATION_RATE = 0.1
        # POPULATION_SIZE = 100
        dic = {'identifier': 'mumod2',
               'process': [(600, 0.025), (600, 0.075), (600, 0.15)]}

        for pop,mut in dic['process']:

            # _population_size = POPULATION_SIZE
            _Population = Population(Chromosome, _population_size=pop, _mutation_rate=mut, assists='1')
            count = 0
            gen = []
            fitness = []
            title = '{}0_P{}_M{}'.format(dic['identifier'], pop, mut)

            while True:
                # Population.fittest < 5:
                try:
                    result, top = _Population.fittest['Fitness'][0], _Population.fittest['Chromosome']
                    if _Population._generations >= 50:
                        print('next!')
                        EvolutionaryAlgorithim.save_generation(_Population)
                        np.save(title + '_gen', np.array(gen))
                        np.save(title + '_fitness', np.array(fitness))
                        break
                    if result >= 2.4:
                        EvolutionaryAlgorithim.save_generation(_Population, title)

                    print('Gen: {}'.format(_Population._generations).ljust(20),'Fitness: {:0.2f}'.format(result).ljust(20),
                          'Chromosome: {}'.format(top))
                    gen.append(_Population._generations)
                    fitness.append(result)
                    EvolutionaryAlgorithim.evolve(_Population)

                    _Population.refine()
                    count += 1

                except (KeyboardInterrupt, SystemExit):
                    print('bye!')
                    EvolutionaryAlgorithim.save_generation(_Population)
                    np.save(title+'_gen', np.array(gen))
                    np.save(title+'_fitness', np.array(fitness))
                    raise

    if to_do[TO_DO] is 'other':
        dir =   '/home/samuel/Local Software/Local Projects/pluto_orbiter_2018/trajectory_tool/'
        # plot_filenames = ['t3_P100_M0.05', 't3_P100_M0.1','t3_P100_M0.15', 't3_P100_M0.2', 't3_P100_M0.3']
        # for fn in plot_filenames:
        #     plt.plot(np.load(dir+fn+'_gen.npy')[0:100],np.load(dir+fn+'_fitness.npy')[0:100], label='P={}, M={}'.format(fn.split('_')[1].replace('P',''), fn.split('_')[2].replace('M','')))
        #
        # plt.legend()
        # plt.show()

        # def process_fitness(fitness_list):

        # plot_filenames = [['t4_P100_M0.1', 't5_P100_M0.1', 't6_P100_M0.1', 't8_P100_M0.1'],
        #                   ['t4_P100_M0.2', 't5_P100_M0.2', 't6_P100_M0.2', 't7_P100_M0.2', 't8_P100_M0.2'],
        #                   ['t4_P100_M0.3', 't5_P100_M0.3', 't6_P100_M0.3', 't7_P100_M0.3', 't8_P100_M0.3']]

        # plot_filenames = [['m1_P200_M0.1', 'm2_P200_M0.1'],
        #                   ['m1_P400_M0.1', 'm2_P400_M0.1'],
        #                   ['m1_P600_M0.1', 'm2_P600_M0.1']]

        # plot_filenames = [['popmod0_P200_M0.1', 'popmod1_P200_M0.1', 'popmod2_P200_M0.1', 'popmod3_P200_M0.1', 'popmod4_P200_M0.1', 'popmod5_P200_M0.1', 'popmod6_P200_M0.1', 'popmod7_P200_M0.1'],
        #                   ['popmod0_P400_M0.1', 'popmod1_P400_M0.1', 'popmod2_P400_M0.1', 'popmod3_P400_M0.1', 'popmod4_P400_M0.1', 'popmod5_P400_M0.1', 'popmod6_P400_M0.1', 'popmod7_P400_M0.1'],
        #                   ['popmod0_P600_M0.1', 'popmod1_P600_M0.1', 'popmod2_P600_M0.1', 'popmod3_P600_M0.1', 'popmod4_P600_M0.1', 'popmod5_P600_M0.1', 'popmod6_P600_M0.1', 'popmod7_P600_M0.1']]

        _type='bar'
        plot_filenames = [['mumod0_P100_M0.025', 'mumod1_P100_M0.025', 'mumod2_P100_M0.025', 'mumod3_P100_M0.025', 'mumod4_P100_M0.025', 'mumod5_P100_M0.025', 'mumod6_P100_M0.025', 'mumod7_P100_M0.025'],
                          ['mumod0_P100_M0.075', 'mumod1_P100_M0.075', 'mumod2_P100_M0.075', 'mumod3_P100_M0.075', 'mumod4_P100_M0.075', 'mumod5_P100_M0.075', 'mumod6_P100_M0.075', 'mumod7_P100_M0.075'],
                          ['mumod0_P100_M0.15', 'mumod1_P100_M0.15', 'mumod2_P100_M0.15', 'mumod3_P100_M0.15', 'mumod4_P100_M0.15', 'mumod5_P100_M0.15', 'mumod6_P100_M0.15', 'mumod7_P100_M0.15']] # plt.style.use('ggplot')
        plt.style.use(['default'])
        colors = ['red', 'green', 'blue']
        hatching = ['.', '/', '\\']
        _linestyle = ['dashed', '-', '-.']
        alpha = [0.20, 0.20, 0.09]
        # plt.subplot(311)]                                                                                           p

        if _type is 'convergence':
            plt.figure(num=None, figsize=(8, 6), dpi=150, facecolor='w', edgecolor='k')
            for idx,fn in enumerate(plot_filenames):
                g1 = np.load(dir + fn[0] + '_gen.npy')[0:50]
                g2 = np.load(dir + fn[1] + '_gen.npy')[0:50]
                g3 = np.load(dir + fn[2] + '_gen.npy')[0:50]
                g4 = np.load(dir + fn[3] + '_gen.npy')[0:50]
                g5 = np.load(dir + fn[4] + '_gen.npy')[0:50]
                g6 = np.load(dir + fn[5] + '_gen.npy')[0:50]
                g7 = np.load(dir + fn[6] + '_gen.npy')[0:50]
                g8 = np.load(dir + fn[7] + '_gen.npy')[0:50]

                f1 = np.load(dir + fn[0] + '_fitness.npy')[0:50]
                f2 = np.load(dir + fn[1] + '_fitness.npy')[0:50]
                f3 = np.load(dir + fn[2] + '_fitness.npy')[0:50]
                f4 = np.load(dir + fn[3] + '_fitness.npy')[0:50]
                f5 = np.load(dir + fn[4] + '_fitness.npy')[0:50]
                f6 = np.load(dir + fn[5] + '_fitness.npy')[0:50]
                f7 = np.load(dir + fn[6] + '_fitness.npy')[0:50]
                f8 = np.load(dir + fn[7] + '_fitness.npy')[0:50]



                f_1 = f1 - np.min(f1)
                f_2 = f2 - np.min(f2)
                f_3 = f3 - np.min(f3)
                f_4 = f4 - np.min(f4)
                f_5 = f5 - np.min(f5)
                f_6 = f6 - np.min(f6)
                f_7 = f7 - np.min(f7)
                f_8 = f8 - np.min(f8)
                # f_3 = f3 - np.min(f3)
                # f_4 = f4 - np.min(f4)

                f_1 = f_1 / np.max(f_1) * 100
                f_2 = f_2 / np.max(f_2) * 100
                f_3 = f_3 / np.max(f_3) * 100
                f_4 = f_4 / np.max(f_4) * 100
                f_5 = f_5 / np.max(f_5) * 100
                f_6 = f_6 / np.max(f_6) * 100
                f_7 = f_7 / np.max(f_7) * 100
                f_8 = f_8 / np.max(f_8) * 100
                # f_3 = f_3 / np.max(f_3) * 100
                # f_4 = f_4 / np.max(f_4) * 100

                # y = np.array([np.mean([i,j,k,m]) for i,j,k,m in zip(f_1, f_2, f_3, f_4)])
                # e = np.array([np.std([f_1[i], f_2[i], f_3[i], f_4[i]]) for i in range(len(g3))])

                y = np.array([np.mean([i,j,l,m,m,o,p,q]) for i,j,l,m,n,o,p,q in zip(f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8)])
                e = np.array([2*np.std([f_1[i], f_2[i], f_3[i], f_4[i], f_5[i], f_6[i], f_7[i], f_8[i]]) for i in range(len(g2))])

                from scipy.interpolate import interp1d

                x_new = np.arange(0, 49, 0.2)
                # plt.subplot(3,1,idx+1)
                print(y-e)
                y_error_minus = interp1d(g2, y-e, kind='cubic')(x_new)

                print(y_error_minus)
                y_error_plus = interp1d(g2, y+e, kind='cubic')(x_new)
                # y = [np.mean([i,j,k]) for i,j,k in (f_1, f_2, f_3)]

                # plt.errorbar(g3, y, yerr=e, fmt='o', color=colors[idx], linewidth=0.5, ls=_linestyle[idx])
                plt.plot(g2, y, color=colors[idx], linestyle=_linestyle[idx])
                plt.fill_between(x_new, y_error_minus, y_error_plus, alpha=alpha[idx], color=colors[idx], hatch=hatching[idx], label='P={}, M={} | 2$\sigma$'.format(fn[0].split('_')[1].replace('P',''), fn[0].split('_')[2].replace('M','')))
                plt.plot(g2, y, linestyle=_linestyle[idx], color=colors[idx], label='P={}, M={} | mean'.format(fn[0].split('_')[1].replace('P',''), fn[0].split('_')[2].replace('M','')))
                # plt.legend()
            plt.grid(b=True, which='both', color='0.65', linestyle='-')
            plt.xlabel('Generation', fontsize=18)
            plt.ylabel('Fitness Ratio [%]', fontsize=18)
            plt.legend(fontsize=15)
            plt.show()

        if _type is 'bar':
            means = []
            stds=[]
            group_names = []
            # plt.figure(num=None, figsize=(8, 6), dpi=150, facecolor='w', edgecolor='k')
            # plt.grid(b=True, which='both', color='0.65', linestyle='-')
            for idx, fn in enumerate(plot_filenames):
                g1 = np.load(dir + fn[0] + '_gen.npy')[0:50]
                g2 = np.load(dir + fn[1] + '_gen.npy')[0:50]
                g3 = np.load(dir + fn[2] + '_gen.npy')[0:50]
                g4 = np.load(dir + fn[3] + '_gen.npy')[0:50]
                g5 = np.load(dir + fn[4] + '_gen.npy')[0:50]
                g6 = np.load(dir + fn[5] + '_gen.npy')[0:50]
                g7 = np.load(dir + fn[6] + '_gen.npy')[0:50]
                g8 = np.load(dir + fn[7] + '_gen.npy')[0:50]

                f1 = np.load(dir + fn[0] + '_fitness.npy')[0:50]
                f2 = np.load(dir + fn[1] + '_fitness.npy')[0:50]
                f3 = np.load(dir + fn[2] + '_fitness.npy')[0:50]
                f4 = np.load(dir + fn[3] + '_fitness.npy')[0:50]
                f5 = np.load(dir + fn[4] + '_fitness.npy')[0:50]
                f6 = np.load(dir + fn[5] + '_fitness.npy')[0:50]
                f7 = np.load(dir + fn[6] + '_fitness.npy')[0:50]
                f8 = np.load(dir + fn[7] + '_fitness.npy')[0:50]

                f_1 = np.max(f1)
                f_2 = np.max(f2)
                f_3 = np.max(f3)
                f_4 = np.max(f4)
                f_5 = np.max(f5)
                f_6 = np.max(f6)
                f_7 = np.max(f7)
                f_8 = np.max(f8)
                # f_3 = f3 - np.min(f3)
                # f_4 = f4 - np.min(f4)
                print(f_1)
                y = np.mean([f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8])
                std2 = np.std([f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8])

                means.append(y)
                stds.append(std2)
                group_names.append('P={}\nM={}'.format(fn[0].split('_')[1].replace('P',''), fn[0].split('_')[2].replace('M','')))

            fig = plt.figure(num=None, figsize=(8, 6), dpi=150, facecolor='w', edgecolor='k')
            plt.grid(b=True, which='both', color='0.65', linestyle='-')
            ax = fig.add_subplot(111)

            ## the data
            N = 3
            # = [18, 35, 30, 35, 27]
            # menStd = [2, 3, 4, 1, 2]
            # womenMeans = [25, 32, 34, 20, 25]
            # womenStd = [3, 5, 2, 3, 3]

            ## necessary variables
            ind = np.arange(N)  # the x locations for the groups
            width = 0.35  # the width of the bars
            print(means, width)
            ## the bars
            rects1 = ax.bar(ind, means, width,
                            color='blue',
                            yerr=stds,
                            error_kw=dict(elinewidth=5, ecolor='red'),
                            alpha=0.6)

            # rects2 = ax.bar(ind + width, womenMeans, width,
            #                 color='red',
            #                 yerr=womenStd,
            #                 error_kw=dict(elinewidth=2, ecolor='black'))

            # axes and labels
            # ax.set_xlim(-width, len(ind) + width)
            # ax.set_ylim(0, 15)
            # ax.set_ylabel('')
            xTickMarks = group_names
            ax.set_xticks(ind)
            xtickNames = ax.set_xticklabels(xTickMarks, fontsize=15)
            # plt.setp(xtickNames, rotation=45, fontsize=10)

            ## add a legend
            # ax.legend((rects1[0]), ('Men'))
            # plt.xlabel('Generation', fontsize=18)
            plt.ylabel('Maximum fitness f(C)', fontsize=18)
            # plt.legend(fontsize=15)
            plt.show()



