from trajectory_tool.mnag_mission_analysis.config import *
from trajectory_tool.mnag_mission_analysis.tools.orbital_mechanics import *
from trajectory_tool.mnag_mission_analysis.tools.plotting import *
from trajectory_tool.mnag_mission_analysis.planetary_flyby import PlanetaryFlyby
from trajectory_tool.mnag_mission_analysis.planetary_node import PlanetaryNode
from trajectory_tool.mnag_mission_analysis.planetary_departure import PlanetaryDeparture
from trajectory_tool.mnag_mission_analysis.planetary_rendezvous import PlanetaryRendezvous
# from trajectory_tool.genetic_algorithim_analysis.genetic_algorithim import Chromosome
from datetime import datetime, timedelta
from poliastro import iod
from tabulate import tabulate
plotly.tools.set_credentials_file(username='GeoffreyGarrett', api_key='9ynAL0uAu3iRs8ET8Zdb')
import pandas as pd


class InterplanetaryTrajectory(object):
    def __init__(self):
        self._id = None
        self._planetary_nodes = []
        self._launch_date = None
        self._durations = None
        self._planetary_departure = None
        self._planetary_flyby = None
        self._planetary_rendezvous = None

    @property
    def planetary_departure(self):
        return self._planetary_departure

    @property
    def planetary_flyby(self):
        return self._planetary_flyby

    @property
    def planetary_rendezvous(self):
        return self._planetary_rendezvous

    def _multi_leg_lambert_solution(self):
        for i in range(len(self._planetary_nodes) - 1):
            node_departure = self._planetary_nodes[i]
            node_arrival = self._planetary_nodes[i+1]
            (v0, v1), = iod.lambert(node_departure.body.parent.k,
                                    r0=node_departure.soi_exit_position_heliocentric,
                                    r=node_arrival.soi_entry_position_heliocentric,
                                    tof=time.Time(node_arrival.epoch_entry, scale='tdb') -
                                    time.Time(node_departure.epoch_exit, scale='tdb'))
            node_departure.v_exit = v0
            node_arrival.v_entry = v1

    def process_itinerary(self, itinerary):
        temp_planetary_nodes = itinerary['planetary_nodes']
        self._launch_date = itinerary['launch_date']
        self._durations = itinerary['durations']
        self._id = itinerary['id']
        epoch_periapses = [self._launch_date] + [self._launch_date + timedelta(days=365*sum(self._durations[:i+1]))
                                                 for i in range(len(self._durations))]
        for i, node in enumerate(temp_planetary_nodes):
            planetary_node = PlanetaryNode()
            planetary_node.body = body_list[node]
            planetary_node.epoch_periapsis = epoch_periapses[i]
            self._planetary_nodes.append(planetary_node)
        self._planetary_departure = PlanetaryDeparture(self._planetary_nodes[0])
        self._planetary_flyby = [PlanetaryFlyby(node) for node in self._planetary_nodes[1:-1]]
        self._planetary_rendezvous = PlanetaryRendezvous(self._planetary_nodes[-1])

    def check_feasibility(self):
        self._multi_leg_lambert_solution()
        for flyby in self.planetary_flyby:
            flyby.check_gravity_assist(flyby.planetary_node.v_entry, flyby.planetary_node.v_exit)

    def basic_analysis(self):
        self._multi_leg_lambert_solution()
        # self.planetary_departure.calculate_basic()
        for flyby in self.planetary_flyby:
            flyby.basic_powered_gravity_assist(flyby.planetary_node.v_entry, flyby.planetary_node.v_exit)
        # self.planetary_rendezvous.calculate_basic()

    def refined_analysis(self):
        boundary_error_velocity = 10                    # random
        boundary_error_position = 10                    # random
        while (boundary_error_velocity >= 10 ** (-7)) or (boundary_error_position >= 10 ** (-5)):    # 1 mm/s / 1 km
            self._multi_leg_lambert_solution()
            # self.planetary_departure.calculate_refined()
            for flyby in self.planetary_flyby:
                flyby.refine_powered_gravity_assist(flyby.planetary_node.v_entry, flyby.planetary_node.v_exit)
            # self.planetary_rendezvous.calculate_refined()
            try:
                boundary_error_velocity = sum([flyby.refined_attributes.error_v for flyby in self.planetary_flyby])
                boundary_error_position = sum([flyby.refined_attributes.error_p for flyby in self.planetary_flyby])
            except TypeError:
                print('hi')
                pass
            print('error_v: {:0.2f} km/s'.ljust(30).format(boundary_error_velocity))
            print('error_p: {:0.2f} km'.ljust(30).format(boundary_error_position))

    def plot3D(self, title=None, flyby=False, interplanetary=False):
        if flyby:
            plot_planetary_flyby(self.planetary_flyby)
        if interplanetary:
            plot_propagation(self)

    @property
    def temp_dv(self):
        launch_dv = self.planetary_departure.delta_v
        flyby_dv = self.planetary_flyby[0].basic_attributes.r_p_dv
        rendezvous = self.planetary_rendezvous.delta_v
        return launch_dv + flyby_dv + rendezvous


if __name__ == '__main__':
    START_EPOCH = datetime(2025, 1, 1, 0, 0, 0, 0)
    planet_mapping = {'1': Mercury,
                      '2': Venus,
                      '3': Earth,
                      '4': Mars,
                      '5': Jupiter,
                      '6': Saturn,
                      '7': Uranus,
                      '8': Neptune,
                      '9': Pluto}

    def mapper(chromosome):
        sequences = chromosome.split(' ')

        assists = list(filter(lambda a: a != '00000', sequences[1:-1]))

        _id = chromosome
        _launch_day = START_EPOCH + timedelta(days=int(sequences[0]))
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


    __itinerary = {'id': 1,
                   'planetary_nodes': ['earth', 'jupiter', 'pluto'],
                   'launch_date': datetime(2027, 12, 1, 0, 0),
                   'durations': [6, 8]}

    iten = mapper('1057 50756 00000 7998')

    # TODO: Fix chromosome mapper for new itinerary format

    ejp = InterplanetaryTrajectory()
    ejp.process_itinerary(iten)
    ejp.basic_analysis()
    ejp.refined_analysis()
    ejp.planetary_flyby[0].guess_powered_gravity_assist(ejp.planetary_flyby[0].planetary_node.v_entry,
                                                        ejp.planetary_flyby[0].planetary_node.v_exit)

    with pd.option_context('display.max_rows', 100, 'display.max_columns', 100, 'display.width', 300):
        np.set_printoptions(precision=3)
        # print(tabulate(ejp.planetary_flyby[0].basic_dataframe, headers='keys', tablefmt='psql', floatfmt=".2f"))
        print(ejp.planetary_flyby[0].guess_dataframe)
        print(ejp.planetary_flyby[0].refined_dataframe)
        # print(ejp.planetary_flyby[0].refined_dataframe)
        # ejp.plot3D(interplanetary=True, flyby=True)


        # data = ejp.planetary_flyby[0].refined_attributes
        #
        # op = OrbitPlotter3D()
        #
        # ss1 = Orbit.from_vectors(Sun, r=ejp.planetary_departure.planetary_node.soi_exit_position_heliocentric,
        #                               v=ejp.planetary_departure.planetary_node.v_exit)
        #
        # ss1.propagate(time.TimeDelta(u.s * (ejp.planetary_flyby[0].planetary_node.epoch_entry - ejp.planetary_departure.planetary_node.epoch_exit).total_seconds()))
        #
        # jup_i = Orbit.from_body_ephem(Jupiter, epoch=time.Time(ejp.planetary_flyby[0].planetary_node.epoch_entry))
        #
        # r_entry_j = ss1.state.r - jup_i.state.r
        # v_entry_j = ss1.state.v - jup_i.state.v
        #
        # ss2 = Orbit.from_vectors(Jupiter, r=r_entry_j, v=v_entry_j)
        # print(data)
        # ss2.propagate(time.TimeDelta(u.s * abs(data.t_i)))
        #
        # from poliastro.maneuver import Maneuver
        #
        # dv = data.v_p_f - data.v_p_i
        #
        # man = Maneuver.impulse(dv)
        #
        # ss2.apply_maneuver(man)
        #
        # ss2.propagate(time.TimeDelta(u.s * data.t_f))
        #
        # jup_f = Orbit.from_body_ephem(Jupiter, epoch=time.Time(ejp.planetary_flyby[0].planetary_node.epoch_exit))
        #
        # r_exit_j = ss2.state.r + jup_f.state.r
        # v_exit_j = ss2.state.v + jup_f.state.v
        #
        # ss3 = Orbit.from_vectors(Sun, r=r_exit_j, v=v_exit_j, epoch=time.Time(ejp.planetary_flyby[0].planetary_node.epoch_exit))
        #
        # ssp = Orbit.from_body_ephem(Pluto, epoch=time.Time(ejp.planetary_rendezvous.planetary_node.epoch_entry))
        #
        # op.plot(ss1)
        # op.plot(ss3)
        # op.plot(jup_f)
        # op.plot(jup_i)
        # op.plot(ssp)
        #
        # layout = go.Layout(title="test", width=800, height=800)
        # fig = go.Figure(data=op._data, layout=layout)
        # plotly.plotly.plot(fig)