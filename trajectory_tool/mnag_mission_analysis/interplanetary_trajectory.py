from trajectory_tool.mnag_mission_analysis.config import *
from trajectory_tool.mnag_mission_analysis.tools.orbital_mechanics import *
from trajectory_tool.mnag_mission_analysis.planetary_flyby import PlanetaryFlyby
from trajectory_tool.mnag_mission_analysis.planetary_node import PlanetaryNode
from trajectory_tool.mnag_mission_analysis.planetary_departure import PlanetaryDeparture
from trajectory_tool.mnag_mission_analysis.planetary_rendezvous import PlanetaryRendezvous
from datetime import datetime, timedelta
from poliastro import iod
from tabulate import tabulate
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
        while (boundary_error_velocity >= 10 ** (-6)) and (boundary_error_position >= 10 ** (-3)):    # 1 mm/s / 1 km
            self._multi_leg_lambert_solution()
            self.planetary_departure.calculate_refined()
            for flyby in self.planetary_flyby:
                flyby.refined_powered_gravity_assist(flyby.planetary_node.v_entry, flyby.planetary_node.v_exit)
            self.planetary_rendezvous.calculate_refined()
            boundary_error_velocity = sum([flyby.refined_attributes.error_v for flyby in self.planetary_flyby])
            boundary_error_position = sum([flyby.refined_attributes.error_p for flyby in self.planetary_flyby])


if __name__ == '__main__':
    __itinerary = {'id': 1,
                   'planetary_nodes': ['earth', 'jupiter', 'pluto'],
                   'launch_date': datetime(2027, 12, 1, 0, 0),
                   'durations': [6, 8]}

    # TODO: Fix chromosome mapper for new itinerary format

    ejp = InterplanetaryTrajectory()
    ejp.process_itinerary(__itinerary)

    # ejp.check_feasibility()

    ejp.basic_analysis()
    # ejp.refined_analysis()

    with pd.option_context('display.max_rows', 100, 'display.max_columns', 100, 'display.width', 10000):
        np.set_printoptions(precision=3)

        print(tabulate(ejp.planetary_flyby[0].basic_dataframe, headers='keys', tablefmt='psql', floatfmt=".2f"))
        # print(ejp.planetary_flyby[0].data_frame)
        # print(ejp.planetary_flyby[0]._basic_attributes.r_p)



