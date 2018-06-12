"""
Created by Alejandro Daniel Noel
"""

from jplephem.spk import SPK  # https://pypi.org/project/jplephem/
import numpy as np
import networkx as nx
from functools import lru_cache
import matplotlib.pyplot as plt
import matplotlib.animation
from mpl_toolkits.mplot3d import Axes3D
import random
import os
import datetime
from trajectory_tool.ale_solar_system.time_utils import datetime_to_jd


class SolarSystem:
    base_path = os.path.dirname(os.path.realpath(__file__))
    naif_ids = {
        'SOLAR_SYSTEM_BARYCENTER': [0, 'plu055.bsp'],
        'MERCURY_BARYCENTER': [1, 'plu055.bsp'],
        'VENUS_BARYCENTER': [2, 'plu055.bsp'],
        'EARTH_BARYCENTER': [3, 'plu055.bsp'],
        'MARS_BARYCENTER': [4, 'plu055.bsp'],
        'JUPITER_BARYCENTER': [5, 'plu055.bsp'],
        'SATURN_BARYCENTER': [6, 'plu055.bsp'],
        'URANUS_BARYCENTER': [7, 'plu055.bsp'],
        'NEPTUNE_BARYCENTER': [8, 'plu055.bsp'],
        'PLUTO_BARYCENTER': [9, 'plu055.bsp'],
        'PLUTO': [9, 'plu055.bsp'],
        'SUN': [10, 'plu055.bsp'],
        'MOON': [301, 'plu055.bsp'],
        'EARTH': [399, 'plu055.bsp'],
        'MERCURY': [199, 'plu055.bsp'],
        'VENUS': [299, 'plu055.bsp'],
        'CHARON': [901, 'de432s.bsp'],
        'NIX': [902, 'de432s.bsp'],
        'HYDRA': [903, 'de432s.bsp'],
        'KERBEROS': [904, 'de432s.bsp'],
        'STYX': [905, 'de432s.bsp'],
        # 'PLUTO': [999, 'de432s.bsp'],
    }

    def __init__(self):
        naif_files = list(set([value[1] for value in self.naif_ids.values()]))
        self.kernels = {naif_file: SPK.open(os.path.join(self.base_path, naif_file)) for naif_file in naif_files}
        self.ref_map = {}
        for kernel_file, kernel in self.kernels.items():
            self.ref_map.update({tuple(pair): kernel_file for pair in kernel.pairs.keys()})
        self.adjacency = nx.Graph()
        self.adjacency.add_edges_from(list(map(tuple, self.ref_map.keys())))

    def id_for_body(self, body_name: str):
        """
        Returns the NAIF id number for the given name of the body/element.
        """
        body_name = body_name.strip().upper().replace(' ', '_')
        if body_name not in self.naif_ids.keys():
            if "BARYCENTER" not in body_name:
                if body_name + '_BARYCENTER' in self.naif_ids.keys():
                    body_name += '_BARYCENTER'
        return self.naif_ids[body_name][0]

    def element_for_id(self, naif_id: int):
        """
        Returns the name of the element for the given NAIF id number.
        """
        for element, id_info in self.naif_ids.items():
            if id_info[0] == naif_id:
                return element
        raise ValueError('NAIF id {} not found'.format(naif_id))

    @lru_cache()
    def __sequence_between_bodies(self, target_body, reference_body):
        """
        Returns the sequence of coordinate transforms from one body/element to another
        using the available ephemeris files.
        """
        return nx.shortest_path(self.adjacency, self.id_for_body(target_body), self.id_for_body(reference_body))

    def coordinates_of(self, body, time, ref_body='sun'):
        """
        Returns coordinates of requested body wrt to ref_body (default to Sun) in km
        """
        index_path = self.__sequence_between_bodies(body, ref_body)
        coordinates = np.zeros(3)
        for i in range(len(index_path) - 1):
            pair = (index_path[i], index_path[i+1])
            if pair in list(self.ref_map.keys()):
                coordinates += self.kernels[self.ref_map[pair]][pair[0], pair[1]].compute(time)[:3]
                continue

            pair = (pair[1], pair[0])
            if pair in list(self.ref_map.keys()):
                coordinates -= self.kernels[self.ref_map[pair]][pair[0], pair[1]].compute(time)[:3]
                continue
            else:
                raise RuntimeError("Pair {} ({}-{}) not found".format(pair, self.element_for_id(pair[0]), self.element_for_id(pair[1])))

        return coordinates

    def distance_between(self, element1, element2, time):
        """
        Returns distance between 2 celestial elements in km at the specified epoch
        """
        return np.linalg.norm(self.coordinates_of(element1, time, ref_body=element2))

    def plot(self, bodies, time, ref_body='sun'):
        """
        Plots a snapshot of the given bodies at the given epoch.
        A reference body or element can also be specified.
        """
        body_positions = np.array([self.coordinates_of(body, time, ref_body=ref_body) for body in bodies]).T
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.scatter(body_positions[0], body_positions[1], body_positions[2])
        plt.show()

    def animate(self, bodies, start_time, end_time, ref_body='sun',
                interval=0.1, title="Solar system", display_radius=None):
        """
        Animates the given bodies in the given interval.
        A reference body or element can also be specified.
        """
        body_positions = np.array([self.coordinates_of(body, start_time, ref_body=ref_body) for body in bodies]).T
        fig = plt.figure()
        # ax = plt.axes(projection='3d')
        ax = Axes3D(fig)
        graph = ax.scatter(body_positions[0], body_positions[1], body_positions[2],
                           c=['#%06X' % random.randint(0, 195**3-1) for _ in range(len(body_positions[0]))])
        title_obj = ax.set_title(title)

        if display_radius:
            ax.set_xlim3d(-display_radius, display_radius)
            ax.set_ylim3d(-display_radius, display_radius)
            ax.set_zlim3d(-display_radius, display_radius)

        def update_plot(time):
            nonlocal body_positions
            body_positions = np.array([self.coordinates_of(body, time, ref_body=ref_body) for body in bodies]).T
            graph._offsets3d = (body_positions[0], body_positions[1], body_positions[2])
            title_obj.set_text(title + ' ' + str(time))

        ani = matplotlib.animation.FuncAnimation(fig, update_plot, frames=np.arange(start_time, end_time, interval), interval=100)
        plt.show()


if __name__ == "__main__":
    solar_system = SolarSystem()

    solar_system.id_for_body('jupiter')

    def unit_vector(vector):
        return vector / np.linalg.norm(vector)


    def angle_between(vector1, vector2):
        v1_u = unit_vector(vector1)
        v2_u = unit_vector(vector2)
        return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

    print(np.linalg.norm(solar_system.coordinates_of('earth', datetime_to_jd(datetime.date(2018, 6, 6)))))
    exit()

    vec1e = solar_system.coordinates_of('earth', datetime_to_jd(datetime.date(2018, 6, 6)))
    vec2e = solar_system.coordinates_of('earth', datetime_to_jd(datetime.date(2018, 7, 6)))
    vec3e = np.cross(vec1e, vec2e)

    vec1p = solar_system.coordinates_of('charon', datetime_to_jd(datetime.date(2018, 6, 6)), ref_body='pluto')
    vec2p = solar_system.coordinates_of('charon', datetime_to_jd(datetime.date(2018, 6, 15)), ref_body='pluto')
    vec3p = np.cross(vec1p, vec2p)


    vec_sp = solar_system.coordinates_of('pluto', datetime_to_jd(datetime.date(2018, 6, 6)))

    print(unit_vector(vec3e))
    print(unit_vector(vec3p))
    print(unit_vector(vec_sp))

    print(np.degrees(angle_between(vec3e, vec3p)))
    print(np.degrees(angle_between(vec_sp, vec3p)))
    print(np.cross((1, 0, 0), (0.1, 1, 0)))
    print(np.cross((1, 0, 0), (-0.1, -1, 0)))


    # solar_system.plot(['earth', 'venus', 'sun'], 2459580.5)
    # solar_system.plot(['earth', 'venus', 'sun'], 2460737.5)
    # print(solar_system.distance_between('earth', 'sun', 2457061.5))
    # solar_system.animate(['pluto', 'charon', 'nix', 'hydra', 'kerberos', 'styx'], 2457061.5, 2457171.5,
    #                      ref_body='pluto barycenter', interval=.4, title="Pluto-Charon system",
    #                      display_radius=solar_system.distance_between('pluto barycenter', 'nix', 2457061.5))

    # solar_system.animate(['sun',
    #                       'mercury',
    #                       'venus',
    #                       'earth',
    #                       'MARS_BARYCENTER',
    #                       'JUPITER_BARYCENTER',
    #                       'SATURN_BARYCENTER',
    #                       'uranus barycenter',
    #                       'neptune barycenter',
    #                       'pluto'], 2457061.5, 2467171.5,
    #                       ref_body='sun', interval=10.0,
    #                       display_radius=solar_system.distance_between('pluto barycenter', 'sun', 2457061.5))

    # solar_system.animate(['earth', 'moon'], 2457061.5, 2457171.5, ref_body='earth barycenter',
    #                      display_radius=solar_system.distance_between('earth', 'moon', 2457061.5))

