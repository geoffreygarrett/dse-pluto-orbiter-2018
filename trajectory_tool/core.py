# IMPORTS -------------------------------------------------------------------------------------------------------------
import matplotlib.pyplot as plt
from matplotlib import cm

from poliastro.plotting import OrbitPlotter
from poliastro.bodies import Earth, Mars, Sun, Pluto, Venus, Jupiter
from poliastro.twobody import Orbit
from poliastro import iod
from astropy.coordinates import solar_system_ephemeris
from astropy import time

import numpy as np
import datetime
from collections import namedtuple
lambert_parameters = namedtuple('lambert_parameters', 'r0 r1 v0 v1 tof attractor epoch0 epoch1 ss0 ss1')

plt.style.use("seaborn")
solar_system_ephemeris.set("jpl")

#
# # DECORATORS ---------------------------------------------------------------------------------------------------------
# def decorators(lambert_hohman=True):
#     class _Decorator(object):
#
#         def __init__(self, decorated):
#             self._decorated = decorated
#
#         def __get__(self, instance, owner):
#             self.cls = owner
#             self.obj = instance
#             return self.__call__
#
#         # Decorator 1
#         def _lambert_hohman(self):
#
#
#         def __call__(self, *args, **kwargs):
#             if lambert_hohman is True:
#                 self._handle_traffic()
#
#             if update_signature is True:
#                 self._update_signature()
#
#     return _Decorator

# def lambert_hohman(lambert_function):
#     soln = lambert_function(*args, **kwargs)


# CORE CLASS -----------------------------------------------------------------------------------------------------------
class TrajectoryTool(object):
    """

    """
    @staticmethod
    def lambert_solve_from_positions(r, r0, epoch0, epoch1, main_attractor=Sun):
        """

        :param r:
        :param r0:
        :param epoch0:
        :param epoch1:
        :param main_attractor:
        :return:
        """
        tof = epoch1 - epoch0
        (v0, v) = iod.lambert(main_attractor.K, r0, r, tof)

        return lambert_parameters(r0=r0, r1=r, v0=v0, v1=v, tof=tof,  attractor=main_attractor, epoch0=epoch0,
                                  epoch1=epoch1, ss0=None, ss1=None)

    @staticmethod
    def lambert_solve_from_orbits(body1_orbit_with_epoch,  body2_orbit_with_epoch,  main_attractor=Sun):
        """
        Solves Lambert's problem using two instantiated Orbit objects from poliastro.twobody.Orbit. These orbits must
         have  an associated EPOCH passed to it during initiation: __init__(*args,**kwargs). (EPOCH 2 !< EPOCH 1)
        :param body1_orbit_with_epoch: The EPOCH passed to body1's orbit is the DEPARTURE from the first body.
        :param body2_orbit_with_epoch: The EPOCH passed to body2's orbit is the ARRIVAL to the second body.
        :param main_attractor: By default the attractor is set to the Sun for interplanetary trajectories.
        :return: (namedtuple) Collection of parameters defining the input and outputs of the Lambert problem.
        """
        ss0 = body1_orbit_with_epoch
        ss1 = body2_orbit_with_epoch

        time_of_flight = ss1.epoch - ss0.epoch                                      # Input parameter TOF
        (r0, r) = zip(ss0.r, ss1.r)                                                 # Input parameters  (r0, r)
        (v0, v) = iod.lambert(main_attractor.k, ss0.r, ss1.r, time_of_flight)       # Output parameters (v0, v)

        return lambert_parameters(r0=r0, r1=r, v0=v0, v1=v, tof=time_of_flight, attractor=main_attractor,
                                  epoch0=ss0.epoch, epoch1=ss1.epoch, ss0=ss0, ss1=ss1)

    @staticmethod
    def lambert_solve_from_bodies(body1, body2, epoch1, epoch2, main_attractor=Sun):
        """
        Solves Lambert's problem using two body objects from poliastro.bodies and respective epochs in time.Time.
        :param body1: (poliastro.body) DEPARTURE body.
        :param body2: (poliastro.body) ARRIVAL/TARGET body.
        :param epoch1: (time.Time) DEPARTURE time.
        :param epoch2: (time.Time) ARRIVAL time.
        :param main_attractor: By default the attractor is set to the Sun for interplanetary trajectories.
        :return: (namedtuple) Collection of parameters defining the input and outputs of the Lambert problem.
        """

        ss0 = Orbit.from_body_ephem(body1, epoch1)
        ss1 = Orbit.from_body_ephem(body2, epoch2)

        time_of_flight = epoch2 - epoch1                                               # Input parameter TOF
        (v0, v), = iod.lambert(main_attractor.k, ss0.r, ss1.r, time_of_flight)         # Output parameters (v0, v)

        return lambert_parameters(r0=ss0.r, r1=ss1.r, v0=v0, v1=v, tof=time_of_flight, attractor=main_attractor,
                                  epoch0=ss0.epoch, epoch1=ss1.epoch, ss0=ss0, ss1=ss1)

    def lambert_hohman_solution(self, *args, **kwargs):
        solution = self.lambert_solve_from_bodies(*args, **kwargs)

        yield

    def __init__(self):
        self.x = 0


if __name__ == '__main__':
    test = TrajectoryTool()

    date_launch = time.Time('2018-05-30 15:02', scale='utc')
    date_arrival = time.Time('2021-05-30 05:17', scale='utc')

    lambert_solution = test.lambert_solve_from_bodies(body1=Earth, body2=Jupiter, epoch1=date_launch,
                                                      epoch2=date_arrival)

    test_orbit = Orbit.from_vectors(lambert_solution.attractor, lambert_solution.r0, lambert_solution.v0,
                                    lambert_solution.epoch0)

    op = OrbitPlotter()
    op.plot(lambert_solution.ss1)
    op.plot(lambert_solution.ss0)
    op.plot(test_orbit)

    plt.show()
