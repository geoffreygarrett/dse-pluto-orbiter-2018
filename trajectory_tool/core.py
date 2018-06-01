# IMPORTS -------------------------------------------------------------------------------------------------------------
import matplotlib.pyplot as plt
from matplotlib import cm
from poliastro.frames import HeliocentricEclipticJ2000
from poliastro.plotting import OrbitPlotter3D, OrbitPlotter2D, OrbitPlotter
from poliastro.bodies import Sun, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto, Body
from poliastro.twobody import Orbit
from poliastro import iod
from astropy.coordinates import solar_system_ephemeris, get_body_barycentric_posvel
import astropy.units as u
from astropy import time
from poliastro.util import time_range
from pprint import pprint

import numpy as np
import datetime
from collections import namedtuple

lambert_parameters = namedtuple('lambert_parameters', 'r0 r1 v0 v1 tof attractor epoch0 epoch1 ss0 ss1, sst')
# r0 = lambert_parameters.r0

plt.style.use("seaborn")
solar_system_ephemeris.set("jpl")

color_earth0 = '#3d4cd5'
color_earthf = '#525fd5'
color_mars0 = '#ec3941'
color_marsf = '#ec1f28'
color_sun = '#ffcc00'
color_orbit = '#888888'
color_trans = '#444444'

colors = ['#ec3941', '#ec1f28']


def color_gen():
    while True:
        for color in colors:
            return color


body_list = {
    'mercury': Mercury,
    'venus': Venus,
    'earth': Earth,
    'mars': Mars,
    'jupiter': Jupiter,
    'saturn': Saturn,
    'uranus': Uranus,
    'neptune': Neptune,
    'pluto': Pluto
}


# body_color = {}

# CORE CLASS -----------------------------------------------------------------------------------------------------------
class TrajectoryTool(object):
    """

    """

    @staticmethod
    def hohmanize_lambert(lambert_function):
        pass

    @staticmethod
    def _lambert_solve_from_positions(r0, r1, epoch0, epoch1, main_attractor):
        time_of_flight = epoch1 - epoch0
        (v0, v1), = iod.lambert(main_attractor.K, r0, r1, time_of_flight)
        return v0, v1, time_of_flight

    @staticmethod
    def _lambert_solve_from_orbits(r0, r1, body0_orbit_with_epoch, body1_orbit_with_epoch, main_attractor):
        time_of_flight = body1_orbit_with_epoch.epoch - body0_orbit_with_epoch.epoch
        (v0, v1), = iod.lambert(main_attractor.K, r0, r1, time_of_flight)
        return v0, v1, time_of_flight

    # @staticmethod
    # def lambert_solve(method='from_positions',r1=None, r0=None, epoch0=None, epoch1=None,  body0_orbit_with_epoch=None,
    #                   body1_orbit_with_epoch=None, body0=None, body1=None, main_attractor=Sun):
    #
    #     if method is 'from_positions':
    #         required_args = (r1 and r0 and epoch0 and epoch1)
    #         assert required_args, 'These arguments are required for the - from_positions - method.'
    #
    #
    #     elif method is 'from_orbits':
    #         required_args = (body0_orbit_with_epoch and body1_orbit_with_epoch)
    #         assert required_args, 'These arguments are required for the - from_orbits - method.'
    #
    #     elif method is 'from_bodies':
    #         required_args = (body0 and body1 and epoch0 and epoch1)
    #         assert required_args, 'These arguments are required for the - from_bodies - method.'
    #
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
        (v0, v), = iod.lambert(main_attractor.K, r0, r, tof)

        return lambert_parameters(r0=r0, r1=r, v0=v0, v1=v, tof=tof, attractor=main_attractor, epoch0=epoch0,
                                  epoch1=epoch1, ss0=None, ss1=None, sst=None)

    @staticmethod
    def lambert_solve_from_orbits(body1_orbit_with_epoch, body2_orbit_with_epoch, main_attractor=Sun):
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

        time_of_flight = ss1.epoch - ss0.epoch  # Input parameter TOF
        (r0, r) = zip(ss0.r, ss1.r)  # Input parameters  (r0, r)
        (v0, v) = iod.lambert(main_attractor.k, ss0.r, ss1.r, time_of_flight)  # Output parameters (v0, v)
        sst = Orbit.from_vectors(main_attractor, ss0.r, v0, epoch=body1_orbit_with_epoch.epoch)

        return lambert_parameters(r0=r0, r1=r, v0=v0, v1=v, tof=time_of_flight, attractor=main_attractor,
                                  epoch0=ss0.epoch, epoch1=ss1.epoch, ss0=ss0, ss1=ss1, sst=sst)

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

        time_of_flight = epoch2 - epoch1  # Input parameter TOF
        (v0, v), = iod.lambert(main_attractor.k, ss0.r, ss1.r, time_of_flight)  # Output parameters (v0, v)
        sst = Orbit.from_vectors(main_attractor, ss0.r, v0, epoch=epoch1)

        return lambert_parameters(r0=ss0.r, r1=ss1.r, v0=v0, v1=v, tof=time_of_flight, attractor=main_attractor,
                                  epoch0=ss0.epoch, epoch1=ss1.epoch, ss0=ss0, ss1=ss1, sst=sst)

    # def lambert_hohman_solution(self, *args, **kwargs):
    #     solution = self.lambert_solve_from_bodies(*args, **kwargs)
    #     yield

    def __init__(self):
        self.N = 100

    def process_itinerary(self, _raw_itinerary, _body_list, _mode='plot'):
        """

        :param _raw_itinerary:     raw_itinerary = {'id': int,
                                                    'launch_data': datetime.datetime,
                                                    'durations': list}
        :param _body_list: (str)   [body1, body2, .... bodyN]
        :param _fast_mode: (boolean)
        :param _mode: ['fast', 'plot', 'full']
        :return: processed_itinerary:
        """
        # GENERATE PROCESSED ITINERARY STRUCTURE -----------------------------------------------------------------------
        _itinerary_data = [None]*(len(_raw_itinerary['durations'])+1)

        _itinerary_data[0] = {'b': body_list[_body_list[0]],
                                'd': time.Time(_raw_itinerary['launch_date']),
                                's': _body_list[0],
                                'v': {}}

        for i in range(len(_raw_itinerary['durations'])):
            _itinerary_data[i + 1] = {'b': body_list[_body_list[i + 1]],
                                           'd': time.Time(_raw_itinerary['launch_date'] +
                                                          datetime.timedelta(
                                                              days=365 * _raw_itinerary['durations'][i])),
                                           's': _body_list[i + 1],
                                           'v': {}}

        # LAMBERT SOLUTIONS --------------------------------------------------------------------------------------------
        for i in range(len(_raw_itinerary['durations'])):
            _itinerary_data[i + 1]['l'] = self.lambert_solve_from_bodies(_itinerary_data[i]['b'],
                                                                              _itinerary_data[i + 1]['b'],
                                                                              _itinerary_data[i]['d'],
                                                                              _itinerary_data[i + 1]['d'])

            # CAN TEST FOR GRAVITY ASSIST FEASIBILITY HERE (ASK GEOFF HOW)
            if i != 0:
                curr_lambert_soln = _itinerary_data[i + 1]['l']

                prev_lambert_soln = _itinerary_data[i]['l']
                v_s = prev_lambert_soln.v1
                v_p = prev_lambert_soln.ss1.state.v

                v_s_p = v_s - v_p
                v_s_p_mag = np.linalg.norm(v_s_p)


        if _mode is 'plot' or 'full':
            # TRAJECTORIES OF LEGS -------------------------------------------------------------------------------------
            for i in range(len(_raw_itinerary['durations'])):
                _itinerary_data[i + 1]['t'] = Orbit.from_vectors(_itinerary_data[i + 1]['l'].attractor,
                                                                      _itinerary_data[i + 1]['l'].r0,
                                                                      _itinerary_data[i + 1]['l'].v0)

            # DEPARTURE BODY DATA --------------------------------------------------------------------------------------
            _itinerary_data[0]['v'] = _itinerary_data[1]['l'].v0
            _itinerary_data[0]['dv'] = (_itinerary_data[0]['v'] - _itinerary_data[1]['l'].ss0.state.v).to(
                u.km / u.s)

            # INTERMEDIATE BODIES --------------------------------------------------------------------------------------
            for i in range(len(_raw_itinerary['durations']) - 1):
                #   # ARRIVAL, PLANET AND DEPARTURE VELOCITY OF BODY i (1)
                _itinerary_data[i + 1]['v']['a'] = _itinerary_data[i + 1]['l'].v1
                _itinerary_data[i + 1]['v']['p'] = _itinerary_data[i + 1]['l'].ss1.state.v
                _itinerary_data[i + 1]['v']['d'] = _itinerary_data[i + 2]['l'].v0

                #   # DELTA V (NO GRAVITY ASSIST) PASSING BODY i (1)
                _itinerary_data[i + 1]['dv'] = np.linalg.norm(((_itinerary_data[i + 1]['v']['d'] -
                                                                     _itinerary_data[i + 1]['v']['p']) -
                                                                    _itinerary_data[i + 1]['v']['a']).to(
                    u.km / u.s))

            # ARRIVAL BODY----------------------------------------------------------------------------------------------
            #   # ARRIVAL, PLANET  VELOCITY OF TARGET BODY i (N)
            _itinerary_data[len(_raw_itinerary['durations'])]['v']['a'] = \
                _itinerary_data[len(_raw_itinerary['durations'])]['l'].v1

            _itinerary_data[len(_raw_itinerary['durations'])]['v']['p'] = \
                _itinerary_data[len(_raw_itinerary['durations'])]['l'].ss1.state.v

            #   # DELTA V (NO GRAVITY ASSIST) PASSING BODY i (N)
            _itinerary_data[len(_raw_itinerary['durations'])]['dv'] = \
                np.linalg.norm(((_itinerary_data[len(_raw_itinerary['durations'])]['v']['p'] -
                                 _itinerary_data[len(_raw_itinerary['durations'])]['v']['a'])).to(u.km / u.s))

        if _mode is 'plot':
            # EXTRA PROCESS FOR PLOTTING -------------------------------------------------------------------------------
            for i in range(len(_raw_itinerary['durations'])):
                _itinerary_data[i + 1]['tv'] = time_range(start=_itinerary_data[i]['d'],
                                                               end=_itinerary_data[i + 1]['d'],
                                                               periods=self.N)

            # GENERATE PROCESSED ITINERARY STRUCTURE -------------------------------------------------------------------
            _itinerary_plot_data = [None]*(len(_raw_itinerary['durations'])+1)
            for i in range(len(_raw_itinerary['durations']) + 1):
                _itinerary_plot_data[i] = {}

            # GENERATE VELOCITY AND POSITION VECTORS OF BODIES
            _itinerary_plot_data[0]['rr'], _itinerary_plot_data[0]['vv'] = \
                get_body_barycentric_posvel(_itinerary_data[0]['s'], _itinerary_data[1]['tv'])

            for i in range(len(_raw_itinerary['durations'])):
                _itinerary_plot_data[i + 1]['rr'], _itinerary_plot_data[i + 1]['vv'] = \
                    get_body_barycentric_posvel(_itinerary_data[i+1]['s'], _itinerary_data[i+1]['tv'])

            for i in range(len(_raw_itinerary['durations'])):
                _itinerary_plot_data[i+1]['tp'] = Orbit.from_vectors(_itinerary_data[i+1]['l'].attractor,
                                                                            _itinerary_data[i+1]['l'].r0,
                                                                            _itinerary_data[i+1]['l'].v0,
                                                                            _itinerary_data[i+1]['l'].epoch0)

            frame = OrbitPlotter()
            frame.set_attractor(Sun)

            for i in range(len(_raw_itinerary['durations'])):
                frame.plot(_itinerary_data[i+1]['l'].ss0, color='0.8')

            frame.plot(_itinerary_data[len(_raw_itinerary['durations'])]['l'].ss1, color='0.8')

            for i in range(len(_raw_itinerary['durations']) + 1):
                frame.plot_trajectory(_itinerary_plot_data[i]['rr'], label=_itinerary_data[i]['b'],
                                      color=color_earth0)

            for i in range(len(_raw_itinerary['durations'])):
                frame.plot_trajectory(
                    _itinerary_plot_data[i+1]['tp'].sample(_itinerary_data[i+1]['tv'])[-1],
                    label="Leg {}".format(i + 1),
                    color=color_trans)

            plt.legend()

            frame._redraw_attractor(0.25 * 10 ** (8) * u.km)  # FIX SUN SIZE
            plt.show()

        return _itinerary_data


if __name__ == '__main__':
    # TEST EVJP
    _test = TrajectoryTool()
    __raw_itinerary1 = ['earth', 'venus', 'jupiter', 'pluto']
    __raw_itinerary2 = {'id': 4332,
                        'launch_date': datetime.datetime(2024, 1, 1, 0, 0),
                        'durations': [0.2102, 1.0114, 11.7784]
                        }

    # TEST EJP
    # __raw_itinerary1 = ['earth', 'jupiter', 'pluto']
    # __raw_itinerary2 = {'id': 4332,
    #                     'launch_date': datetime.datetime(2026, 1, 1, 0, 0),
    #                     'durations': [2.5114, 19.7784]
    #                     }

    # TEST EJ
    # __raw_itinerary1 = ['earth', 'jupiter']
    # __raw_itinerary2 = {'id': 4332,
    #                     'launch_date': datetime.datetime(2026, 1, 1, 0, 0),
    #                     'durations': [2.5114]
    #                     }

    processed = _test.process_itinerary(__raw_itinerary2, __raw_itinerary1, _mode='fast')
    pprint(processed)

