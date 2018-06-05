import datetime
from collections import namedtuple
from copy import deepcopy

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy import time
from astropy.coordinates import solar_system_ephemeris, get_body_barycentric_posvel
from poliastro import iod
from poliastro.bodies import Sun, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto
from poliastro.plotting import OrbitPlotter
from poliastro.twobody import Orbit
from poliastro.util import time_range
from trajectory_tool.helper import body_d_domain

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

epo = {
    'body':Earth,
    'alt': 180 #km
}


# body_color = {}

# CORE CLASS -----------------------------------------------------------------------------------------------------------
class TrajectoryTool(object):
    """

    """

    @staticmethod
    def rotation_z(theta):
        return np.array([[np.cos(theta), -np.sin(theta), 0],
                         [np.sin(theta), np.cos(theta), 0],
                         [0, 0, 1]])

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

    def __init__(self):
        self.N = 100

    def unit_vector(self, vector):
        """ Returns the unit vector of the vector.  """
        return vector / np.linalg.norm(vector)

    def angle_between(self, v1, v2):
        """ Returns the angle in radians between vectors 'v1' and 'v2'::

                >>> self.angle_between((1, 0, 0), (0, 1, 0))
                1.5707963267948966
                >>> self.angle_between((1, 0, 0), (1, 0, 0))
                0.0
                >>> self.angle_between((1, 0, 0), (-1, 0, 0))
                3.141592653589793
        """
        v1_u = self.unit_vector(v1)
        v2_u = self.unit_vector(v2)
        return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

    def angle_between_cc(self, v1, v2):
        """ Returns the angle dependent on defined direction from v1 to v2.
        :param v1:
        :param v2:
        :param direction: Default is Clockwise.
        :return:
        """
        a_i = self.angle_between(v1, v2)

        # rotation check
        v2_n = np.matmul(self.rotation_z(0.01), deepcopy(v2))
        a_new = self.angle_between(v1, v2_n)
        if a_new > a_i:  # Clockwise measurement
            angle = a_i
        else:
            angle = 2 * np.pi - a_i  # Adjust for clockwise measurement
        return angle

    def lambert_solve_sequence(self, method='from_positions', r1=None, r0=None, epoch0=None, epoch1=None,
                               body0_orbit_with_epoch=None,
                               body1_orbit_with_epoch=None, body0=None, body1=None, main_attractor=Sun):

        if method is 'from_positions':
            required_args = (r1 and r0 and epoch0 and epoch1)
            assert required_args, 'These arguments are required for the - from_positions - method.'
            return self._lambert_solve_from_positions(r0, r1, epoch0, epoch1, main_attractor)

        elif method is 'from_orbits':
            required_args = (body0_orbit_with_epoch and body1_orbit_with_epoch)
            assert required_args, 'These arguments are required for the - from_orbits - method.'
            return self._lambert_solve_from_orbits(body0_orbit_with_epoch, body1_orbit_with_epoch, main_attractor)

        elif method is 'from_bodies':
            required_args = (body0 and body1 and epoch0 and epoch1)
            assert required_args, 'These arguments are required for the - from_bodies - method.'
            ss0 = Orbit.from_body_ephem(body0, epoch0)
            ss1 = Orbit.from_body_ephem(body1, epoch1)
            return self._lambert_solve_from_positions(ss0.r, ss1.r, epoch0, epoch1, main_attractor)

    def optimise_gravity_assist(self, v_s_i, v_s_f, v_p, body, epoch, plot=True):

        # Parse body name in lower string form for trajectory plotter.
        body_str = body.__str__().split(' ')[0].lower()

        # Calculate the s/c relative velocity to planet.
        v_s_p_i = (v_s_i - v_p).to(u.km / u.s)
        v_s_p_f = (v_s_f - v_p).to(u.km / u.s)

        # Calculate incoming and outgoing angles relative to planet velocity vector.
        theta_i = self.angle_between_cc(v_p, v_s_p_i)
        theta_f = self.angle_between_cc(v_p, v_s_p_f)

        # The change from incoming to outgoing.
        delta = theta_f - theta_i

        # Classical orbital parameters for hyperbolic trajectory.
        e = abs(1 / np.sin(delta / 2))
        a = -body.k.to(u.km ** 3 / u.s ** 2) / (np.dot(deepcopy(v_s_p_i), deepcopy(v_s_p_i))) / (u.km ** 2 / u.s ** 2)

        # Other parameters, possibly not required at all.
        # b = np.sqrt(np.square(a)*(e**2 - 1))
        # d = b/np.sin(theta_i)

        # Closest approach distance.
        r_p = -a * (e - 1)

        # TODO: Incorporate inspection for feasible r_p

        # Sphere of influence for body.
        r_soi = body_d_domain[body_str]['upper'] * (u.km)
        r_atm = body_d_domain[body_str]['lower'] * (u.km)

        # Transform velocity back to heliocentric reference frame.
        v_out = np.linalg.norm(v_s_p_i) * self.unit_vector(v_s_p_f) + v_p

        # Check to see if closest approach is below set body limit.
        if r_p > r_atm:
            dv_extra = np.linalg.norm(v_s_f - v_out)

        else:
            raise ValueError("The gravity assist is not possible...\n"
                             "r_p: {:0.2f} < r_min: {:0.2f}".format(r_p, r_atm))
            # dv_extra = None


        # # Determine hypberbolic orbit from vectors for GMAT verification
        b = np.sqrt(np.square(a) * (e ** 2 - 1))
        v_in_u = self.unit_vector(v_s_p_i)
        r_soi = r_soi
        x_mag = np.sqrt(np.square(r_soi) - np.square(b))
        x_vec = x_mag * np.negative(v_in_u)
        b_vec = b * np.negative(np.matmul(self.rotation_z(np.pi / 2), self.unit_vector(v_s_p_i)))
        r_ent = r_soi.value * self.unit_vector(b_vec.value + x_vec.value) * (u.km)



        # print(r_ent)
        if plot:
            ss = OrbitPlotter()
            ss_soi = Orbit.circular(body, alt=r_soi, epoch=epoch)
            ss.plot(ss_soi, label=str(body) + ' SOI', color='red')
            # print(e)
            ss_hyp2 = Orbit.from_classical(body, a, e * (u.km / u.km), inc=0 * u.rad, raan=0 * u.rad, argp=0 * u.rad,
                                          nu=-0.5 * u.rad, epoch=epoch)

            ss_hyp = Orbit.from_vectors(body, r=r_ent, v=v_s_p_i, epoch=epoch)
            # print(ss_hyp.e)
            # assert ss_hyp.e == e
            tv = time_range(start=epoch, periods=150, end=epoch + time.TimeDelta(100 * u.day))
            ss.plot(ss_hyp2, label='test', color='0.8')
            ss.plot_trajectory(ss_hyp.sample(tv)[-1], label=str(body), color='green')

            # val = 1.0 * np.linalg.norm(a)
            # plt.xlim(-val, val)
            # plt.ylim(-val, val)

            plt.show()

        return dv_extra

    def process_itinerary(self, _raw_itinerary, _body_list, _mode='plot', _grav_ass=False):
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
        _itinerary_data = [None] * (len(_raw_itinerary['durations']) + 1)

        _itinerary_data[0] = {'b': body_list[_body_list[0]],
                              'd': time.Time(_raw_itinerary['launch_date'], scale='tdb'),
                              's': _body_list[0],
                              'v': {}}

        print('\n')
        print('-' * 40 + '-' * len(' ID: {}'.format(_raw_itinerary['id'])))
        print('Initializing...'.ljust(40) + ' ID: {}\n'.format(_raw_itinerary['id']))
        for i in range(len(_raw_itinerary['durations'])):
            _itinerary_data[i + 1] = {'b': body_list[_body_list[i + 1]],
                                      'd': time.Time(_raw_itinerary['launch_date'] +
                                                     datetime.timedelta(
                                                         days=365 * sum(_raw_itinerary['durations'][:i + 1])),
                                                     scale='tdb'),
                                      's': _body_list[i + 1],
                                      'v': {}}

        # LAMBERT SOLUTIONS --------------------------------------------------------------------------------------------
        print('Solving Lambert multi-leg problem...'.ljust(40) + ' ID: {}\n'.format(_raw_itinerary['id']))
        for i in range(len(_raw_itinerary['durations'])):
            _itinerary_data[i + 1]['l'] = self.lambert_solve_from_bodies(_itinerary_data[i]['b'],
                                                                         _itinerary_data[i + 1]['b'],
                                                                         _itinerary_data[i]['d'],
                                                                         _itinerary_data[i + 1]['d'])

        if (_mode is 'delta_v' or 'plot' or 'full'):
            # DEPARTURE BODY DATA --------------------------------------------------------------------------------------
            _itinerary_data[0]['v']['p'] = _itinerary_data[1]['l'].ss0.state.v.to(u.km / u.s)
            _itinerary_data[0]['v']['d'] = _itinerary_data[1]['l'].v0.to(u.km / u.s)

            _itinerary_data[0]['dv'] = np.linalg.norm((_itinerary_data[0]['v']['d'] - _itinerary_data[0]['v']['p']).to(
                u.km / u.s).value)


            # INTERMEDIATE BODIES --------------------------------------------------------------------------------------
            for i in range(len(_raw_itinerary['durations']) - 1):
                #   # ARRIVAL, PLANET AND DEPARTURE VELOCITY OF BODY i (1)
                _itinerary_data[i + 1]['v']['a'] = _itinerary_data[i + 1]['l'].v1.to(u.km / u.s)
                _itinerary_data[i + 1]['v']['p'] = _itinerary_data[i + 1]['l'].ss1.state.v.to(u.km / u.s)
                _itinerary_data[i + 1]['v']['d'] = _itinerary_data[i + 2]['l'].v0.to(u.km / u.s)

                print('Optimising gravity assist...'.ljust(40) + ' ID: {}\n'.format(_raw_itinerary['id']))
                #   # DELTA V (NO GRAVITY ASSIST) PASSING BODY i (1)
                _itinerary_data[i + 1]['dv'] = self.optimise_gravity_assist(v_s_i=_itinerary_data[i + 1]['v']['a'],
                                                                            v_s_f=_itinerary_data[i + 1]['v']['d'],
                                                                            v_p=_itinerary_data[i + 1]['v']['p'],
                                                                            body=_itinerary_data[i + 1]['b'],
                                                                            epoch=_itinerary_data[i + 1]['d'],
                                                                            plot=False)

            # ARRIVAL BODY----------------------------------------------------------------------------------------------
            #   # ARRIVAL, PLANET  VELOCITY OF TARGET BODY i (N)
            _itinerary_data[len(_raw_itinerary['durations'])]['v']['a'] = \
                _itinerary_data[len(_raw_itinerary['durations'])]['l'].v1

            _itinerary_data[len(_raw_itinerary['durations'])]['v']['p'] = \
                _itinerary_data[len(_raw_itinerary['durations'])]['l'].ss1.state.v.to(u.km / u.s)

            #   # DELTA V (NO GRAVITY ASSIST) PASSING BODY i (N)
            _itinerary_data[len(_raw_itinerary['durations'])]['dv'] = \
                np.linalg.norm(((_itinerary_data[len(_raw_itinerary['durations'])]['v']['p'] -
                                 _itinerary_data[len(_raw_itinerary['durations'])]['v']['a'])).to(u.km / u.s))


            print('Delta-v result: {:0.2f} km/s'.format(
                sum([_itinerary_data[i]['dv'] for i in range(len(_itinerary_data))])).ljust(40)
                  + ' ID: {}\n'.format(_raw_itinerary['id']))

        if (_mode is 'plot' or 'full'):
            # TRAJECTORIES OF LEGS -------------------------------------------------------------------------------------
            for i in range(len(_raw_itinerary['durations'])):
                _itinerary_data[i + 1]['t'] = Orbit.from_vectors(_itinerary_data[i + 1]['l'].attractor,
                                                                 _itinerary_data[i + 1]['l'].r0,
                                                                 _itinerary_data[i + 1]['l'].v0)

        if _mode is 'plot':
            print('Plotting...'.ljust(40) + ' ID: {}\n'.format(_raw_itinerary['id']))
            # EXTRA PROCESS FOR PLOTTING -------------------------------------------------------------------------------
            for i in range(len(_raw_itinerary['durations'])):
                _itinerary_data[i + 1]['tv'] = time_range(start=_itinerary_data[i]['d'],
                                                          end=_itinerary_data[i + 1]['d'],
                                                          periods=self.N)

            # GENERATE PROCESSED ITINERARY STRUCTURE -------------------------------------------------------------------
            _itinerary_plot_data = [None] * (len(_raw_itinerary['durations']) + 1)
            for i in range(len(_raw_itinerary['durations']) + 1):
                _itinerary_plot_data[i] = {}

            # GENERATE VELOCITY AND POSITION VECTORS OF BODIES
            _itinerary_plot_data[0]['rr'], _itinerary_plot_data[0]['vv'] = \
                get_body_barycentric_posvel(_itinerary_data[0]['s'], _itinerary_data[1]['tv'])

            for i in range(len(_raw_itinerary['durations'])):
                _itinerary_plot_data[i + 1]['rr'], _itinerary_plot_data[i + 1]['vv'] = \
                    get_body_barycentric_posvel(_itinerary_data[i + 1]['s'], _itinerary_data[i + 1]['tv'])

            for i in range(len(_raw_itinerary['durations'])):
                _itinerary_plot_data[i + 1]['tp'] = Orbit.from_vectors(_itinerary_data[i + 1]['l'].attractor,
                                                                       _itinerary_data[i + 1]['l'].r0,
                                                                       _itinerary_data[i + 1]['l'].v0,
                                                                       _itinerary_data[i + 1]['l'].epoch0)

            frame = OrbitPlotter()
            frame.set_attractor(Sun)

            for i in range(len(_raw_itinerary['durations'])):
                frame.plot(_itinerary_data[i + 1]['l'].ss0, color='0.8')

            frame.plot(_itinerary_data[len(_raw_itinerary['durations'])]['l'].ss1, color='0.8')

            for i in range(len(_raw_itinerary['durations']) + 1):
                frame.plot_trajectory(_itinerary_plot_data[i]['rr'], label=_itinerary_data[i]['b'],
                                      color=color_earth0)

            for i in range(len(_raw_itinerary['durations'])):
                frame.plot_trajectory(
                    _itinerary_plot_data[i + 1]['tp'].sample(_itinerary_data[i + 1]['tv'])[-1],
                    label="Leg {}".format(i + 1),
                    color=color_trans)

            plt.legend()
            frame._redraw_attractor(0.25 * 10 ** (8) * u.km)  # FIX SUN SIZE
            print('Displaying plot!'.ljust(40) + ' ID: {}\n'.format(_raw_itinerary['id']))
            plt.show()

        print('Complete!'.ljust(40) + ' ID: {}'.format(_raw_itinerary['id']))
        print('-' * 40 + '-' * len(' ID: {}\n'.format(_raw_itinerary['id'])))
        return _itinerary_data


if __name__ == '__main__':

    ####################################################################################################################
    test = True
    _test = TrajectoryTool()
    ####################################################################################################################
    if test:
        for i in range(100):
            # TEST EJP -------------------------------------------------------------------------------------------------
            __raw_itinerary1 = ['earth', 'jupiter', 'pluto']
            __raw_itinerary2 = {'id': i,
                                'launch_date': datetime.datetime(2028, 12, 20, 0, 0),
                                'durations': [1.67397, 22.96712]
                                }
            # ----------------------------------------------------------------------------------------------------------

            processed = _test.process_itinerary(__raw_itinerary2, __raw_itinerary1, _mode='delta_v', _grav_ass=True)
<<<<<<< Updated upstream


    ####################################################################################################################
