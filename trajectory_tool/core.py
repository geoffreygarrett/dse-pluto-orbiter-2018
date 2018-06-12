
import datetime
from collections import namedtuple
from copy import deepcopy
# import networkx
import plotly.plotly as py
import plotly.graph_objs as go

import plotly.plotly as py
import plotly.figure_factory as FF
import plotly.graph_objs as go

import numpy as np
from scipy.spatial import Delaunay

key = 'kHLfFnsPiyxyAfWXgLN6'
user = 'Jones1311'

import plotly
plotly.tools.set_credentials_file(username=user, api_key=key)
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy import time
from astropy.coordinates import solar_system_ephemeris, get_body_barycentric_posvel
from poliastro import iod
from poliastro.bodies import Sun, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto
from poliastro.plotting import OrbitPlotter, OrbitPlotter3D
from poliastro.twobody import Orbit
from poliastro.util import time_range
from trajectory_tool.helper import body_d_domain

lambert_parameters = namedtuple('lambert_parameters', 'r0 r1 v0 v1 tof attractor epoch0 epoch1 ss0 ss1, sst')

plt.style.use("seaborn")
solar_system_ephemeris.set("jpl")

color_earth0 = '#3d4cd5'
color_earthf = '#525fd5'
color_mars0 = '#ec3941'
color_marsf = '#ec1f28'
color_sun = '#ffcc00'
color_orbit = '#888888'
color_trans = '#444444'

color_legs = ['#00FF66', '#00FFFF', '#FF00FF']
color_legs = color_legs + [color_legs[0]]

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

ins = {
    'sma': 2175  # km
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

    # @staticmethod
    def hohmanize_lambert(self, body0, body1, epoch0, main_attractor=Sun):
        flight_times = np.arange(0.25, 20, 0.25)
        lambert_solutions = []
        v_planet_list = []
        epoch0time = time.Time(epoch0, scale='tdb')
        ss0 = Orbit.from_body_ephem(body0, epoch0time)

        x = []
        y = []
        z = []
        temp = np.arange(0.25, 20, 0.25)
        epochee = [epoch0 + time.TimeDelta(365*u.d*j) for j in temp]

        for idx, epoch in enumerate(epochee):
            epoch0 = epoch
            x.append(temp[idx])

            for tof in flight_times:
                epoch1 = epoch0 + time.TimeDelta(tof*u.d*365, scale='tdb')
                _tof = epoch1-epoch0
                ss1 = Orbit.from_body_ephem(body1, epoch1)
                (v0, v), = iod.lambert(main_attractor.k, ss0.r, ss1.r, _tof)
                lambert_solutions.append((v0, v))
                print(v0, v)
                v_planet_list.append(ss1.state.v.to(u.km/u.s))
                y.append(tof)

            for lambert, v_planet in zip(lambert_solutions, v_planet_list):
                print(np.linalg.norm(v_planet.to(u.km/u.s) - lambert[1].to(u.km/u.s)))
                z.append(np.linalg.norm(v_planet.to(u.km/u.s) - lambert[1].to(u.km/u.s)))



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
        (v0, v), = iod.lambert(main_attractor.k, ss0.r, ss1.r, time_of_flight)  # Output parameters (v0, v)
        sst = Orbit.from_vectors(main_attractor, ss0.r, v0, epoch=body1_orbit_with_epoch.epoch)

        return lambert_parameters(r0=ss0.r, r1=ss1.r, v0=v0, v1=v, tof=time_of_flight, attractor=main_attractor,
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
        self.N = 350

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

    def refined_gravity_assist(self, v_s_i_initial, v_s_f_initial, v_planet_initial, body_assisting, body_next, epoch_assist,
                               epoch_next_body, epoch_previous_body=None, previous_body=None, mode='fast',
                               verification=False):

        _return = []
        v_s_f = v_s_f_initial
        v_s_i = v_s_i_initial
        delta_t_assist = 0
        v_planet_entry = v_planet_initial
        v_planet_exit = v_planet_initial
        delta_delta_t = 200
        iterations = 0

        # Parse body name in lower string form for trajectory plotter.
        body_assist_str = body_assisting.__str__().split(' ')[0].lower()
        body_next_str = body_next.__str__().split(' ')[0].lower()

        # Sphere of influence for body.
        r_soi = body_d_domain[body_assist_str]['upper'] * (u.km)
        r_atm = body_d_domain[body_assist_str]['lower'] * (u.km)

        # Lambert iteration TO body.
        if mode is 'full':
            assert (previous_body and epoch_previous_body), "Previous Lambert iteration requires these arguments."
            initial_epoch_entry = epoch_assist

        while delta_delta_t >= 100:

            # Calculate the s/c relative velocity to planet.
            v_s_p_i = (v_s_i - v_planet_entry).to(u.km / u.s)
            v_s_p_f = (v_s_f - v_planet_exit).to(u.km / u.s)

            # Calculate incoming and outgoing angles relative to planet velocity vector.
            theta_i = self.angle_between_cc(v_planet_entry, v_s_p_i)
            theta_f = self.angle_between_cc(v_planet_exit, v_s_p_f)

            # The change from incoming to outgoing.
            delta = theta_f - theta_i

            # Classical orbital parameters for hyperbolic trajectory.
            e = abs(1 / np.sin(delta / 2))
            a = -body_assisting.k.to(u.km ** 3 / u.s ** 2) / (np.linalg.norm(deepcopy(v_s_p_i)))**2 / \
                (u.km ** 2 / u.s ** 2)
            b = np.sqrt(np.square(a) * (e ** 2 - 1))
            r_p = -a * (e - 1)

            # Transform velocity back to heliocentric reference frame.
            v_out = np.linalg.norm(v_s_p_i) * self.unit_vector(v_s_p_f) + v_planet_exit

            # Definition of plane of gravity assist
            n = np.cross(v_s_p_i, v_s_p_f)
            b_unit_v_enter = self.unit_vector(np.cross(v_s_p_i, n))

            # Check which direction the b_vector should go in.
            if (self.angle_between_cc(b_unit_v_enter, v_s_p_f) > np.pi / 2 and self.angle_between_cc(v_s_p_f,
                                                                                        b_unit_v_enter) > np.pi / 2):
                b_vec_enter = b_unit_v_enter * b.value * (u.km)
            else:
                b_vec_enter = -b_unit_v_enter * b.value * (u.km)

            # Calculating the position on the SOI for entrance.
            x_mag_ent = np.sqrt(r_soi.value ** 2 - b.value ** 2)
            x_vec_ent = self.unit_vector(np.negative(v_s_p_i)) * x_mag_ent
            r_ent = x_vec_ent.value + b_vec_enter.value

            b_unit_v_exit = self.unit_vector(np.cross(v_s_p_f, n))

            # Check which direction the b_vector should go in.
            if (self.angle_between_cc(b_unit_v_exit, np.negative(v_s_p_i)) > np.pi / 2 and self.angle_between_cc(np.negative(v_s_p_i),
                                                                                        b_unit_v_exit) > np.pi / 2):
                b_vec_exit = b_unit_v_exit * b.value * (u.km)
            else:
                b_vec_exit = -b_unit_v_exit * b.value * (u.km)

            # Calculating the position on the SOI for exit
            x_mag_exit = np.sqrt(r_soi.value ** 2 - b.value ** 2)
            x_vec_exit = self.unit_vector(v_s_p_f) * x_mag_exit
            r_ext = x_vec_exit.value + b_vec_exit.value

            # Velocity at perigee
            v_p = np.sqrt(body_assisting.k / r_p) * (1 / e)

            # Calculating the hyperbolic anomaly
            def H(theta):
                return np.arccosh(abs((e + np.cos(theta)) / (1 + e * np.cos(theta))))

            # Calculating the time since passage of perigee
            def t_p(H):
                return np.sqrt((-a) ** 3 / (body_assisting.k)) * (e * np.sinh(H) - H)

            # Calculating the angle of the exit.
            theta_rsoi_exit = (2 * np.pi - self.angle_between(r_ent, r_ext)) / 2

            # Time taken for entire manoeuvre.
            delta_t_previous = delta_t_assist

            half_delta_t_assist = t_p(H(theta_rsoi_exit)).to(u.s)
            delta_t_assist = 2 * half_delta_t_assist

            delta_delta_t = abs(delta_t_assist - delta_t_previous).to(u.s).value

            # Epoch during exit
            epoch_exit = epoch_assist + time.TimeDelta(half_delta_t_assist)


            # Epoch during entry
            epoch_entry = epoch_assist - time.TimeDelta(half_delta_t_assist)

            ss_assist = Orbit.from_body_ephem(body_assisting, epoch_entry)
            try:
                # New Lambert solution for Previous Leg
                (v_launch, v_s_i_new), = iod.lambert(Sun.k, Orbit.from_body_ephem(previous_body, epoch_previous_body).r,
                                                     ss_assist.r + r_ent * (u.km), epoch_entry-epoch_previous_body)
                v_planet_entry = ss_assist.state.v
                v_s_i = v_s_i_new

            except AssertionError:
                raise ValueError("The gravity assist is not possible...\n"
                                 "Divergent time iteration.")
                # print('e_prev: ', epoch_previous_body)
                # print('e_entr: ', epoch_entry)
                # print('e_assi: ', epoch_assist)

            # print('Clean ------------------------')
            # print('e_prev: ', epoch_previous_body)
            # print('e_entr: ', epoch_entry)
            # print('e_assi: ', epoch_assist)
            # print('e_exit: ', epoch_exit)
            # print('\n')

            try:
                 ss_assist_exit = Orbit.from_body_ephem(body_assisting, epoch_exit)
                 (v_s_f_new, v0_arrival), = iod.lambert(Sun.k, ss_assist_exit.r + r_ext*(u.km),
                                                               Orbit.from_body_ephem(body_next, epoch_next_body).r,
                                                               epoch_next_body-epoch_exit)
                 v_planet_exit = Orbit.from_body_ephem(body_assisting, epoch_exit).state.v
                 v_s_f = v_s_f_new

            except AssertionError:
                raise ValueError("The gravity assist is not possible...\n"
                                 "Divergent time iteration.")

            if iterations >= 200:
                raise ValueError("The gravity assist is not possible...\n"
                                 "Divergent time iteration.")

            iterations += 1

        # Check to see if closest approach is below set body limit.
        if (r_p > r_atm) and (r_p <= 0.7*r_soi):
            dv_extra_in = np.linalg.norm(v_s_i - v_s_i_initial)
            dv_extra_out = np.linalg.norm(v_s_f - v_out)
            _return += [dv_extra_out, dv_extra_in]

        else:
            raise ValueError("The gravity assist is not possible...\n"
                             "r_p: {:0.2f} < r_min: {:0.2f}".format(r_p, r_atm))
        return _return

    def optimise_gravity_assist(self, v_s_i, v_s_f, v_p, body, epoch, mode='fast', verification=False):
        _return = []
        possible_modes = set(['fast', 'plot2D', 'plot3D', 'verification'])

        if set([mode]) < possible_modes:
            pass
        else:
            raise ValueError("Mode not recognised for gravity assist.")

        # Parse body name in lower string form for trajectory plotter.
        body_str = body.__str__().split(' ')[0].lower()

        # Sphere of influence for body.
        r_soi = body_d_domain[body_str]['upper'] * (u.km)
        r_atm = body_d_domain[body_str]['lower'] * (u.km)

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
        b = np.sqrt(np.square(a)*(e**2 - 1))
        # d = b/np.sin(theta_i)

        # Closest approach distance.
        r_p = -a * (e - 1)

        # Transform velocity back to heliocentric reference frame.
        v_out = np.linalg.norm(v_s_p_i) * self.unit_vector(v_s_p_f) + v_p

        # Check to see if closest approach is below set body limit.
        if r_p > r_atm:
            dv_extra = np.linalg.norm(v_s_f - v_out)
            _return.append(dv_extra)

        else:
            raise ValueError("The gravity assist is not possible...\n"
                             "r_p: {:0.2f} < r_min: {:0.2f}".format(r_p, r_atm))
            # dv_extra = None

        # Definition of plane of gravity assist
        n = np.cross(v_s_p_i, v_s_p_f)

        b_unit_v1 = self.unit_vector(np.cross(v_s_p_f, n))

        # Check which direction the b_vector should go in.
        if (self.angle_between_cc(b_unit_v1, v_s_p_f) > np.pi/2 and self.angle_between_cc(v_s_p_f, b_unit_v1) > np.pi/2):
            b_vec1 = b_unit_v1 * b.value * (u.km)
        else:
            b_vec1 = -b_unit_v1 * b.value * (u.km)

        x_mag1 = np.sqrt(r_soi.value**2 - b.value**2)
        x_vec1 = self.unit_vector(np.negative(v_s_p_f)) *  x_mag1
        r_ext = x_vec1.value+b_vec1.value

        b_unit_v = self.unit_vector(np.cross(v_s_p_i, n))

        # Check which direction the b_vector should go in.
        if (self.angle_between_cc(b_unit_v, v_s_p_f) > np.pi/2 and self.angle_between_cc(v_s_p_f, b_unit_v) > np.pi/2):
            b_vec = b_unit_v * b.value * (u.km)
        else:
            b_vec = -b_unit_v * b.value * (u.km)

        x_mag = np.sqrt(r_soi.value**2 - b.value**2)
        x_vec = self.unit_vector(np.negative(v_s_p_i)) *  x_mag
        r_ent = x_vec.value+b_vec.value

        v_per = np.sqrt(body.k/r_p)*(1/e)
        # t_p = np.sqrt((-a)**3/(body.k))*(e*np.sinh(H_inf)-H_inf)

        def H(theta):
            return np.arccosh(abs((e + np.cos(theta)) / (1 + e * np.cos(theta))))

        def t_p(H):
            return np.sqrt((-a) ** 3 / (body.k)) * (e * np.sinh(H) - H)

        theta_rsoi_exit = (2*np.pi - self.angle_between(r_ent, r_ext))/2
        delta_t_assist = 2*t_p(H(theta_rsoi_exit)).to(u.s)

        if verification:
            _return = _return + [r_ent] + [r_p] + [e] + [a]
            print('r_entry = ', r_ent)
            print('v_entry = ', v_s_p_i)
            print('r_exit  =', r_ext)
            print('v_exit  =', v_s_p_f)
            print('r_p =', r_p)

        if set([mode]) < set(['plot2D', 'plot3D']):

            if mode is 'plot2D':
                ss = OrbitPlotter()
            elif mode is 'plot3D':
                ss = OrbitPlotter3D()

            ss_soi = Orbit.circular(body, alt=r_soi, epoch=epoch)
            ss.plot(ss_soi, label=str(body) + ' SOI', color='red')
            ss_hyp2 = Orbit.from_classical(body, a, e * (u.km / u.km), inc=0 * u.rad, raan=0 * u.rad, argp=0 * u.rad,
                                          nu=-0.5 * u.rad, epoch=epoch)

            ss_hyp = Orbit.from_vectors(body, r=r_ent*(u.km), v=v_s_p_i, epoch=epoch)
            tv = time_range(start=epoch, periods=150, end=epoch + time.TimeDelta(100 * u.day))
            ss.plot(ss_hyp2, label='test', color='0.8')
            ss.plot_trajectory(ss_hyp.sample(tv)[-1], label=str(body), color='green')

            if mode is 'plot2D':
                plt.show()
            elif mode is 'plot3D':
                ss.set_view(30 * u.deg, 260 * u.deg, distance=3 * u.km)
                ss.show(title="EJP Example")
                ss.savefig("EJPExample.png", title="MSL Mission: from Earth to Mars")



        return _return

    def process_itinerary(self, _raw_itinerary, _body_list, _mode='fast', _grav_ass=False, verbose=False):
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

        if verbose:
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
        if verbose:
            print('Solving Lambert multi-leg problem...'.ljust(40) + ' ID: {}\n'.format(_raw_itinerary['id']))
        for i in range(len(_raw_itinerary['durations'])):
            _itinerary_data[i + 1]['l'] = self.lambert_solve_from_bodies(_itinerary_data[i]['b'],
                                                                         _itinerary_data[i + 1]['b'],
                                                                         _itinerary_data[i]['d'],
                                                                         _itinerary_data[i + 1]['d'])

        if _mode is 'delta_v' or 'plot2D' or 'full' or 'plot3D' or 'dv':
            # DEPARTURE BODY DATA --------------------------------------------------------------------------------------
            _itinerary_data[0]['v']['p'] = _itinerary_data[1]['l'].ss0.state.v.to(u.km / u.s)
            _itinerary_data[0]['v']['d'] = _itinerary_data[1]['l'].v0.to(u.km / u.s)

            v_inf = np.linalg.norm((_itinerary_data[0]['v']['d'] - _itinerary_data[0]['v']['p']).to(
                    u.km / u.s).value)

            body = _itinerary_data[0]['b']

            r_0 = epo['alt'] + body.R.to(u.km).value
            v_0 = np.sqrt(np.square(v_inf*(u.km/u.s))+2*body.k.to(u.km**3/u.s**2)/(r_0*(u.km)))
            _itinerary_data[0]['dv'] = (v_0 - np.sqrt(body.k.to(u.km**3/u.s**2)/(r_0*(u.km)))).value

            # INTERMEDIATE BODIES --------------------------------------------------------------------------------------
            for i in range(len(_raw_itinerary['durations']) - 1):
                #   # ARRIVAL, PLANET AND DEPARTURE VELOCITY OF BODY i (1)
                _itinerary_data[i + 1]['v']['a'] = _itinerary_data[i + 1]['l'].v1.to(u.km / u.s)
                _itinerary_data[i + 1]['v']['p'] = _itinerary_data[i + 1]['l'].ss1.state.v.to(u.km / u.s)
                _itinerary_data[i + 1]['v']['d'] = _itinerary_data[i + 2]['l'].v0.to(u.km / u.s)

                if verbose:
                    print('Optimising gravity assist...'.ljust(40) + ' ID: {}\n'.format(_raw_itinerary['id']))
                #   # DELTA V (NO GRAVITY ASSIST) PASSING BODY i (1)

                if _mode is "dv":
                    _itinerary_data[i + 1]['dv'] = sum(
                        self.refined_gravity_assist(v_s_i_initial=_itinerary_data[i + 1]['v']['a'],
                                                    v_s_f_initial=_itinerary_data[i + 1]['v']['d'],
                                                    v_planet_initial=_itinerary_data[i + 1]['v']['p'],
                                                    body_assisting=_itinerary_data[i + 1]['b'],
                                                    body_next = _itinerary_data[i+2]['b'],
                                                    epoch_assist=_itinerary_data[i + 1]['d'],
                                                    epoch_next_body=_itinerary_data[i+2]['d'],
                                                    epoch_previous_body=_itinerary_data[i]['d'],
                                                    previous_body=_itinerary_data[i]['b'],
                                                    mode='full')
                )

                if _mode is "delta_v":
                    _itinerary_data[i+1]['dv'] = \
                    self.refined_gravity_assist(v_s_i_initial=_itinerary_data[i + 1]['v']['a'],
                                                v_s_f_initial=_itinerary_data[i + 1]['v']['d'],
                                                v_planet_initial=_itinerary_data[i + 1]['v']['p'],
                                                body_assisting=_itinerary_data[i + 1]['b'],
                                                body_next=_itinerary_data[i + 2]['b'],
                                                epoch_assist=_itinerary_data[i + 1]['d'],
                                                epoch_next_body=_itinerary_data[i + 2]['d'],
                                                epoch_previous_body=_itinerary_data[i]['d'],
                                                previous_body=_itinerary_data[i]['b'],
                                                mode='full')

            # def refined_gravity_assist(self, v_s_i, v_s_f_initial, v_planet_initial, body_assisting, body_next,
            #                            epoch_entry, epoch_next_body, mode='fast', verification=False):
            # ARRIVAL BODY----------------------------------------------------------------------------------------------
            #   # ARRIVAL, PLANET  VELOCITY OF TARGET BODY i (N)
            _itinerary_data[len(_raw_itinerary['durations'])]['v']['a'] = \
                _itinerary_data[len(_raw_itinerary['durations'])]['l'].v1

            _itinerary_data[len(_raw_itinerary['durations'])]['v']['p'] = \
                _itinerary_data[len(_raw_itinerary['durations'])]['l'].ss1.state.v.to(u.km / u.s)

            #   # DELTA V (NO GRAVITY ASSIST) PASSING BODY i (N)
            idx_arrival = len(_raw_itinerary['durations'])

            _v_0 = _itinerary_data[idx_arrival]['v']['p'].to(u.km/u.s).value- _itinerary_data[idx_arrival]['v']['a'].to(u.km/u.s).value
            _v_orbit = np.sqrt(_itinerary_data[idx_arrival]['b'].k.to(u.km**3/u.s**2).value/ins['sma'])

            _itinerary_data[idx_arrival]['dv'] = \
                np.linalg.norm( _v_0 - _v_orbit)


                # np.linalg.norm(((_itinerary_data[len(_raw_itinerary['durations'])]['v']['p'] -
                #                  _itinerary_data[len(_raw_itinerary['durations'])]['v']['a'])).to(u.km / u.s))

            if verbose:
                print('Delta-v result: {:0.2f} km/s'.format(
                    sum([_itinerary_data[i]['dv'] for i in range(len(_itinerary_data))])).ljust(40)
                      + ' ID: {}\n'.format(_raw_itinerary['id']))

        if (_mode is 'plot' or 'full'):
            # TRAJECTORIES OF LEGS -------------------------------------------------------------------------------------
            for i in range(len(_raw_itinerary['durations'])):
                _itinerary_data[i + 1]['t'] = Orbit.from_vectors(_itinerary_data[i + 1]['l'].attractor,
                                                                 _itinerary_data[i + 1]['l'].r0,
                                                                 _itinerary_data[i + 1]['l'].v0)

        if set([_mode]) < set(['plot2D', 'plot3D']):
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

            if _mode is 'plot3D':
                frame = OrbitPlotter3D()
            elif _mode is 'plot2D':
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
                    color=color_legs[i],
                    )

            #TEMP
            # frame.plot_trajectory(Orbit.from_vectors(Jupiter, np.array([-47468291.44350722, 7154251.0880544, 4337971.73842018])*(u.km)+_itinerary_data[2]['l'].r0, np.array([-6.9512146, 0.87033231, 0.31982218])*(u.km/u.s)+_itinerary_data[2]['v']['p'], _itinerary_data[2]['l'].epoch0).sample(_itinerary_data[1]['tv'])[-1])

            frame._redraw_attractor(0.25 * 10 ** (8) * u.km)
            print('Displaying plot!'.ljust(40) + ' ID: {}\n'.format(_raw_itinerary['id']))

            if _mode is 'plot2D':
                plt.legend()
                plt.show()

            else:
                frame.set_view(30 * u.deg, 260 * u.deg, distance=3 * u.km)
                # frame.show(title="EJP Example")
                frame.savefig("EJPExample.png", title="EJP Example trajectory sequence")

        if verbose:
            print('Complete!'.ljust(40) + ' ID: {}'.format(_raw_itinerary['id']))
            print('-' * 40 + '-' * len(' ID: {}\n'.format(_raw_itinerary['id'])))
        return _itinerary_data



if __name__ == '__main__':

    ####################################################################################################################
    test = True
    _test = TrajectoryTool()
    ####################################################################################################################
    if test:
        # for i in range(100):
        # TEST EJP -------------------------------------------------------------------------------------------------
        __raw_itinerary1 = ['earth', 'jupiter', 'pluto']
        __raw_itinerary2 = {'id': 1,
                            'launch_date': datetime.datetime(2027, 12, 1, 0, 0),
                            'durations': [2.115, 22.852]
                            }
        # ----------------------------------------------------------------------------------------------------------
        processed = _test.process_itinerary(__raw_itinerary2, __raw_itinerary1, _mode='fast', _grav_ass=True)

        total = [24-float(i) for i in np.linspace(0, 9, 100)]
        # print(total)
        # def optimise_gravity_assist(self, v_s_i, v_s_f, v_p, body, epoch, plot=False, verification=False):

        # res = _test.optimise_gravity_assist(v_s_i=processed[1]['v']['a'],
        #                                       v_s_f=processed[1]['v']['d'],
        #                                       v_p  =processed[1]['v']['p'],
        #                                       body =processed[1]['b'],
        #                                       epoch=processed[1]['d'],
        #                                       mode='plot3D',
        #                                       verification=True
        #                                       )

        # def refined_gravity_assist(self, v_s_i, v_s_f_initial, v_planet, body_assisting, body_next, epoch_entry,
        #                            epoch_next_body, mode='fast', verification=False):

        # -----------------------------------------------------------------------------------
        # print(processed[1]['d'])
        # ref = _test.refined_gravity_assist(v_s_i=processed[1]['v']['a'],
        #                                    v_s_f_initial=processed[1]['v']['d'],
        #                                    v_planet_initial=processed[1]['v']['p'],
        #                                    body_assisting=processed[1]['b'],
        #                                    epoch_entry=processed[1]['d'],
        #                                    body_next= processed[2]['b'],
        #                                    epoch_next_body=processed[2]['d'],
        #                                    )
        #
        # _test.hohmanize_lambert(processed[1]['b'],
        #                         processed[2]['b'],
        #                         processed[1]['d'])
        #
        # print([processed[i]['dv'] for i in range(len(processed))])
        # print(processed[1]['d'])
        # ------------------------------------------------------------------------------------
        # print(processed[1]['v']['p'] - processed[1]['v']['a'])

        # print(res)

    ####################################################################################################################


    ####################################################################################################################
