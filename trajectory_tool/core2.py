########################################################################################################################
# IMPORTS                                                                                                              #
########################################################################################################################

# GENERAL & PROJECT-MADE
import datetime
import numpy as np
import math
import matplotlib.pyplot as plt
from trajectory_tool.helper import *
from trajectory_tool.plotting import *
from copy import deepcopy, copy
from scipy import optimize

#TEMP
import plotly
import plotly.graph_objs as go

# ASTROPHYSICS & ORBITAL MECHANICS
import astropy.units as u
from astropy import time
from astropy.coordinates import solar_system_ephemeris, get_body_barycentric_posvel
from poliastro import iod
from poliastro.bodies import Sun, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto
from poliastro.plotting import OrbitPlotter, OrbitPlotter3D
from poliastro.twobody import Orbit
from poliastro.util import time_range
from poliastro.maneuver import Maneuver

########################################################################################################################
# CONFIGURATION                                                                                                        #
########################################################################################################################

# GLOBAL CONFIGURATION
plt.style.use("seaborn")
solar_system_ephemeris.set("jpl")
color_legs = ['#00FF66', '#00FFFF', '#FF00FF']
color_legs2 = ['#00FF66', '#33c7ff', '#FF00FF', '#00FF66']
color_legs = color_legs + [color_legs[0]]
colors = ['#ec3941', '#ec1f28']

# BODY OBJECT DICTIONARY FROM NAME
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

########################################################################################################################
# MISSION SPECIFIC PARAMETERS                                                                                          #
########################################################################################################################

# EPO INJECTION PARAMETERS
epo = {
    'body': Earth,
    'alt': 180  # km
}

# PLUTO INSERTION PARAMETERS
ins = {
    'sma': 2175  # km
}

########################################################################################################################
# TRAJECTORY TOOL CLASS CONSTRUCTOR                                                                                    #
########################################################################################################################


class TrajectoryTool(object):
    """

    """
    @staticmethod
    def _body_string_lower(body_object):
        return body_object.__str__().split(' ')[0].lower()

    @staticmethod
    def rotation_matrix(axis, theta):
        """
        Return the rotation matrix associated with counterclockwise rotation about
        the given axis by theta radians.
        """
        axis = np.asarray(axis)
        axis = axis / math.sqrt(np.dot(axis, axis))
        a = math.cos(theta / 2.0)
        b, c, d = -axis * math.sin(theta / 2.0)
        aa, bb, cc, dd = a * a, b * b, c * c, d * d
        bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
        return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                         [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                         [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

    # 3D GEOMETRY METHODS ##############################################################################################
    @staticmethod
    def rotation_z(theta):
        return np.array([[np.cos(theta), -np.sin(theta), 0],
                         [np.sin(theta), np.cos(theta), 0],
                         [0, 0, 1]])

    @staticmethod
    def unit_vector(vector):
        """ Returns the unit vector of the vector.  """
        return vector / np.linalg.norm(vector)

    def angle_between(self, v1, v2):
        """ Returns the angle in radians between vectors 'v1' and 'v2' """
        v1_u = self.unit_vector(v1)
        v2_u = self.unit_vector(v2)
        return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

    def angle_between_cc(self, v1, v2):
        """ Returns the angle clockwise in z from v1 to v2."""
        a_i = self.angle_between(v1, v2)
        # rotation check
        v2_n = np.matmul(self.rotation_z(0.01), deepcopy(v2))
        a_new = self.angle_between(v1, v2_n)
        if a_new > a_i:  # Clockwise measurement
            angle = a_i
        else:
            angle = 2 * np.pi - a_i  # Adjust for clockwise measurement
        return angle

    # 3D ORBIT METHODS #################################################################################################
    @staticmethod
    def coplanar_2_j2000(aop, lan, inc):
        taop = np.array([[  np.cos(aop), np.sin(aop), 0],
                        [- np.sin(aop), np.cos(aop), 0],
                        [            0,           0, 1]])

        tinc = np.array([[1,             0,            0],
                        [0,   np.cos(inc), np.sin(inc)],
                        [0, - np.sin(inc), np.cos(inc)]])

        tlan = np.array([[  np.cos(lan), np.sin(lan), 0],
                        [- np.sin(lan), np.cos(lan), 0],
                        [            0,           0, 1]])
        return np.matmul(taop, np.matmul(tinc, tlan))

    @staticmethod
    def polytime_2_datetime(_time):
        temp = copy(_time)
        temp.format = 'datetime'
        return temp.value

    # LAMBERT SOLUTION METHOD ##########################################################################################
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
        self.N = 500

    # EQUATIONS DEFINED BY ELIZABETTA IORFIDA  TODO: FINISH REFERENCE [21]
    @staticmethod
    def _alpha(v_inf_i, v_inf_f, r_p, mu_body):
        return np.arcsin(1/(1+(r_p*np.square(np.linalg.norm(v_inf_i.to(u.km/u.s))))*(u.km/u.s)**2/mu_body)) + \
               np.arcsin(1/(1+(r_p*np.square(np.linalg.norm(v_inf_f.to(u.km/u.s))))*(u.km/u.s)**2/mu_body))

    # EQUATIONS DEFINED BY ELIZABETTA IORFIDA  TODO: FINISH REFERENCE [21]
    @staticmethod
    def _f_e_i(e_i_n, a_i, a_f, alpha_req):
        return np.arcsin(1/e_i_n) + np.arcsin(1/(1-(a_i/a_f)*(1-e_i_n))) - alpha_req

    @staticmethod
    def _d_f_e_i(e_i_n, a_i, a_f):
        def const_form(a, b, c):
            return - a / (np.sqrt(1 - (b / np.square(c))) * np.square(c))
        t1 = const_form(1, 1, e_i_n)
        t2 = const_form(a_i * a_f, np.square(a_f), a_f + a_i * (e_i_n - 1))
        return t1 + t2

    @staticmethod
    def hyp_a(v_inf, mu):
        return (- mu / (np.square(np.linalg.norm(v_inf.to(u.km/u.s)))*(u.km/u.s)**2)).to(u.km)

    def newton_rhapson_pga(self, v_inf_i, v_inf_f, mu_body, alpha_required):
        e_i_0 = 1.1
        a_i = self.hyp_a(v_inf_i, mu_body)
        a_f = self.hyp_a(v_inf_f, mu_body)

        alpha_required = alpha_required

        def func(_e_i):
            return self._f_e_i(_e_i, a_i=a_i.to(u.km).value, a_f=a_f.to(u.km).value,
                               alpha_req=alpha_required)

        def _fprime(_e_i):
            return self._d_f_e_i(_e_i, a_i=a_i.to(u.km).value, a_f=a_f.to(u.km).value)

        e_i = optimize.newton(func, e_i_0, _fprime)
        r_p = a_i * (1 - e_i)
        e_f = - (r_p / a_f - 1)
        v_p_i = np.sqrt(np.square(np.linalg.norm(v_inf_i)) * (e_i + 1) / (e_i - 1))
        v_p_f = np.sqrt(np.square(np.linalg.norm(v_inf_i)) * (e_f + 1) / (e_f - 1))
        dv = np.linalg.norm(v_p_f - v_p_i)
        return e_i, e_f, a_i, a_f, v_p_i, v_p_f, dv, r_p

    def powered_gravity_assist(self, _itinerary_data_indexed, mode='scalar'):
        rsoi    = body_d_domain[self._body_string_lower(_itinerary_data_indexed['b'])]['upper']
        mu_body = _itinerary_data_indexed['b'].k
        v_inf_i = _itinerary_data_indexed['v']['a']
        v_inf_f = _itinerary_data_indexed['v']['d']
        alpha_required = self.angle_between(v_inf_i, v_inf_f)
        e_i, e_f, a_i, a_f, v_p_i, v_p_f, rp_dv, r_p = self.newton_rhapson_pga(v_inf_i, v_inf_f, mu_body, alpha_required)

        # _itinerary_data_indexed['dv'] = rp_dv  TODO: Implement properly

        if mode is 'scalar_evaluation':
            _gravass_params = gravass_parameters(type='scalar',
                                                 a_i_mag=a_i.to(u.km),
                                                 a_f_mag=a_f.to(u.km),
                                                 e_i_mag=e_i,
                                                 e_i_vec=None,
                                                 e_f_vec=None,
                                                 e_f_mag=e_f,
                                                 v_inf_i_vec=v_inf_i.to(u.km / u.s),
                                                 v_inf_f_vec=v_inf_f.to(u.km / u.s),
                                                 v_planet_i_vec=_itinerary_data_indexed['v']['p'].to(u.km / u.s),
                                                 v_planet_f_vec=_itinerary_data_indexed['v']['p'].to(u.km / u.s),
                                                 r_p_dv_mag=rp_dv * (u.km / u.s),
                                                 v_p_i_vec=v_p_i * (u.km / u.s),
                                                 v_p_f_vec=v_p_f * (u.km / u.s),
                                                 t_p_i=0,
                                                 t_p_f=0,
                                                 aop=None,
                                                 lan=None,
                                                 inc=None)

        elif mode is 'vector_evaluation':

            # Orbital plane.
            n_vec_orbital = self.unit_vector(np.cross(v_inf_i, v_inf_f))

            # Rotation about orbital plane normal vector for v_p_unit_vec with angle of d_i.
            d_i = 2 * np.arcsin(1 / e_i)
            v_p_unit_vec = self.unit_vector(np.dot(self.rotation_matrix(axis=n_vec_orbital,
                                                                        theta=d_i), v_inf_i))

            # v_p_i_vec and v_p_f_vec
            v_p_i_vec = v_p_i * v_p_unit_vec
            v_p_f_vec = v_p_f * v_p_unit_vec

            # r_p_unit_vec and r_p_vec
            r_p_unit_vec = self.unit_vector(np.dot(self.rotation_matrix(axis=n_vec_orbital,
                                                                        theta=-np.pi / 2), v_p_i_vec).value)
            r_p_vec = r_p * r_p_unit_vec

            # eccentricity vectors
            e_i_vec = np.cross(v_p_i_vec, np.cross(r_p_vec, v_p_i_vec)) / \
                      _itinerary_data_indexed['b'].k.to(u.km ** 3 / u.s ** 2).value - r_p_unit_vec
            e_f_vec = np.cross(v_p_f_vec, np.cross(r_p_vec, v_p_f_vec)) / \
                      _itinerary_data_indexed['b'].k.to(u.km ** 3 / u.s ** 2).value - r_p_unit_vec

            # Classical orbit parameters
            inclination = np.arccos(
                np.dot(np.array([0, 0, 1]), n_vec_orbital) / (np.linalg.norm(n_vec_orbital)))

            n_vec = np.cross(np.array([0, 0, 1]), np.cross(r_p_vec, v_p_i_vec))
            lan = np.arccos(np.dot(np.array([1, 0, 0]), n_vec) / (np.linalg.norm(n_vec) * 1))

            if n_vec[1] < 0:
                lan = 2 * np.pi - lan

            aop = np.arccos(np.dot(n_vec, e_i_vec) / np.linalg.norm(n_vec) / np.linalg.norm(e_i_vec))

            if e_i_vec[-1] < 0:
                aop = 2 * np.pi - aop

            theta_inf_i = np.arccos(-1 / e_i)
            theta_inf_f = np.arccos(-1 / e_f).value

            if n_vec_orbital[-1]<0:
                theta_inf_i = -theta_inf_i
                theta_inf_f = -theta_inf_f

            H_rsoi_i = np.arcsinh(rsoi
                                  * np.sin(theta_inf_i)/(a_i.value * np.sqrt(e_i**2 - 1)))
            H_rsoi_f = np.arcsinh(rsoi
                                  * np.sin(theta_inf_f) / (a_f.value * np.sqrt(e_f.value ** 2 - 1))) * -1

            mu = _itinerary_data_indexed['b'].k.to(u.km ** 3 / u.s ** 2).value

            t_rsoi_i = np.sqrt((-a_i)**3/mu).value * (e_i * np.sinh(H_rsoi_i)-H_rsoi_i)
            t_rsoi_f = np.sqrt((-a_f)**3/mu).value * (e_f * np.sinh(H_rsoi_f)-H_rsoi_f)

            epoch_rp = _itinerary_data_indexed['d']
            epoch_entry = _itinerary_data_indexed['d'] + time.TimeDelta(t_rsoi_i * u.s)
            epoch_exit = _itinerary_data_indexed['d'] + time.TimeDelta(t_rsoi_f * u.s)

            _itinerary_data_indexed['d'] = {'rp': epoch_rp}
            _itinerary_data_indexed['d']['i'] = epoch_entry
            _itinerary_data_indexed['d']['f'] = epoch_exit

            op = OrbitPlotter3D()
            op.set_attractor(Jupiter)
            ss_i = Orbit.from_classical(attractor=Jupiter, a=a_i, ecc=e_i * u.one, inc=inclination * u.rad,
                                          raan=lan * u.rad, argp=aop * u.rad, nu=0*u.rad, epoch=epoch_rp)

            ss_f = Orbit.from_classical(attractor=Jupiter, a=a_f, ecc=e_f*u.one, inc=inclination*u.rad,
                                          raan=lan*u.rad, argp=aop*u.rad, nu=0*u.rad, epoch=epoch_rp)

            epoch_entry_dt = self.polytime_2_datetime(_itinerary_data_indexed['d']['i'])
            epoch_exit_dt = self.polytime_2_datetime(_itinerary_data_indexed['d']['f'])
            epoch_rp_dt = self.polytime_2_datetime(_itinerary_data_indexed['d']['rp'])

            tv_ent = time_range(epoch_entry_dt, periods=100, spacing=None, end=epoch_rp_dt)
            tv_ext = time_range(epoch_rp_dt, periods=100, spacing=None, end=epoch_exit_dt)

            op.plot_trajectory(ss_i.sample(tv_ent)[-1], label='Entry')
            op.plot_trajectory(ss_f.sample(tv_ext)[-1], label='Exit')

            op._data.append(create_soi(rsoi))
            op.plot(ss_i)
            op.plot(ss_f)

            plt.xlim(-0.5*rsoi, 0.5*rsoi)
            plt.ylim(-0.5*rsoi, 0.5*rsoi)

            op.set_view(30 * u.deg, 260 * u.deg, distance=3 * u.km)
            op.savefig("EJPExample.png", title="EJP Optimal trajectory sequence")

            _gravass_params = gravass_parameters(type='vector',
                                                 a_i_mag=a_i.to(u.km),
                                                 a_f_mag=a_f.to(u.km),
                                                 e_i_mag=e_i,
                                                 e_i_vec=e_i_vec,
                                                 e_f_mag=e_f,
                                                 e_f_vec=e_f_vec,
                                                 v_inf_i_vec=v_inf_i.to(u.km / u.s),
                                                 v_inf_f_vec=v_inf_f.to(u.km / u.s),
                                                 v_planet_i_vec=Orbit.from_body_ephem(
                                                     _itinerary_data_indexed['b'].state.v, epoch_entry),
                                                 v_planet_f_vec=Orbit.from_body_ephem(
                                                     _itinerary_data_indexed['b'].state.v, epoch_entry),
                                                 r_p_dv_mag=rp_dv * (u.km / u.s),
                                                 v_p_i_vec=v_p_i * (u.km / u.s),
                                                 v_p_f_vec=v_p_f * (u.km / u.s),
                                                 t_p_i=t_rsoi_i,
                                                 t_p_f=t_rsoi_f,
                                                 aop=aop,
                                                 lan=lan,
                                                 inc=inclination)
        else:
            raise AttributeError("Mode is not recognised.")
        return _gravass_params

    # PRELIMINARY STATIONARY ANALYSIS
    def stationary_process_itinerary(self, _raw_itinerary, _body_list, mode='fast', verbose=False):
        """

        :param _raw_itinerary:
        :param _body_list:
        :param verbose
        :return:
        """

        # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        if verbose:
            print('\n')
            print('-' * 40 + '-' * len(' ID: {}'.format(_raw_itinerary['id'])))
            print('Initializing...'.ljust(40) + ' ID: {}\n'.format(_raw_itinerary['id']))
        # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

        # GENERATE PROCESSED ITINERARY STRUCTURE -----------------------------------------------------------------------
        _itinerary_data = [None] * (len(_raw_itinerary['durations']) + 1)

        _itinerary_data[0] = {'b': body_list[_body_list[0]],
                              'd': time.Time(_raw_itinerary['launch_date'], scale='tdb'),
                              's': _body_list[0],
                              'v': {}}

        # SET EPOCH FOR EACH BODY ENCOUNTER ACCORDING TO LEG DISTRIBUTION ----------------------------------------------
        for i in range(len(_raw_itinerary['durations'])):
            _itinerary_data[i + 1] = {'b': body_list[_body_list[i + 1]],
                                      'd': time.Time(_raw_itinerary['launch_date'] +
                                                     datetime.timedelta(
                                                         days=365 * sum(_raw_itinerary['durations'][:i + 1])),
                                                     scale='tdb'),
                                      's': _body_list[i + 1],
                                      'v': {}}

        # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        if verbose:                                                                                                  # $
            print('Solving Lambert multi-leg problem...'.ljust(40) + ' ID: {}\n'.format(_raw_itinerary['id']))       # $
        # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

        ################################################################################################################
        # MULTIPLE-LEG LAMBERT SOLUTION -------------------------------------------------------------------------------#
        ################################################################################################################
        for i in range(len(_raw_itinerary['durations'])):
            _itinerary_data[i + 1]['l'] = self.lambert_solve_from_bodies(_itinerary_data[i]['b'],
                                                                         _itinerary_data[i + 1]['b'],
                                                                         _itinerary_data[i]['d'],
                                                                         _itinerary_data[i + 1]['d'])

        ################################################################################################################
        # GRAVITY ASSIST FEASIBILITY ----------------------------------------------------------------------------------#
        ################################################################################################################
        for i in range(len(_raw_itinerary['durations'])-1):
                _itinerary_data[i + 1]['v']['a'] = _itinerary_data[i + 1]['l'].v1.to(u.km / u.s)
                _itinerary_data[i + 1]['v']['p'] = _itinerary_data[i + 1]['l'].ss1.state.v.to(u.km / u.s)
                _itinerary_data[i + 1]['v']['d'] = _itinerary_data[i + 2]['l'].v0.to(u.km / u.s)

                v_inf_i = (_itinerary_data[i + 1]['v']['a'] - _itinerary_data[i + 1]['v']['p']).to(u.km / u.s)
                v_inf_f = (_itinerary_data[i + 1]['v']['d'] - _itinerary_data[i + 1]['v']['p']).to(u.km / u.s)

                alpha_required = self.angle_between(v_inf_i, v_inf_f)
                alpha_max = self._alpha(v_inf_i=v_inf_i,
                                        v_inf_f=v_inf_f,
                                        r_p=body_d_domain[
                                                self._body_string_lower(_itinerary_data[i+1]['b'])]['lower']*(u.km),
                                        mu_body=_itinerary_data[i+1]['b'].k)

                # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if alpha_required * u.rad <= alpha_max:
                    pass
                else:
                    raise ValueError("Gravity assist required bending angle is not possible to meet given the \n"
                                     "requirement of alpha: {} <= alpha_max: {}".format(alpha_required, alpha_max))
                # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                if (mode is 'scalar_evaluation') or (mode is 'vector_evaluation'):
                    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                    if verbose:                                                                                      # $
                        print('Calculating gravity assist parameters...'.ljust(40) + ' ID: {}\n'.format(             # $
                            _raw_itinerary['id']))                                                                   # $
                    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



                    e_i, e_f, a_i, a_f, v_p_i, v_p_f, rp_dv, r_p = self.newton_rhapson_pga(v_inf_i, v_inf_f, mu_body)

                    _itinerary_data[i + 1]['dv'] = rp_dv

                    if mode is 'scalar_evaluation':
                        _gravass_params = gravass_parameters(type='scalar',
                                                             a_i_mag=a_i.to(u.km),
                                                             a_f_mag=a_f.to(u.km),
                                                             e_i_mag=e_i,
                                                             e_i_vec=None,
                                                             e_f_vec=None,
                                                             e_f_mag=e_f,
                                                             v_inf_i_vec=v_inf_i.to(u.km/u.s),
                                                             v_inf_f_vec=v_inf_f.to(u.km/u.s),
                                                             v_planet_i_vec=_itinerary_data[i + 1]['v']['p'].to(u.km / u.s),
                                                             v_planet_f_vec=_itinerary_data[i + 1]['v']['p'].to(u.km / u.s),
                                                             r_p_dv_mag=rp_dv*(u.km / u.s),
                                                             v_p_i_vec=v_p_i*(u.km / u.s),
                                                             v_p_f_vec=v_p_f*(u.km / u.s),
                                                             t_p_i=0,
                                                             t_p_f=0,
                                                             aop=None,
                                                             lan=None,
                                                             inc=None)

                        #
                        _itinerary_data[i + 1]['ga'] = _gravass_params

                    if mode is 'vector_evaluation':
                        I = np.array([1, 0, 0])
                        J = np.array([0, 1, 0])
                        K = np.array([1, 0, 1])

                        # Orbital plane.
                        n_vec_orbital = self.unit_vector(np.cross(v_inf_i, v_inf_f))

                        # Rotation about orbital plane normal vector for v_p_unit_vec with angle of d_i.
                        d_i = 2 * np.arcsin(1 / e_i)
                        v_p_unit_vec = self.unit_vector(np.dot(self.rotation_matrix(axis=n_vec_orbital,
                                                                                    theta=d_i), v_inf_i))

                        # v_p_i_vec and v_p_f_vec
                        v_p_i_vec = v_p_i * v_p_unit_vec
                        v_p_f_vec = v_p_f * v_p_unit_vec

                        # r_p_unit_vec and r_p_vec
                        r_p_unit_vec = self.unit_vector(np.dot(self.rotation_matrix(axis=n_vec_orbital,
                                                                                    theta=-np.pi / 2), v_p_i_vec).value)
                        r_p_vec = r_p * r_p_unit_vec

                        # eccentricity vectors
                        e_i_vec = np.cross(v_p_i_vec, np.cross(r_p_vec, v_p_i_vec)) / \
                                  _itinerary_data[i + 1]['b'].k.to(u.km ** 3 / u.s ** 2).value - r_p_unit_vec
                        e_f_vec = np.cross(v_p_i_vec, np.cross(r_p_vec, v_p_f_vec)) / \
                                  _itinerary_data[i + 1]['b'].k.to(u.km ** 3 / u.s ** 2).value - r_p_unit_vec

                        # Classical orbit parameters
                        inclination = np.arccos(
                            np.dot(np.array([0, 0, 1]), n_vec_orbital) / (np.linalg.norm(n_vec_orbital)))

                        n_vec = np.cross(np.array([0, 0, 1]), np.cross(r_p_vec, v_p_i_vec))
                        lan = np.arccos(np.dot(np.array([1, 0, 0]), n_vec) / (np.linalg.norm(n_vec) * 1))

                        if n_vec[1] < 0:
                            lan = 2 * np.pi - lan

                        aop = np.arccos(np.dot(n_vec, e_i_vec) / np.linalg.norm(n_vec) / np.linalg.norm(e_i_vec))

                        if e_i_vec[-1] < 0:
                            aop = 2 * np.pi - aop

                        _gravass_params = gravass_parameters(type='vector',
                                                             a_i_mag=a_i.to(u.km),
                                                             a_f_mag=a_f.to(u.km),
                                                             e_i_mag=e_i,
                                                             e_i_vec=e_i_vec,
                                                             e_f_mag=e_f,
                                                             e_f_vec=e_f_vec,
                                                             v_inf_i_vec=v_inf_i.to(u.km / u.s),
                                                             v_inf_f_vec=v_inf_f.to(u.km / u.s),
                                                             v_planet_i_vec=_itinerary_data[i + 1]['v']['p'].to(
                                                                 u.km / u.s),
                                                             v_planet_f_vec=_itinerary_data[i + 1]['v']['p'].to(
                                                                 u.km / u.s),
                                                             r_p_dv_mag=rp_dv * (u.km / u.s),
                                                             v_p_i_vec=v_p_i * (u.km / u.s),
                                                             v_p_f_vec=v_p_f * (u.km / u.s),
                                                             t_p_i=0,
                                                             t_p_f=0,
                                                             aop=aop,
                                                             lan=lan,
                                                             inc=inclination)
                        _itinerary_data[i + 1]['ga'] = _gravass_params

                    # op = OrbitPlotter()
                    #
                    # # @u.quantity_input(a=u.m, ecc=u.one, inc=u.rad, raan=u.rad, argp=u.rad, nu=u.rad)
                    # # def from_classical(cls, attractor, a, ecc, inc, raan, argp, nu, epoch=J2000):
                    # op.plot(Orbit.from_classical(Jupiter, a_f, e_f, inc=inclination*u.rad, raan=lan*u.rad, argp=aop*u.rad, nu=0*u.rad))
                    # plt.show()

        if (mode is 'scalar_evaluation') or (mode is 'vector_evaluation') or (mode is 'plot2D') or (mode is 'plot3D'):
            # DEPARTURE BODY DATA --------------------------------------------------------------------------------------
            _itinerary_data[0]['v']['p'] = _itinerary_data[1]['l'].ss0.state.v.to(u.km / u.s)
            _itinerary_data[0]['v']['d'] = _itinerary_data[1]['l'].v0.to(u.km / u.s)

            v_inf = np.linalg.norm((_itinerary_data[0]['v']['d'] - _itinerary_data[0]['v']['p']).to(
                    u.km / u.s).value)

            body = _itinerary_data[0]['b']

            r_0 = epo['alt'] + body.R.to(u.km).value
            v_0 = np.sqrt(np.square(v_inf*(u.km/u.s))+2*body.k.to(u.km**3/u.s**2)/(r_0*(u.km)))
            _itinerary_data[0]['dv'] = (v_0 - np.sqrt(body.k.to(u.km**3/u.s**2)/(r_0*(u.km)))).value

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

        ################################################################################################################
        # PLOTTING                                                                                                     #
        ################################################################################################################
        if (mode is 'plot2D') or (mode is 'plot3D'):

            # TRAJECTORIES OF LEGS -------------------------------------------------------------------------------------
            for i in range(len(_raw_itinerary['durations'])):
                _itinerary_data[i + 1]['t'] = Orbit.from_vectors(_itinerary_data[i + 1]['l'].attractor,
                                                                 _itinerary_data[i + 1]['l'].r0,
                                                                 _itinerary_data[i + 1]['l'].v0)

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

            if mode is 'plot3D':
                frame = OrbitPlotter3D()
                frame.set_attractor(Sun)
            elif mode is 'plot2D':
                frame = OrbitPlotter()
                frame.set_attractor(Sun)
            else:
                raise ValueError

            for i in range(len(_raw_itinerary['durations'])):
                frame.plot(_itinerary_data[i + 1]['l'].ss0, color='0.8')

            frame.plot(_itinerary_data[len(_raw_itinerary['durations'])]['l'].ss1, color='0.8')

            for i in range(len(_raw_itinerary['durations']) + 1):
                frame.plot_trajectory(_itinerary_plot_data[i]['rr'], label=_itinerary_data[i]['b'],
                                      color=color_legs2[i])

            for i in range(len(_raw_itinerary['durations'])):
                frame.plot_trajectory(
                    _itinerary_plot_data[i + 1]['tp'].sample(_itinerary_data[i + 1]['tv'])[-1],
                    label="Leg {}".format(i + 1),
                    color=color_legs[i+1],
                )

            frame._redraw_attractor(0.25 * 10 ** (8) * u.km)
            print('Displaying plot!'.ljust(40) + ' ID: {}\n'.format(_raw_itinerary['id']))

            if mode is 'plot2D':
                plt.legend()
                plt.show()

            else:
                frame.set_view(30 * u.deg, 260 * u.deg, distance=3 * u.km)
                # frame.show(title="EJP Example")
                frame.savefig("EJPExample.png", title="EJP Optimal trajectory sequence")

        # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        if verbose:
            print('Complete!'.ljust(40) + ' ID: {}'.format(_raw_itinerary['id']))
            print('-' * 40 + '-' * len(' ID: {}\n'.format(_raw_itinerary['id'])))
        # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

        return _itinerary_data

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
                                      color='green')

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
        processed = _test.stationary_process_itinerary(__raw_itinerary2, __raw_itinerary1, mode='full')

        # print(processed)

        _test.powered_gravity_assist(processed[1], mode='vector_evaluation')

        # total = [24-float(i) for i in np.linspace(0, 9, 100)]

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
