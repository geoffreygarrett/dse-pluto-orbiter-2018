import numpy as np
import astropy.units as u
from scipy import optimize
from .geometry import *
from trajectory_tool.mnag_mission_analysis.data_structures import *
import pandas as pd
from datetime import datetime, timedelta
from sympy import Eq, var, solve
import sympy

def alpha(v_inf_i, v_inf_f, r_p, mu_body):
    return np.arcsin(
        1 / (1 + (r_p * np.square(np.linalg.norm(v_inf_i.to(u.km / u.s)))) * (u.km / u.s) ** 2 / mu_body)) + \
           np.arcsin(1 / (1 + (r_p * np.square(np.linalg.norm(v_inf_f.to(u.km / u.s)))) * (u.km / u.s) ** 2 / mu_body))


# EQUATIONS DEFINED BY ELIZABETTA IORFIDA  TODO: FINISH REFERENCE [21]
def f_e_i(e_i_n, a_i, a_f, alpha_req):
    return np.arcsin(1 / e_i_n) + np.arcsin(1 / (1 - (a_i / a_f) * (1 - e_i_n))) - alpha_req


def d_f_e_i(e_i_n, a_i, a_f):
    def const_form(a, b, c):
        return - a / (np.sqrt(1 - (b / np.square(c))) * np.square(c))

    t1 = const_form(1, 1, e_i_n)
    t2 = const_form(a_i * a_f, np.square(a_f), a_f + a_i * (e_i_n - 1))
    return t1 + t2


def _hyp_sma(v_inf, mu):
    return (- mu / (np.square(np.linalg.norm(v_inf.to(u.km / u.s))) * (u.km / u.s) ** 2)).to(u.km)


# def _theta_rsoi(r_p, r_soi, e):
#     return np.arccos((1/e)* ((r_p*(1+e))/(r_soi))-1)

def _theta_rsoi(r_p, r_soi, e):
    return np.arccos((2*r_p)/(e*r_soi) - 1/e)


# VERIFIED
def F(e, theta):
    num = (np.sqrt(e+1) + np.sqrt(e-1) * np.tan(theta/2))
    den = (np.sqrt(e+1) - np.sqrt(e-1) * np.tan(theta/2))
    return np.log(num/den)


# VERIFIED
def Mn(e, F):
    return e*np.sinh(F) - F


# VERIFIED
def tp(Mn, h, e, mu):
    return Mn * ((h**3)/(mu**2 * (e**2 - 1) ** (3/2)))


def _newton_rhapson(v_inf_i, v_inf_f, mu_body, alpha_required):
    e_i_0 = 1.1
    sma_i = _hyp_sma(v_inf_i, mu_body)
    sma_f = _hyp_sma(v_inf_f, mu_body)
    def func(_e_i):
        return f_e_i(_e_i, a_i=sma_i.to(u.km).value, a_f=sma_f.to(u.km).value,
                     alpha_req=alpha_required)
    def _fprime(_e_i):
        return d_f_e_i(_e_i, a_i=sma_i.to(u.km).value, a_f=sma_f.to(u.km).value)
    ecc_i = optimize.newton(func, e_i_0, _fprime)
    r_p = sma_i * (1 - ecc_i)
    ecc_f = - (r_p / sma_f - 1)
    v_p_i = np.sqrt(np.square(np.linalg.norm(v_inf_i)) * (ecc_i + 1) / (ecc_i - 1))
    v_p_f = np.sqrt(np.square(np.linalg.norm(v_inf_f)) * (ecc_f + 1) / (ecc_f - 1))
    r_p_dv = np.linalg.norm(v_p_f - v_p_i)
    return ecc_i, ecc_f, sma_i, sma_f, v_p_i, v_p_f, r_p_dv, r_p


def basic_flyby(_base_gravity_assist, planetary_node):
    basic_attributes = flyby_basic
    basic_attributes.v_i, basic_attributes.v_f, basic_attributes.v_planet_i, basic_attributes.v_planet_f, \
        basic_attributes.v_inf_i, basic_attributes.v_inf_f, basic_attributes.a_req = _base_gravity_assist
    args = _newton_rhapson(basic_attributes.v_inf_i, basic_attributes.v_inf_f,
                           planetary_node.body.k, basic_attributes.a_req)
    basic_attributes.ecc_i = args[0]
    basic_attributes.ecc_f = args[1]
    basic_attributes.sma_i = args[2]
    basic_attributes.sma_f = args[3]
    basic_attributes.v_p_i = args[4]
    basic_attributes.v_p_f = args[5]
    basic_attributes.r_p_dv = args[6]
    basic_attributes.r_p = args[7]
    return basic_attributes


def basic_flyby_df(_flyby_basic: flyby_basic, unit=False):
    columns = ['IN_OUT', 'v', 'v_inf', 'v_planet', 'ecc', 'sma', 'v_p', 'r_p', 'dv']
    df_in = pd.DataFrame(columns=columns)
    df_out = pd.DataFrame(columns=columns)

    if unit:
        df_in.IN_OUT = ['IN']
        df_in.v = [np.linalg.norm(_flyby_basic.v_i.to(u.km/u.s)) * (u.km/u.s)]
        df_in.v_inf = [np.linalg.norm(_flyby_basic.v_inf_i.to(u.km/u.s)) * (u.km/u.s)]
        df_in.v_planet = [np.linalg.norm(_flyby_basic.v_planet_i.to(u.km/u.s)) * (u.km/u.s)]
        df_in.ecc = [_flyby_basic.ecc_i]
        df_in.sma = [_flyby_basic.sma_i.to(u.km)]
        df_in.v_p = [_flyby_basic.v_p_i * u.km/u.s]
        df_in.r_p = [_flyby_basic.r_p.to(u.km)]
        df_in.dv = [_flyby_basic.r_p_dv * u.km/u.s]

        df_out.IN_OUT = ['OUT']
        df_out.v = [np.linalg.norm(_flyby_basic.v_f.to(u.km/u.s)) * (u.km/u.s)]
        df_out.v_inf = [np.linalg.norm(_flyby_basic.v_inf_f.to(u.km/u.s)) * (u.km/u.s)]
        df_out.v_planet = [np.linalg.norm(_flyby_basic.v_planet_f.to(u.km/u.s)) * (u.km/u.s)]
        df_out.ecc = [_flyby_basic.ecc_f]
        df_out.sma = [_flyby_basic.sma_f.to(u.km)]
        df_out.v_p = [_flyby_basic.v_p_f * u.km/u.s]
        df_out.r_p = [_flyby_basic.r_p.to(u.km)]
        df_out.dv = [_flyby_basic.r_p_dv * u.km/u.s]

    else:
        df_in.IN_OUT = ['IN']
        df_in.v = [np.linalg.norm(_flyby_basic.v_i.to(u.km / u.s))]
        df_in.v_inf = [np.linalg.norm(_flyby_basic.v_inf_i.to(u.km / u.s))]
        df_in.v_planet = [np.linalg.norm(_flyby_basic.v_planet_i.to(u.km / u.s))]
        df_in.ecc = [_flyby_basic.ecc_i]
        df_in.sma = [_flyby_basic.sma_i.to(u.km).value]
        df_in.v_p = [_flyby_basic.v_p_i]
        df_in.r_p = [_flyby_basic.r_p.to(u.km).value]
        df_in.dv = [_flyby_basic.r_p_dv]

        df_out.IN_OUT = ['OUT']
        df_out.v = [np.linalg.norm(_flyby_basic.v_f.to(u.km / u.s))]
        df_out.v_inf = [np.linalg.norm(_flyby_basic.v_inf_f.to(u.km / u.s))]
        df_out.v_planet = [np.linalg.norm(_flyby_basic.v_planet_f.to(u.km / u.s))]
        df_out.ecc = [_flyby_basic.ecc_f]
        df_out.sma = [_flyby_basic.sma_f.to(u.km).value]
        df_out.v_p = [_flyby_basic.v_p_f]
        df_out.r_p = [_flyby_basic.r_p.to(u.km).value]
        df_out.dv = [_flyby_basic.r_p_dv]
    return df_in.append(df_out)


def guess_flyby(basic_attributes, planetary_node):
    guess_attributes = flyby_guess
    args = pga_scalar_2_vector(basic_attributes.v_inf_i, basic_attributes.v_inf_f, basic_attributes.v_p_i,
                               basic_attributes.v_p_f, basic_attributes.sma_i, basic_attributes.sma_f,
                               basic_attributes.ecc_i, basic_attributes.ecc_f, basic_attributes.r_p,
                               planetary_node.body, planetary_node.soi_periapsis_magnitude, planetary_node)
    guess_attributes.v_i = basic_attributes.v_i
    guess_attributes.v_f = basic_attributes.v_f
    guess_attributes.v_inf_i = basic_attributes.v_inf_i
    guess_attributes.v_inf_f = basic_attributes.v_inf_f
    guess_attributes.v_planet_i = planetary_node.v_planet_i
    guess_attributes.v_planet_f = planetary_node.v_planet_f
    guess_attributes.sma_i = basic_attributes.sma_i
    guess_attributes.sma_f = basic_attributes.sma_f
    guess_attributes.r_p = basic_attributes.r_p
    guess_attributes.r_p_dv = basic_attributes.r_p_dv
    guess_attributes.v_p_i = args[0]
    guess_attributes.v_p_f = args[1]
    guess_attributes.raan = args[2]
    guess_attributes.argp = args[3]
    guess_attributes.inc = args[4]
    guess_attributes.t_i = args[5]
    guess_attributes.t_f = args[6]
    guess_attributes.ecc_i = args[7]
    guess_attributes.ecc_f = args[8]
    guess_attributes.r_entry = args[9]
    guess_attributes.r_exit = args[10]
    return guess_attributes


def refine_flyby(_base_gravity_assist, planetary_node):

    basic_attributes = flyby_basic
    basic_attributes.v_i, basic_attributes.v_f, basic_attributes.v_planet_i, basic_attributes.v_planet_f, \
        basic_attributes.v_inf_i, basic_attributes.v_inf_f, basic_attributes.a_req = _base_gravity_assist
    args = _newton_rhapson(basic_attributes.v_inf_i, basic_attributes.v_inf_f,
                           planetary_node.body.k, basic_attributes.a_req)
    basic_attributes.ecc_i = args[0]
    basic_attributes.ecc_f = args[1]
    basic_attributes.sma_i = args[2]
    basic_attributes.sma_f = args[3]
    basic_attributes.v_p_i = args[4]
    basic_attributes.v_p_f = args[5]
    basic_attributes.r_p_dv = args[6]
    basic_attributes.r_p = args[7]

    _refined_attributes = flyby_refined
    args = pga_scalar_2_vector(basic_attributes.v_inf_i, basic_attributes.v_inf_f, basic_attributes.v_p_i,
                               basic_attributes.v_p_f, basic_attributes.sma_i, basic_attributes.sma_f,
                               basic_attributes.ecc_i, basic_attributes.ecc_f, basic_attributes.r_p,
                               planetary_node.body, planetary_node.soi_periapsis_magnitude, planetary_node)

    if _refined_attributes.error_v == 0:
        _refined_attributes.error_v = 1
    if _refined_attributes.error_p == 0:
        _refined_attributes.error_p = 1

    _refined_attributes.v_i = basic_attributes.v_i
    _refined_attributes.v_f = basic_attributes.v_f
    _refined_attributes.v_inf_i = basic_attributes.v_inf_i
    _refined_attributes.v_inf_f = basic_attributes.v_inf_f
    _refined_attributes.v_planet_i = planetary_node.v_planet_i
    _refined_attributes.v_planet_f = planetary_node.v_planet_f
    _refined_attributes.sma_i = basic_attributes.sma_i
    _refined_attributes.sma_f = basic_attributes.sma_f
    _refined_attributes.r_p = basic_attributes.r_p
    _refined_attributes.r_p_dv = basic_attributes.r_p_dv
    _refined_attributes.v_p_i = args[0]
    _refined_attributes.v_p_f = args[1]
    _refined_attributes.raan = args[2]
    _refined_attributes.argp = args[3]
    _refined_attributes.inc = args[4]
    _refined_attributes.t_i = args[5]
    _refined_attributes.t_f = args[6]
    _refined_attributes.ecc_i = args[7]
    _refined_attributes.ecc_f = args[8]
    _refined_attributes.r_entry = args[9]
    _refined_attributes.r_exit = args[10]
    return _refined_attributes


def guess_flyby_df(_flyby_guess: flyby_guess, unit=False):
    columns = ['IN_OUT', 'v_x', 'v_y', 'v_z', 'v_inf_x', 'v_inf_y', 'v_inf_z', 'v_planet_x', 'v_planet_y', 'v_planet_z',
               'ecc', 'ecc_x', 'ecc_y', 'ecc_z', 'sma', 'argp', 'inc', 'raan', 't_days', 'r_soi_x', 'r_soi_y', 'r_soi_z',
               'v_p_x', 'v_p_y', 'v_p_z', 'r_p', 'dv']
    df_in = pd.DataFrame(columns=columns)
    df_out = pd.DataFrame(columns=columns)

    if unit:
        df_in.IN_OUT = ['IN']
        df_in.v_x = [_flyby_guess.v_i[0]]
        df_in.v_y = [_flyby_guess.v_i[1]]
        df_in.v_z = [_flyby_guess.v_i[2]]
        df_in.v_inf_x = [_flyby_guess.v_inf_i[0]]
        df_in.v_inf_y = [_flyby_guess.v_inf_i[1]]
        df_in.v_inf_z = [_flyby_guess.v_inf_i[2]]
        df_in.v_planet_x = [_flyby_guess.v_planet_i[0]]
        df_in.v_planet_y = [_flyby_guess.v_planet_i[1]]
        df_in.v_planet_z = [_flyby_guess.v_planet_i[2]]
        df_in.ecc = [np.linalg.norm(_flyby_guess.ecc_i)]
        df_in.ecc_x = [_flyby_guess.ecc_i[0]]
        df_in.ecc_y = [_flyby_guess.ecc_i[1]]
        df_in.ecc_z = [_flyby_guess.ecc_i[2]]
        df_in.sma = [_flyby_guess.sma_i.to(u.km)]

        df_in.argp = [_flyby_guess.argp]
        df_in.inc = [_flyby_guess.inc]
        df_in.raan = [_flyby_guess.raan]

        df_in.v_p_x = [_flyby_guess.v_p_i[0]]
        df_in.v_p_y = [_flyby_guess.v_p_i[1]]
        df_in.v_p_z = [_flyby_guess.v_p_i[2]]
        df_in.r_p = [_flyby_guess.r_p.to(u.km)]
        df_in.dv = [_flyby_guess.r_p_dv * u.km/u.s]

        df_in.r_soi_x = [_flyby_guess.r_entry[0]]
        df_in.r_soi_y = [_flyby_guess.r_entry[1]]
        df_in.r_soi_z = [_flyby_guess.r_entry[2]]

        df_in.t_days = [_flyby_guess.t_i/3600/24]

        df_out.IN_OUT = ['OUT']
        df_out.v_x = [_flyby_guess.v_f[0]]
        df_out.v_y = [_flyby_guess.v_f[1]]
        df_out.v_z = [_flyby_guess.v_f[2]]
        df_out.v_inf_x = [_flyby_guess.v_inf_f[0]]
        df_out.v_inf_y = [_flyby_guess.v_inf_f[1]]
        df_out.v_inf_z = [_flyby_guess.v_inf_f[2]]
        df_out.v_planet_x = [_flyby_guess.v_planet_f[0]]
        df_out.v_planet_y = [_flyby_guess.v_planet_f[1]]
        df_out.v_planet_z = [_flyby_guess.v_planet_f[2]]
        df_out.ecc = [np.linalg.norm(_flyby_guess.ecc_f)]
        df_out.ecc_x = [_flyby_guess.ecc_f[0]]
        df_out.ecc_y = [_flyby_guess.ecc_f[1]]
        df_out.ecc_z = [_flyby_guess.ecc_f[2]]
        df_out.sma = [_flyby_guess.sma_f.to(u.km)]

        df_out.argp = [_flyby_guess.argp]
        df_out.inc = [_flyby_guess.inc]
        df_out.raan = [_flyby_guess.raan]

        df_out.v_p_x = [_flyby_guess.v_p_f[0]]
        df_out.v_p_y = [_flyby_guess.v_p_f[1]]
        df_out.v_p_z = [_flyby_guess.v_p_f[2]]
        df_out.r_p = [_flyby_guess.r_p.to(u.km)]
        df_out.dv = [_flyby_guess.r_p_dv * u.km/u.s]

        df_out.r_soi_x = [_flyby_guess.r_exit[0]]
        df_out.r_soi_y = [_flyby_guess.r_exit[1]]
        df_out.r_soi_z = [_flyby_guess.r_exit[2]]

        df_out.t_days = [_flyby_guess.t_f/3600/24]

    else:
        df_in.IN_OUT = ['IN']
        df_in.v_x = [_flyby_guess.v_i[0]]
        df_in.v_y = [_flyby_guess.v_i[1]]
        df_in.v_z = [_flyby_guess.v_i[2]]
        df_in.v_inf_x = [_flyby_guess.v_inf_i[0]]
        df_in.v_inf_y = [_flyby_guess.v_inf_i[1]]
        df_in.v_inf_z = [_flyby_guess.v_inf_i[2]]
        df_in.v_planet_x = [_flyby_guess.v_planet_i[0]]
        df_in.v_planet_y = [_flyby_guess.v_planet_i[1]]
        df_in.v_planet_z = [_flyby_guess.v_planet_i[2]]
        df_in.ecc = [np.linalg.norm(_flyby_guess.ecc_i)]
        df_in.ecc_x = [_flyby_guess.ecc_i[0]]
        df_in.ecc_y = [_flyby_guess.ecc_i[1]]
        df_in.ecc_z = [_flyby_guess.ecc_i[2]]
        df_in.sma = [_flyby_guess.sma_i.to(u.km).value]

        df_in.argp = [_flyby_guess.argp]
        df_in.inc = [_flyby_guess.inc]
        df_in.raan = [_flyby_guess.raan]

        df_in.v_p_x = [_flyby_guess.v_p_i[0].value]
        df_in.v_p_y = [_flyby_guess.v_p_i[1].value]
        df_in.v_p_z = [_flyby_guess.v_p_i[2].value]
        df_in.r_p = [_flyby_guess.r_p.to(u.km).value]
        df_in.dv = [_flyby_guess.r_p_dv]

        df_in.r_soi_x = [_flyby_guess.r_entry[0]]
        df_in.r_soi_y = [_flyby_guess.r_entry[1]]
        df_in.r_soi_z = [_flyby_guess.r_entry[2]]

        df_in.t_days = [_flyby_guess.t_i/3600/24]

        df_out.IN_OUT = ['OUT']
        df_out.v_x = [_flyby_guess.v_f[0]]
        df_out.v_y = [_flyby_guess.v_f[1]]
        df_out.v_z = [_flyby_guess.v_f[2]]
        df_out.v_inf_x = [_flyby_guess.v_inf_f[0]]
        df_out.v_inf_y = [_flyby_guess.v_inf_f[1]]
        df_out.v_inf_z = [_flyby_guess.v_inf_f[2]]
        df_out.v_planet_x = [_flyby_guess.v_planet_f[0]]
        df_out.v_planet_y = [_flyby_guess.v_planet_f[1]]
        df_out.v_planet_z = [_flyby_guess.v_planet_f[2]]
        df_out.ecc = [np.linalg.norm(_flyby_guess.ecc_f)]
        df_out.ecc_x = [_flyby_guess.ecc_f[0]]
        df_out.ecc_y = [_flyby_guess.ecc_f[1]]
        df_out.ecc_z = [_flyby_guess.ecc_f[2]]
        df_out.sma = [_flyby_guess.sma_f.to(u.km).value]

        df_out.argp = [_flyby_guess.argp]
        df_out.inc = [_flyby_guess.inc]
        df_out.raan = [_flyby_guess.raan]

        df_out.v_p_x = [_flyby_guess.v_p_f[0].value]
        df_out.v_p_y = [_flyby_guess.v_p_f[1].value]
        df_out.v_p_z = [_flyby_guess.v_p_f[2].value]
        df_out.r_p = [_flyby_guess.r_p.to(u.km).value]
        df_out.dv = [_flyby_guess.r_p_dv]

        df_out.r_soi_x = [_flyby_guess.r_exit[0]]
        df_out.r_soi_y = [_flyby_guess.r_exit[1]]
        df_out.r_soi_z = [_flyby_guess.r_exit[2]]

        df_out.t_days = [_flyby_guess.t_f/3600/24]
    return df_in.append(df_out)


def refined_flyby_df(_flyby_guess: flyby_refined, unit=False):
    columns = ['IN_OUT', 'v_x', 'v_y', 'v_z', 'v_inf_x', 'v_inf_y', 'v_inf_z', 'v_planet_x', 'v_planet_y', 'v_planet_z',
               'ecc', 'ecc_x', 'ecc_y', 'ecc_z', 'sma', 'argp', 'inc', 'raan', 't_days', 'r_soi_x', 'r_soi_y', 'r_soi_z',
               'v_p_x', 'v_p_y', 'v_p_z', 'r_p', 'dv', 'error_v', 'error_p']
    df_in = pd.DataFrame(columns=columns)
    df_out = pd.DataFrame(columns=columns)

    if unit:
        df_in.IN_OUT = ['IN']
        df_in.v_x = [_flyby_guess.v_i[0].to(u.km/u.s)]
        df_in.v_y = [_flyby_guess.v_i[1].to(u.km/u.s)]
        df_in.v_z = [_flyby_guess.v_i[2].to(u.km/u.s)]
        df_in.v_inf_x = [_flyby_guess.v_inf_i[0].to(u.km/u.s)]
        df_in.v_inf_y = [_flyby_guess.v_inf_i[1].to(u.km/u.s)]
        df_in.v_inf_z = [_flyby_guess.v_inf_i[2].to(u.km/u.s)]
        df_in.v_planet_x = [_flyby_guess.v_planet_i[0].to(u.km/u.s)]
        df_in.v_planet_y = [_flyby_guess.v_planet_i[1].to(u.km/u.s)]
        df_in.v_planet_z = [_flyby_guess.v_planet_i[2].to(u.km/u.s)]
        df_in.ecc = [np.linalg.norm(_flyby_guess.ecc_i)]
        df_in.ecc_x = [_flyby_guess.ecc_i[0]]
        df_in.ecc_y = [_flyby_guess.ecc_i[1]]
        df_in.ecc_z = [_flyby_guess.ecc_i[2]]
        df_in.sma = [_flyby_guess.sma_i.to(u.km)]

        df_in.argp = [_flyby_guess.argp]
        df_in.inc = [_flyby_guess.inc]
        df_in.raan = [_flyby_guess.raan]

        df_in.v_p_x = [_flyby_guess.v_p_i[0]]
        df_in.v_p_y = [_flyby_guess.v_p_i[1]]
        df_in.v_p_z = [_flyby_guess.v_p_i[2]]
        df_in.r_p = [_flyby_guess.r_p.to(u.km)]
        df_in.dv = [_flyby_guess.r_p_dv * u.km/u.s]

        df_in.r_soi_x = [_flyby_guess.r_entry[0]]
        df_in.r_soi_y = [_flyby_guess.r_entry[1]]
        df_in.r_soi_z = [_flyby_guess.r_entry[2]]

        df_in.t_days = [_flyby_guess.t_i/3600/24]

        df_in.error_v = [_flyby_guess.error_v]
        df_in.error_p = [_flyby_guess.error_p]

        df_out.IN_OUT = ['OUT']
        df_out.v_x = [_flyby_guess.v_f[0].to(u.km/u.s)]
        df_out.v_y = [_flyby_guess.v_f[1].to(u.km/u.s)]
        df_out.v_z = [_flyby_guess.v_f[2].to(u.km/u.s)]
        df_out.v_inf_x = [_flyby_guess.v_inf_f[0].to(u.km/u.s)]
        df_out.v_inf_y = [_flyby_guess.v_inf_f[1].to(u.km/u.s)]
        df_out.v_inf_z = [_flyby_guess.v_inf_f[2].to(u.km/u.s)]
        df_out.v_planet_x = [_flyby_guess.v_planet_f[0].to(u.km/u.s)]
        df_out.v_planet_y = [_flyby_guess.v_planet_f[1].to(u.km/u.s)]
        df_out.v_planet_z = [_flyby_guess.v_planet_f[2].to(u.km/u.s)]
        df_out.ecc = [np.linalg.norm(_flyby_guess.ecc_f)]
        df_out.ecc_x = [_flyby_guess.ecc_f[0]]
        df_out.ecc_y = [_flyby_guess.ecc_f[1]]
        df_out.ecc_z = [_flyby_guess.ecc_f[2]]
        df_out.sma = [_flyby_guess.sma_f.to(u.km)]

        df_out.argp = [_flyby_guess.argp]
        df_out.inc = [_flyby_guess.inc]
        df_out.raan = [_flyby_guess.raan]

        df_out.v_p_x = [_flyby_guess.v_p_f[0]]
        df_out.v_p_y = [_flyby_guess.v_p_f[1]]
        df_out.v_p_z = [_flyby_guess.v_p_f[2]]
        df_out.r_p = [_flyby_guess.r_p.to(u.km)]
        df_out.dv = [_flyby_guess.r_p_dv * u.km/u.s]

        df_out.r_soi_x = [_flyby_guess.r_exit[0]]
        df_out.r_soi_y = [_flyby_guess.r_exit[1]]
        df_out.r_soi_z = [_flyby_guess.r_exit[2]]

        df_out.t_days = [_flyby_guess.t_f/3600/24]

        df_out.error_v = [_flyby_guess.error_v]
        df_out.error_p = [_flyby_guess.error_p]

    else:
        df_in.IN_OUT = ['IN']
        df_in.v_x = [_flyby_guess.v_i[0].to(u.km/u.s)]
        df_in.v_y = [_flyby_guess.v_i[1].to(u.km/u.s)]
        df_in.v_z = [_flyby_guess.v_i[2].to(u.km/u.s)]
        df_in.v_inf_x = [_flyby_guess.v_inf_i[0].to(u.km/u.s)]
        df_in.v_inf_y = [_flyby_guess.v_inf_i[1].to(u.km/u.s)]
        df_in.v_inf_z = [_flyby_guess.v_inf_i[2].to(u.km/u.s)]
        df_in.v_planet_x = [_flyby_guess.v_planet_i[0].to(u.km/u.s)]
        df_in.v_planet_y = [_flyby_guess.v_planet_i[1].to(u.km/u.s)]
        df_in.v_planet_z = [_flyby_guess.v_planet_i[2].to(u.km/u.s)]
        df_in.ecc = [np.linalg.norm(_flyby_guess.ecc_i)]
        df_in.ecc_x = [_flyby_guess.ecc_i[0]]
        df_in.ecc_y = [_flyby_guess.ecc_i[1]]
        df_in.ecc_z = [_flyby_guess.ecc_i[2]]
        df_in.sma = [_flyby_guess.sma_i.to(u.km).value]

        df_in.argp = [_flyby_guess.argp]
        df_in.inc = [_flyby_guess.inc]
        df_in.raan = [_flyby_guess.raan]

        df_in.v_p_x = [_flyby_guess.v_p_i[0].value]
        df_in.v_p_y = [_flyby_guess.v_p_i[1].value]
        df_in.v_p_z = [_flyby_guess.v_p_i[2].value]
        df_in.r_p = [_flyby_guess.r_p.to(u.km).value]
        df_in.dv = [_flyby_guess.r_p_dv]

        df_in.r_soi_x = [_flyby_guess.r_entry[0]]
        df_in.r_soi_y = [_flyby_guess.r_entry[1]]
        df_in.r_soi_z = [_flyby_guess.r_entry[2]]

        df_in.t_days = [_flyby_guess.t_i/3600/24]

        df_in.error_v = [_flyby_guess.error_v]
        df_in.error_p = [_flyby_guess.error_p]

        df_out.IN_OUT = ['OUT']
        df_out.v_x = [_flyby_guess.v_f[0].to(u.km/u.s)]
        df_out.v_y = [_flyby_guess.v_f[1].to(u.km/u.s)]
        df_out.v_z = [_flyby_guess.v_f[2].to(u.km/u.s)]
        df_out.v_inf_x = [_flyby_guess.v_inf_f[0].to(u.km/u.s)]
        df_out.v_inf_y = [_flyby_guess.v_inf_f[1].to(u.km/u.s)]
        df_out.v_inf_z = [_flyby_guess.v_inf_f[2].to(u.km/u.s)]
        df_out.v_planet_x = [_flyby_guess.v_planet_f[0].to(u.km/u.s)]
        df_out.v_planet_y = [_flyby_guess.v_planet_f[1].to(u.km/u.s)]
        df_out.v_planet_z = [_flyby_guess.v_planet_f[2].to(u.km/u.s)]
        df_out.ecc = [np.linalg.norm(_flyby_guess.ecc_f)]
        df_out.ecc_x = [_flyby_guess.ecc_f[0]]
        df_out.ecc_y = [_flyby_guess.ecc_f[1]]
        df_out.ecc_z = [_flyby_guess.ecc_f[2]]
        df_out.sma = [_flyby_guess.sma_f.to(u.km).value]

        df_out.argp = [_flyby_guess.argp]
        df_out.inc = [_flyby_guess.inc]
        df_out.raan = [_flyby_guess.raan]

        df_out.v_p_x = [_flyby_guess.v_p_f[0].value]
        df_out.v_p_y = [_flyby_guess.v_p_f[1].value]
        df_out.v_p_z = [_flyby_guess.v_p_f[2].value]
        df_out.r_p = [_flyby_guess.r_p.to(u.km).value]
        df_out.dv = [_flyby_guess.r_p_dv]

        df_out.r_soi_x = [_flyby_guess.r_exit[0]]
        df_out.r_soi_y = [_flyby_guess.r_exit[1]]
        df_out.r_soi_z = [_flyby_guess.r_exit[2]]

        df_out.t_days = [_flyby_guess.t_f/3600/24]

        df_out.error_v = [_flyby_guess.error_v]
        df_out.error_p = [_flyby_guess.error_p]
    return df_in.append(df_out)


eq_hyperbola = Eq(np.square(var('x'))/np.square(var('a')) - np.square(var('y'))/np.square(var('b')), 1)


def pga_scalar_2_vector(v_inf_i, v_inf_f, v_p_i, v_p_f, a_i, a_f, e_i, e_f, r_p, body, rsoi, planetary_node):
    I = np.array([1, 0, 0])
    J = np.array([0, 1, 0])
    K = np.array([0, 0, 1])

    mu = body.k.to(u.km ** 3 / u.s ** 2).value
    e_i = np.linalg.norm(e_i)
    e_f = np.linalg.norm(e_f)
    # Orbital plane.
    n_vec_orbital = unit_vector(np.cross(v_inf_i, v_inf_f))

    # Rotation about orbital plane normal vector for v_p_unit_vec with angle of d_i.
    alpha_i = np.arcsin(1 / np.linalg.norm(e_i))
    v_p_unit_vec = unit_vector(np.dot(rotation_matrix(axis=n_vec_orbital,
                                                      theta=alpha_i), v_inf_i))

    # v_p_i_vec and v_p_f_vec
    v_p_i_vec = v_p_i * v_p_unit_vec
    v_p_f_vec = v_p_f * v_p_unit_vec

    # r_p_unit_vec and r_p_vec
    r_p_unit_vec = unit_vector(np.dot(rotation_matrix(axis=n_vec_orbital,
                                                      theta=-np.pi / 2), v_p_i_vec).value)
    r_p_vec = r_p * r_p_unit_vec

    print('r_p_vec:')
    print(r_p_vec)

    # FROM ORBITAL MECHANICS FOR ENGINEERS #################################################################
    # eccentricity vectors
    e_i_vec = np.cross(v_p_i_vec, np.cross(r_p_vec, v_p_i_vec)) / \
              body.k.to(u.km ** 3 / u.s ** 2).value - r_p_unit_vec

    assert abs(np.linalg.norm(e_i_vec) - e_i)/e_i <= 0.001

    e_f_vec = np.cross(v_p_f_vec, np.cross(r_p_vec, v_p_f_vec)) / \
              body.k.to(u.km ** 3 / u.s ** 2).value - r_p_unit_vec

    assert abs(np.linalg.norm(e_f_vec) - e_f)/e_f <= 0.001

    # Classical orbit parameters
    inclination = np.arccos(
        np.dot(K, n_vec_orbital) / (np.linalg.norm(n_vec_orbital)))

    n_vec = np.cross(K, np.cross(r_p_vec, v_p_i_vec))
    lan = np.arccos(np.dot(I, n_vec) / (np.linalg.norm(n_vec) * 1))

    if n_vec[1] < 0:
        lan = 2 * np.pi - lan

    aop = np.arccos(np.dot(n_vec, e_i_vec) / np.linalg.norm(n_vec) / np.linalg.norm(e_i_vec))

    if e_i_vec[-1] < 0:
        aop = 2 * np.pi - aop

    theta_inf_i = np.arccos(-1 / e_i)
    theta_inf_f = np.arccos(-1 / e_f)

    # def f_t_i(adjustment):
    #     th = _theta_rsoi(r_p, adjustment*rsoi, e_i)
    #     # r_entry = Orbit.from_classical(attractor=body, a=a_i, ecc=e_i * u.one, inc=inclination * u.rad,
    #     #                                raan=lan * u.rad, argp=aop * u.rad, nu=0*u.rad, epoch=time.Time(planetary_node.epoch_periapsis, scale='tdb'))
    #     #
    #
    # def f_t_f(adjustment):
    #     th = _theta_rsoi(r_p, adjustment*rsoi, e_f)
    #


    # r_exit = Orbit.from_classical(attractor=body, a=a_f, ecc=e_f * u.one, inc=inclination * u.rad,
    #                               raan=lan * u.rad, argp=aop * u.rad, nu=0*u.rad, epoch=time.Time(planetary_node.epoch_periapsis,scale='tdb'))
    #
    # r_entry = r_entry.propagate(time.TimeDelta(tp_i*u.s)).r
    # r_exit = r_exit.propagate(time.TimeDelta(tp_f*u.s)).r

    # adj_i = optimize.newton(f_t_i, 0)
    # adj_f = optimize.newton(f_t_f, 0)

    theta_rsoi_i = _theta_rsoi(r_p, 1.051*rsoi, e_i)
    theta_rsoi_f = _theta_rsoi(r_p, 1.05256*rsoi, e_f)

    print('true_anom_entry')
    print(theta_rsoi_i)

    print('\ntrue_anom_exit')
    print(theta_rsoi_f)

    if n_vec_orbital[-1] > 0:
        theta_inf_i = -theta_inf_i
        theta_rsoi_i = -theta_rsoi_i

    else:
        theta_inf_f = -theta_inf_f
        theta_rsoi_f = -theta_rsoi_f

    F_i = F(e_i, theta_rsoi_i.value)
    F_f = F(e_f, theta_rsoi_f.value)

    Mh_i = Mn(e_i, F_i)
    Mh_f = Mn(e_f, F_f)

    h_i = np.sqrt(2*planetary_node.body.k.to(u.km**3/u.s**2).value * r_p.value)
    h_f = np.sqrt(2*planetary_node.body.k.to(u.km**3/u.s**2).value * r_p.value)

    tp_i = tp(Mh_i, h_i, e_i, planetary_node.body.k.to(u.km**3/u.s**2).value)
    tp_f = tp(Mh_f, h_f, e_f, planetary_node.body.k.to(u.km**3/u.s**2).value)

    r_entry = Orbit.from_classical(attractor=body, a=a_i, ecc=e_i * u.one, inc=inclination * u.rad,
                                   raan=lan * u.rad, argp=aop * u.rad, nu=0*u.rad, epoch=time.Time(planetary_node.epoch_periapsis, scale='tdb'))

    r_exit = Orbit.from_classical(attractor=body, a=a_f, ecc=e_f * u.one, inc=inclination * u.rad,
                                  raan=lan * u.rad, argp=aop * u.rad, nu=0*u.rad, epoch=time.Time(planetary_node.epoch_periapsis,scale='tdb'))

    r_entry = r_entry.propagate(time.TimeDelta(tp_i*u.s)).r
    r_exit = r_exit.propagate(time.TimeDelta(tp_f*u.s)).r
    print(np.linalg.norm(r_exit))
    print(np.linalg.norm(r_entry))
    print(rsoi)
    return v_p_i_vec, v_p_f_vec, lan, aop, inclination, tp_i, tp_f, e_i_vec, e_f_vec, r_entry, r_exit

