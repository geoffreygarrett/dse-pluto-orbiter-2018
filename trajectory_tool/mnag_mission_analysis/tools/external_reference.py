import numpy as np
import astropy.units as u
from scipy import optimize


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


def hyp_a(v_inf, mu):
    return (- mu / (np.square(np.linalg.norm(v_inf.to(u.km / u.s))) * (u.km / u.s) ** 2)).to(u.km)


def newton_rhapson_pga(v_inf_i, v_inf_f, mu_body, alpha_required):
    e_i_0 = 1.1
    a_i = hyp_a(v_inf_i, mu_body)
    a_f = hyp_a(v_inf_f, mu_body)

    alpha_required = alpha_required

    def func(_e_i):
        return f_e_i(_e_i, a_i=a_i.to(u.km).value, a_f=a_f.to(u.km).value,
                     alpha_req=alpha_required)

    def _fprime(_e_i):
        return d_f_e_i(_e_i, a_i=a_i.to(u.km).value, a_f=a_f.to(u.km).value)

    e_i = optimize.newton(func, e_i_0, _fprime)
    r_p = a_i * (1 - e_i)
    e_f = - (r_p / a_f - 1)
    v_p_i = np.sqrt(np.square(np.linalg.norm(v_inf_i)) * (e_i + 1) / (e_i - 1))
    v_p_f = np.sqrt(np.square(np.linalg.norm(v_inf_i)) * (e_f + 1) / (e_f - 1))
    dv = np.linalg.norm(v_p_f - v_p_i)
    return e_i, e_f, a_i, a_f, v_p_i, v_p_f, dv, r_p
