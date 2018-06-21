import numpy as np


def _e(r_p, v_inf, mu):
    return 1 + (r_p * v_inf ** 2)/mu


def _h(r_p, v_inf, mu):
    return r_p * np.sqrt(np.square(v_inf) + 2 * (mu)/(r_p))


def _vp(h, r_p):
    return h/r_p


def _vc(mu, r_p):
    return np.sqrt(mu/r_p)


def _beta(r_p, v_inf, mu):
    return np.arccos(1/(1+(r_p * np.square(v_inf))/(mu)))


def refined_departure(v_exit, r_p, planetary_node):
    v_inf = v_exit - planetary_node.v_planet_i
    e = _e(r_p, v_inf, planetary_node.body.k.value)
    h = _h(r_p, v_inf, planetary_node.body.k.value)
    vp = _vp(h, r_p)
    vc = _vc(planetary_node.body.k.value, r_p)
    dv = vp - vc
    beta = _beta(r_p, v_inf, planetary_node.body.k.value)
