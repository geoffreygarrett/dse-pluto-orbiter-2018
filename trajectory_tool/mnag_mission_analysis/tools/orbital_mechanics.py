from poliastro.twobody import Orbit
import numpy as np
from astropy import time


def soi(body, epoch):
    parent_body = body.parent
    mu_small = body.k
    mu_large = parent_body.k
    r_between = (Orbit.from_body_ephem(parent_body, time.Time(epoch, scale='tdb')).r -
                 Orbit.from_body_ephem(body, time.Time(epoch, scale='tdb')).r)
    r_value, r_unit = np.linalg.norm(r_between), r_between.unit
    return r_value * (mu_small/mu_large) ** (2/5) * r_unit


