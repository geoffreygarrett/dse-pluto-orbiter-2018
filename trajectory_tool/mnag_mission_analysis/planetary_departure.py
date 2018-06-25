from .planetary_node import PlanetaryNode
from trajectory_tool.mnag_mission_analysis.departure_helper import *
from .planetary_node import PlanetaryNode
import numpy as np
import astropy.units as u


class PlanetaryDeparture(object):

    def __init__(self, planetary_node: PlanetaryNode):
        self.planetary_node = planetary_node
        self._body = planetary_node.body
        self._epoch_periapsis = planetary_node.epoch_periapsis
        self._parking_rp = 180 + 6371
    # def calculate_refined(self, v_exit, r_p):
    #     return

    @property
    def v_inf_f(self):
        return np.linalg.norm(self.planetary_node.v_exit - self.planetary_node.v_planet_f)

    @property
    def v_p_hyp(self):
        return np.sqrt(np.square(self.v_inf_f) + (2* self.planetary_node.body.k.to(u.km**3/u.s**2).value)/self._parking_rp)

    @property
    def v_p_cir(self):
        return np.sqrt(self.planetary_node.body.k.to(u.km**3/u.s**2).value/ self._parking_rp)

    @property
    def delta_v(self):
        return np.linalg.norm(self.v_p_hyp-self.v_p_cir)

