from .planetary_node import PlanetaryNode
import numpy as np
import astropy.units as u


class PlanetaryRendezvous(object):
    def __init__(self, planetary_node: PlanetaryNode):
        self.planetary_node = planetary_node
        self._body = planetary_node.body
        self._epoch_periapsis = planetary_node.epoch_periapsis

        self._capture_rp = 1588
        self._capture_e = 0.25

    @property
    def v_inf_i(self):
        return np.linalg.norm(self.planetary_node.v_entry - self.planetary_node.v_planet_i)

    @property
    def v_p_hyp(self):
        v_inf_i = self.v_inf_i
        return np.sqrt(np.square(v_inf_i) + (2*self.planetary_node.body.k.to(u.km**3/u.s**2).value)/self._capture_rp)

    @property
    def v_p_capture(self):
        return np.sqrt(self.planetary_node.body.k.to(u.km**3/u.s**2).value*(1-self._capture_e)/self._capture_rp)

    @property
    def delta_v(self):
        return np.linalg.norm(self.v_p_hyp-self.v_p_capture)