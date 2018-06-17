from .planetary_node import PlanetaryNode
from .tools.external_reference import *
from .tools.geometry import *
from poliastro.twobody import Orbit
from astropy import time


class PlanetaryFlyby(object):
    def __init__(self, planetary_node:PlanetaryNode):
        self._order_of_states = ['checked', 'rough', 'refined', 'plot']

        self._planetary_node = planetary_node
        self._body = planetary_node.body
        self._epoch_periapsis = planetary_node.epoch_periapsis
        self.__state = None

        self._v_i = None
        self._v_inf_i = None
        self._v_planet_i = None
        self._ecc_i = None
        self._sma_i = None

        self._vp_i = None
        self._rp = None
        self._rp_DV = None
        self._vp_f = None

        self._v_f = None
        self._v_inf_f = None
        self._v_planet_f = None
        self._ecc_f = None
        self._sma_f = None

        self._aop = None
        self._inc = None
        self._raan = None

    def __repr__(self):
        return str(self._body) + ' Planetary Flyby | Periapsis epoch: ' + str(self._epoch_periapsis)

    @property
    def planetary_node(self):
        return self._planetary_node

    @property
    def v_i(self):
        return self._v_i

    @property
    def v_f(self):
        return self._v_f

    @v_i.setter
    def v_i(self, arg):
        self._v_i = arg

    @v_f.setter
    def v_f(self, arg):
        self._v_i = arg

    def _base_gravity_assist(self, v_i, v_f):
        self._v_i = v_i
        self._v_f = v_f
        v_planet_entry = Orbit.from_body_ephem(self.planetary_node.body,
                                               time.Time(self.planetary_node.epoch_entry, scale='tdb')).state.v
        v_planet_exit = Orbit.from_body_ephem(self.planetary_node.body,
                                              time.Time(self.planetary_node.epoch_exit, scale='tdb')).state.v
        self._v_inf_i = self._v_i - v_planet_entry
        self._v_inf_f = self._v_f - v_planet_exit
        self._alpha_required = angle_between(self._v_inf_i, self._v_inf_f)

    def check_gravity_assist(self):
        self._base_gravity_assist(self.planetary_node.v_entry, self.planetary_node.v_exit)
        alpha_max = alpha(self._v_inf_i, self._v_inf_f, self.planetary_node.periapsis_minimum,
                          self.planetary_node.body.k)
        try:
            assert self._alpha_required*u.rad <= alpha_max
        except AssertionError:
            raise NotImplementedError("Gravity assist bending angle required (α_req) exceeds possible bending angle "
                                      "(α_max)\nwith a single impulse at periapsis. Multiple impulses are not "
                                      "implemented.\n{} > {}".format(self._alpha_required, alpha_max))

    def rough_powered_gravity_assist(self):
        self.check_gravity_assist()
        # TODO: Add first epoch_entry and epoch_exit for rough solution.
        self._ecc_i, self._ecc_f, self._sma_i, self._sma_f, self._vp_i, self._vp_f, self._rp_DV, self._rp = \
            newton_rhapson_pga(self._v_inf_i, self._v_inf_f, self.planetary_node.body.k, self._alpha_required)

    # def refined_powered_gravity_assist(self):





