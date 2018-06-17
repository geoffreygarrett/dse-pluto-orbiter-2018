from .planetary_node import PlanetaryNode
from .tools.external_reference import *
from .tools.geometry import *
from poliastro.twobody import Orbit

class Flyby(object):
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

    @property
    def planetary_node(self):
        return self._planetary_node

    @property
    def v_i(self):
        return self._v_inf_i

    @property
    def v_f(self):
        return self._v_inf_f

    @v_i.setter
    def v_i(self, arg):
        self._v_i = arg

    @v_f.setter
    def v_f(self, arg):
        self._v_i = arg

    @property
    def state(self):
        return self._state

    @property
    def _state(self):
        return self.__state[-1]

    @_state.setter
    def _state(self, arg):
        _state_error = "Maintain the order of states: \n {}".format(self._order_of_states)
        try:
            if len(self.__state) is 0:
                assert arg == self._order_of_states[0], _state_error
                self.__state = []
            elif len(self.__state) is 1:
                assert arg == self._order_of_states[1], _state_error
            elif len(self.__state) is 2:
                assert arg == self._order_of_states[2], _state_error
        except AssertionError:
            raise EnvironmentError(_state_error)
        self.__state.append(arg)

    def _base_gravity_assist(self, v_i, v_f):
        self.v_i = v_i
        self.v_f = v_f
        v_planet_periapsis = Orbit.from_body_ephem(self.planetary_node.body, self.planetary_node.epoch_periapsis)
        self._v_inf_i = self.v_i - v_planet_periapsis
        self._v_inf_f = self.v_f - v_planet_periapsis
        self._alpha_required = angle_between(self._v_inf_i, self._v_inf_f)

    def check_gravity_assist(self, v_i, v_f):
        self._base_gravity_assist(v_i, v_f)
        self._state = 'checked'
        alpha_max = alpha(self._v_inf_i, self._v_inf_f, self.planetary_node.periapsis_minimum,
                          self.planetary_node.body.k)
        try:
            assert self._alpha_required <= alpha_max
        except AssertionError:
            raise NotImplementedError("Gravity assist bending angle required (α_req) exceeds possible bending angle "
                                      "(α_max)\nwith a single impulse at periapsis. Multiple impulses are not "
                                      "implemented.\n{} > {}".format(self._alpha_required, alpha_max))

    def rough_powered_gravity_assist(self, v_i=None, v_f=None):
        """
        :param v_i: [x, y, z] * astropy.unit
        :param v_f: [x, y, z] * astropy.unit
        :param mode: (str) | check / scalar / vector
        :return:
        """
        if self._alpha_required is None:
            assert v_i is not None and v_f is not None, "v_i and v_f must be provided in this current state. \n" \
                                                        "Current state: {}".format(self.state)
            self.check_gravity_assist(v_i, v_f)
        self._state= 'rough'
        # TODO: Add first epoch_entry and epoch_exit for rough solution.
        self._ecc_i, self._ecc_f, self._sma_i, self._sma_f, self._vp_i, self._vp_f, self._rp_DV , self._rp = \
            newton_rhapson_pga(self._v_inf_i, self._v_inf_f, self.planetary_node.body.k, self._alpha_required)

    def refined_powered_gravity_assist(self, v_i=None, v_f=None):





