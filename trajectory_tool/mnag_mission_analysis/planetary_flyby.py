from .planetary_node import PlanetaryNode
from .tools.external_reference import *

class Flyby(object):
    def __init__(self, planetary_node):
        self._order_of_states = ['checked', 'rough', 'refined']

        self._body = planetary_node.body
        self._epoch_rp = planetary_node.epoch_rp
        self._state = None

        self._v_epoch_i = None
        self._v_i = None
        self._v_inf_i = None
        self._v_planet_i = None
        self._ecc_i = None
        self._sma_i = None

        self._vp_i = None
        self._rp = None
        self._rp_DV = None
        self._vp_f = None

        self._v_epoch_f = None
        self._v_f = None
        self._v_inf_f = None
        self._v_planet_f = None
        self._ecc_f = None
        self._sma_f = None

        self._aop = None
        self._inc = None
        self._raan = None

    @property
    def state(self):
        return self._state[-1]

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

    @state.setter
    def state(self, arg):
        _state_error = "Maintain the order of states: \n {}".format(self._order_of_states)
        if len(self._state) is 0:
            assert arg == self._order_of_states[0], _state_error
        elif len(self._state) is 1:
            assert arg == self._order_of_states[1], _state_error
        elif len(self._state) is 2:
            assert arg == self._order_of_states[2], _state_error
        self._state.append(arg)

    def _base_gravity_assist(self, v_i, v_f):
        self.v_i = v_i
        self.v_f = v_f

    def check_gravity_assist(self, v_i, v_f):
        self.state = 'checked'

    def rough_powered_gravity_assist(self, v_i, v_f):
        """

        :param v_i: [x, y, z] * astropy.unit
        :param v_f: [x, y, z] * astropy.unit
        :param mode: (str) | check / scalar / vector
        :return:
        """
        self.state= 'rough'


    def refined_powered_gravity_assist(self, v_inf_i, v_inf_f):


