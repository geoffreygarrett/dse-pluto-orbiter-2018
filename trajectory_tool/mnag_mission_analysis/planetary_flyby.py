from .planetary_node import PlanetaryNode
from .tools.flyby_helper import *
from .config import *
from .tools.geometry import *
from .tools.plotting import *
from poliastro.twobody import Orbit
from astropy import time
from .data_structures import *
import pandas as pd


class PlanetaryFlyby(object):
    def __init__(self, planetary_node: PlanetaryNode):
        # General Attributes
        self._planetary_node = planetary_node
        self._body = planetary_node.body
        self._epoch_periapsis = planetary_node.epoch_periapsis
        self._state = None

        # Attribute data structures
        self._basic_attributes = None
        self._guess_attributes = None
        self._refined_attributes = None

        # Attribute data frames
        self._basic_dataframe = None
        self._guess_dataframe = None
        self._refined_dataframe = None

    def __repr__(self):
        return str(self._body) + ' Planetary Flyby | Periapsis epoch: ' + str(self._epoch_periapsis)

    @property
    def planetary_node(self):
        return self._planetary_node

    @property
    def state(self):
        return self._state

    @property
    def guess_attributes(self):
        return self._guess_attributes

    @property
    def basic_attributes(self):
        return self._basic_attributes

    @property
    def refined_attributes(self):
        if self._refined_attributes is None:
            return self._guess_attributes
        else:
            return self._refined_attributes

    @basic_attributes.setter
    def basic_attributes(self, arg: flyby_basic):
        self.basic_dataframe = arg
        self._basic_attributes = arg

    @guess_attributes.setter
    def guess_attributes(self, arg: flyby_guess):
        self.guess_dataframe = arg
        self._guess_attributes = arg

    @refined_attributes.setter
    def refined_attributes(self, arg: flyby_refined):
        self.refined_dataframe = arg

    @property
    def basic_dataframe(self):
        return self._basic_dataframe

    @property
    def guess_dataframe(self):
        return self._guess_dataframe

    @property
    def refined_dataframe(self):
        return self._refined_dataframe

    @basic_dataframe.setter
    def basic_dataframe(self, arg: flyby_basic):
        self._basic_dataframe = basic_flyby_df(arg, unit=False)

    @guess_dataframe.setter
    def guess_dataframe(self, arg: flyby_guess):
        self._guess_dataframe = guess_flyby_df(arg, unit=False)

    @refined_dataframe.setter
    def refined_dataframe(self, arg: flyby_refined):
        self._refined_dataframe = refined_flyby_df(arg, unit=False)

    @state.setter
    def state(self, arg):
        self._state = arg

    def plot(self):
        plot_planetary_flyby(self)

    def _base_gravity_assist(self, v_i, v_f):
        v_planet_i = Orbit.from_body_ephem(self.planetary_node.body, time.Time(self.planetary_node.epoch_entry,
                                                                               scale='tdb')).state.v
        v_planet_f = Orbit.from_body_ephem(self.planetary_node.body, time.Time(self.planetary_node.epoch_exit,
                                                                               scale='tdb')).state.v
        v_inf_i = v_i - v_planet_i
        v_inf_f = v_f - v_planet_f
        a_req = angle_between(v_inf_i, v_inf_f)
        a_max = alpha(v_inf_i, v_inf_f, self.planetary_node.periapsis_minimum, self.planetary_node.body.k)
        try:
            assert a_req * u.rad <= a_max
        except AssertionError:
            raise NotImplementedError("Gravity assist bending angle required (α_req) exceeds possible bending angle "
                                      "(α_max)\nwith a single impulse at periapsis. Multiple impulses are not "
                                      "implemented.\n{} > {}".format(a_req, a_max))
        return v_i, v_f, v_planet_i, v_planet_f, v_inf_i, v_inf_f, a_req

    def check_gravity_assist(self, v_i, v_f):
        self.state = 'checked'
        self._base_gravity_assist(v_i, v_f)

    def _basic_powered_gravity_assist(self, v_i, v_f):
        if self.state is 'checked':
            pass
        else:
            self.check_gravity_assist(v_i, v_f)
        return basic_flyby(self._base_gravity_assist(v_i, v_f), self.planetary_node)

    def basic_powered_gravity_assist(self, v_i, v_f):
        self.basic_attributes = self._basic_powered_gravity_assist(v_i, v_f)

    def guess_powered_gravity_assist(self, v_i, v_f):
        if self.state is 'checked':
            pass
        else:
            self.check_gravity_assist(v_i, v_f)
        print(guess_flyby(self._basic_powered_gravity_assist(v_i, v_f), self.planetary_node))
        self.guess_attributes = guess_flyby(self._basic_powered_gravity_assist(v_i, v_f), self.planetary_node)

    def refine_powered_gravity_assist(self, v_i, v_f):
        if self.state is 'checked':
            pass
        else:
            self.check_gravity_assist(v_i, v_f)
        if self.guess_attributes is None:
            self._guess_attributes = guess_flyby(self._basic_powered_gravity_assist(v_i, v_f), self.planetary_node)
        # self._refined_attributes = refine_flyby(self._basic_powered_gravity_assist(v_i, v_f), self.planetary_node,
        #                                         self.refined_attributes)

        self._vp_i, self._vp_f, self._raan, self._aop, self._inc, self._t_i, self._t_f, self._ecc_i, self._ecc_f, \
        self.planetary_node.soi_entry_position_body_ecliptic, self.planetary_node.soi_exit_position_body_ecliptic \
            = \
            pga_scalar_2_vector(self._v_inf_i, self._v_inf_f, self._vp_i, self._vp_f, self._sma_i, self._sma_f,
                                self._ecc_i, self._ecc_f, self._rp, self._body,
                                self.planetary_node.soi_periapsis_magnitude, self._epoch_periapsis)
