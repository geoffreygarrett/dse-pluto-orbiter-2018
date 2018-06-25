from .planetary_node import PlanetaryNode
from .tools.flyby_helper import *
from .config import *
from .tools.geometry import *
from .tools.plotting import *
from poliastro.twobody import Orbit
from astropy import time
from .data_structures import *
import pandas as pd
# from trajectory_tool.genetic_algorithim_analysis.genetic_algorithim import Chromosome


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
        self.count = 0

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
        self._refined_attributes = arg

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
        self._guess_dataframe = guess_flyby_df(arg, unit=True)

    @refined_dataframe.setter
    def refined_dataframe(self, arg: flyby_refined):
        self._refined_dataframe = refined_flyby_df(arg, unit=True)

    @state.setter
    def state(self, arg):
        self._state = arg

    def plot(self):
        plot_planetary_flyby(self)

    def _base_gravity_assist(self, v_i, v_f):

        # TODO: TRANSFORM V_INF_I TO GCRS
        v_inf_i = v_i - self.planetary_node.v_planet_i
        v_inf_f = v_f - self.planetary_node.v_planet_f
        a_req = angle_between(v_inf_i, v_inf_f)
        a_max = alpha(v_inf_i, v_inf_f, self.planetary_node.periapsis_minimum, self.planetary_node.body.k)
        try:
            assert a_req * u.rad <= a_max
        except AssertionError:
            raise NotImplementedError("Gravity assist bending angle required (α_req) exceeds possible bending angle "
                                      "(α_max)\nwith a single impulse at periapsis. Multiple impulses are not "
                                      "implemented.\n{} > {}".format(a_req, a_max))
        return v_i, v_f, self.planetary_node.v_planet_i, self.planetary_node.v_planet_f, v_inf_i, v_inf_f, a_req

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

    def _guess_powered_gravity_assist(self, v_i, v_f):
        if self.state is 'checked':
            pass
        else:
            self.check_gravity_assist(v_i, v_f)
        # print(guess_flyby(self._basic_powered_gravity_assist(v_i, v_f), self.planetary_node))
        return guess_flyby(self._basic_powered_gravity_assist(v_i, v_f), self.planetary_node)

    def guess_powered_gravity_assist(self, v_i, v_f):
        self.guess_attributes = self._guess_powered_gravity_assist(v_i, v_f)

    def refine_powered_gravity_assist(self, v_i, v_f):
        if self.state is 'checked':
            pass
        else:
            self.check_gravity_assist(v_i, v_f)

        self.refined_attributes = refine_flyby(self._base_gravity_assist(v_i, v_f), self.planetary_node)

        v_planet_i_old = self.planetary_node.v_planet_i
        v_planet_f_old = self.planetary_node.v_planet_f

        self.planetary_node.epoch_entry = self.planetary_node.epoch_periapsis + timedelta(seconds=self._refined_attributes.t_i)
        self.planetary_node.epoch_exit = self.planetary_node.epoch_periapsis + timedelta(seconds=self._refined_attributes.t_f)

        self._refined_attributes.error_v = np.linalg.norm(self.planetary_node.v_planet_i.to(u.km/u.s) - v_planet_i_old.to(u.km/u.s)) + \
                                           np.linalg.norm(self.planetary_node.v_planet_f.to(u.km/u.s) - v_planet_f_old.to(u.km/u.s))

        self._refined_attributes.error_p = np.linalg.norm(self.planetary_node.soi_entry_position_body_ecliptic - self._refined_attributes.r_entry) + \
                                           np.linalg.norm(self.planetary_node.soi_exit_position_body_ecliptic - self._refined_attributes.r_exit)

        self.planetary_node.soi_entry_position_body_ecliptic = self._refined_attributes.r_entry
        self.planetary_node.soi_exit_position_body_ecliptic = self._refined_attributes.r_exit







