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
        self._order_of_states = ['checked', 'rough', 'refined']

        self._basic_attributes = flyby_basic
        self._guess_attributes = flyby_guess
        self._refined_attributes = flyby_refined

        # General Attributes
        self._planetary_node = planetary_node
        self._body = planetary_node.body
        self._epoch_periapsis = planetary_node.epoch_periapsis
        self._state = None

        # Basic Attributes (Entry)
        self._v_i = None
        self._v_inf_i = None
        self._v_planet_i = None
        self._ecc_i = None
        self._sma_i = None
        self._vp_i = None

        # Basic Attributes (Periapsis)
        self._rp = None
        self._rp_DV = None

        # Basic Attributes (Exit)
        self._vp_f = None
        self._v_f = None
        self._v_inf_f = None
        self._v_planet_f = None
        self._ecc_f = None
        self._sma_f = None

        # Refined Update Basic Attributes

        # Refined Attributes
        self._aop = None
        self._inc = None
        self._raan = None
        self._t_i = None
        self._t_f = None

        self._df_entry_rough = None
        self._df_exit_rough = None

        self._df_entry_refined = None
        self._df_exit_refined = None

        self._rough = None
        self._refined = None

    def __repr__(self):
        return str(self._body) + ' Planetary Flyby | Periapsis epoch: ' + str(self._epoch_periapsis)

    @property
    def data_frame(self):
        # TODO: Clean this up eventually, ONLY prints for rough currently.
        rough = ['ENTRY_EXIT', 'V_heliocentric_x', 'V_heliocentric_y', 'V_heliocentric_z', 'V_infinity_x',
                 'V_infinity_y', 'V_infinity_z', 'V_planet_x', 'V_planet_y', 'V_planet_z', 'V_periapsis', 'SMA', 'ECC',
                 'rp', 'DV']
        refined = []
        df = pd.DataFrame(columns=rough + refined)
        df = df.astype('object')
        df.ENTRY_EXIT = ['ENTRY']
        df.V_heliocentric_x = [self._v_i.to(u.km / u.s)[0]]
        df.V_heliocentric_y = [self._v_i.to(u.km / u.s)[1]]
        df.V_heliocentric_z = [self._v_i.to(u.km / u.s)[2]]
        df.V_infinity_x = [self._v_inf_i.to(u.km / u.s)[0]]
        df.V_infinity_y = [self._v_inf_i.to(u.km / u.s)[1]]
        df.V_infinity_z = [self._v_inf_i.to(u.km / u.s)[2]]
        df.V_periapsis = [self._vp_i * (u.km / u.s)]
        df.V_planet_x = [self._v_planet_i.to(u.km / u.s)[0]]
        df.V_planet_y = [self._v_planet_i.to(u.km / u.s)[1]]
        df.V_planet_z = [self._v_planet_i.to(u.km / u.s)[2]]
        df.SMA = [self._sma_i]
        df.ECC = [self._ecc_i]
        df.rp = [self._rp]
        df.DV = [self._rp_DV]

        df2 = pd.DataFrame(columns=rough + refined)
        df2 = df2.astype('object')
        df2.ENTRY_EXIT = ['EXIT']
        df2.V_heliocentric_x = [self._v_f.to(u.km / u.s)[0]]
        df2.V_heliocentric_y = [self._v_f.to(u.km / u.s)[1]]
        df2.V_heliocentric_z = [self._v_f.to(u.km / u.s)[2]]
        df2.V_infinity_x = [self._v_inf_f.to(u.km / u.s)[0]]
        df2.V_infinity_y = [self._v_inf_f.to(u.km / u.s)[1]]
        df2.V_infinity_z = [self._v_inf_f.to(u.km / u.s)[2]]
        df2.V_periapsis = [self._vp_f * (u.km / u.s)]
        df2.V_planet_x = [self._v_planet_f.to(u.km / u.s)[0]]
        df2.V_planet_y = [self._v_planet_f.to(u.km / u.s)[1]]
        df2.V_planet_z = [self._v_planet_f.to(u.km / u.s)[2]]
        df2.SMA = [self._sma_f]
        df2.SMA = df2.SMA
        df2.ECC = [self._ecc_f]
        df2.rp = self._rp
        df2.DV = self._rp_DV

        np.set_printoptions(precision=3)
        df_all = df.append(df2)
        for par in ['SMA', 'ECC', 'DV', 'rp', 'V_heliocentric_x', 'V_heliocentric_y', 'V_heliocentric_z',
                    'V_infinity_x', 'V_infinity_y', 'V_infinity_z', 'V_planet_x', 'V_planet_y', 'V_planet_z',
                    'V_periapsis']:
            df_all[par] = df_all[par].map(lambda x: '{0:.5}'.format(x))

        return df_all

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

    @property
    def state(self):
        return self._state

    @state.setter
    def state(self, arg):
        if arg not in ['checked', 'rough', 'refined']:
            raise EnvironmentError("Order of states not satisfied.")
        else:
            self._state = arg

    @property
    def rough(self):
        return self._rough

    @property
    def refined(self):
        return self._refined

    # @rough.setter
    # def rough(self, arg:list):
    #
    # @refined.setter
    # def refined(self, arg:list):

    def plot(self):
        plot_planetary_flyby(self)

    def _base_gravity_assist(self, v_i, v_f):
        self._v_i = v_i
        self._v_f = v_f
        self._v_planet_i = v_planet_entry = Orbit.from_body_ephem(self.planetary_node.body,
                                                                  time.Time(self.planetary_node.epoch_entry,
                                                                            scale='tdb')).state.v
        self._v_planet_f = v_planet_exit = Orbit.from_body_ephem(self.planetary_node.body,
                                                                 time.Time(self.planetary_node.epoch_exit,
                                                                           scale='tdb')).state.v
        self._v_inf_i = self._v_i - v_planet_entry
        self._v_inf_f = self._v_f - v_planet_exit
        self._alpha_required = angle_between(self._v_inf_i, self._v_inf_f)

        return v_i, v_f, v_planet_i, v_planet_f, v_inf_i, v_inf_f, a_req

    def check_gravity_assist(self):
        self.state = 'checked'
        self._base_gravity_assist(self.planetary_node.v_entry, self.planetary_node.v_exit)
        alpha_max = alpha(self._v_inf_i, self._v_inf_f, self.planetary_node.periapsis_minimum,
                          self.planetary_node.body.k)
        try:
            assert self._alpha_required * u.rad <= alpha_max
        except AssertionError:
            raise NotImplementedError("Gravity assist bending angle required (α_req) exceeds possible bending angle "
                                      "(α_max)\nwith a single impulse at periapsis. Multiple impulses are not "
                                      "implemented.\n{} > {}".format(self._alpha_required, alpha_max))

    def basic_powered_gravity_assist(self):
        self.state = 'rough'
        self.check_gravity_assist()
        # TODO: Add first epoch_entry and epoch_exit for rough solution.
        self._ecc_i, self._ecc_f, self._sma_i, self._sma_f, self._vp_i, self._vp_f, self._rp_DV, self._rp = \
            newton_rhapson_pga(self._v_inf_i, self._v_inf_f, self.planetary_node.body.k, self._alpha_required)

    def guess_powered_gravity_assist(self):
        self.state = 'refined'
        self.check_gravity_assist()
        self._vp_i, self._vp_f, self._raan, self._aop, self._inc, self._t_i, self._t_f, self._ecc_i, self._ecc_f, \
        self.planetary_node.soi_entry_position_body_ecliptic, self.planetary_node.soi_exit_position_body_ecliptic \
            = \
            pga_scalar_2_vector(self._v_inf_i, self._v_inf_f, self._vp_i, self._vp_f, self._sma_i, self._sma_f,
                                self._ecc_i, self._ecc_f, self._rp, self._body,
                                self.planetary_node.soi_periapsis_magnitude, self._epoch_periapsis)

    def refined_powered_gravity_assist(self):
        self.state = 'refined'
        self.check_gravity_assist()
        self._vp_i, self._vp_f, self._raan, self._aop, self._inc, self._t_i, self._t_f, self._ecc_i, self._ecc_f, \
        self.planetary_node.soi_entry_position_body_ecliptic, self.planetary_node.soi_exit_position_body_ecliptic \
            = \
            pga_scalar_2_vector(self._v_inf_i, self._v_inf_f, self._vp_i, self._vp_f, self._sma_i, self._sma_f,
                                self._ecc_i, self._ecc_f, self._rp, self._body,
                                self.planetary_node.soi_periapsis_magnitude, self._epoch_periapsis)













