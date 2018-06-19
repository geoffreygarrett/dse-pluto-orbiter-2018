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
    def state(self):
        return self._state

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
        self._basic_attributes = self._base_gravity_assist(v_i, v_f)

    @property
    def guess_attributes(self):
        return self._guess_attributes

    @property
    def basic_attributes(self):
        return self._basic_attributes

    @property
    def refined_attributes(self):
        return self._refined_attributes

    @basic_attributes.setter
    def basic_attributes(self, arg: flyby_basic):
        self.basic_dataframe = arg

    @guess_attributes.setter
    def guess_attributes(self, arg: flyby_guess):
        self.guess_dataframe = arg

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
        columns = ['IN_OUT', 'v', 'v_inf', 'v_planet', 'ecc', 'sma', 'v_p', 'r_p', 'dv']
        df_in  = pd.DataFrame(columns)
        df_out = pd.DataFrame(columns)
        
        self._basic_dataframe = df()

    @guess_dataframe.setter
    def guess_dataframe(self, arg: flyby_guess):
        df = pd.DataFrame()
        self._guess_dataframe = df

    @refined_dataframe.setter
    def refined_dataframe(self, arg: flyby_refined):
        self._refined_dataframe = df

    def guess_powered_gravity_assist(self, v_i, v_f):
        if self.state is 'checked':
            pass
        else:
            self.check_gravity_assist(v_i, v_f)
        self._guess_attributes = guess_flyby(self._basic_powered_gravity_assist(v_i, v_f), self.planetary_node)

    def refined_powered_gravity_assist(self):
        self.state = 'refined'
        self.check_gravity_assist()
        self._vp_i, self._vp_f, self._raan, self._aop, self._inc, self._t_i, self._t_f, self._ecc_i, self._ecc_f, \
        self.planetary_node.soi_entry_position_body_ecliptic, self.planetary_node.soi_exit_position_body_ecliptic \
            = \
            pga_scalar_2_vector(self._v_inf_i, self._v_inf_f, self._vp_i, self._vp_f, self._sma_i, self._sma_f,
                                self._ecc_i, self._ecc_f, self._rp, self._body,
                                self.planetary_node.soi_periapsis_magnitude, self._epoch_periapsis)
