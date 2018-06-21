from poliastro.bodies import Body
from datetime import datetime
from .tools.orbital_mechanics import *
from .tools.misc import *
from .config import *
import astropy.units as u
from astropy import time
from poliastro.twobody import Orbit


class PlanetaryNode(object):
    """
    Defines a planetary node that is used to defined a interplanetary manoeuvre.
    """
    def __init__(self):
        self._order_of_states = ['basic', 'rough', 'refined']
        self.__state = []

        # Basic attributes
        self._body = None
        self._epoch_periapsis = None
        self._periapsis_minimum = None
        self._soi_periapsis_magnitude = None
        self._orbit_object_periapsis_epoch = None

        # Rough attributes
        self._v_entry = None
        self._epoch_entry = None
        self._soi_entry_magnitude = None
        self._soi_entry_position_heliocentric = None
        self._v_exit = None
        self._epoch_exit = None
        self._soi_exit_magnitude = None
        self._soi_exit_position_heliocentric = None

        # Refined attributes
        self._soi_entry_position_body_ecliptic = None
        self._soi_exit_position_body_ecliptic = None
        self._orbit_object_entry_epoch = None
        self._orbit_object_exit_epoch = None

        # States
        self._state_basic = False
        self._state_rough = False
        self._state_refined = False

    def __repr__(self):
        return str(self._body) + ' PlanetaryNode'

    @property
    def body(self):
        return self._body

    @property
    def periapsis_minimum(self):
        return self._periapsis_minimum

    @property
    def v_entry(self):
        return self._v_entry

    @property
    def v_exit(self):
        return self._v_exit

    @property
    def epoch_periapsis(self):
        return self._epoch_periapsis

    @property
    def epoch_entry(self):
        if self._epoch_entry is None:
            return self._epoch_periapsis
        else:
            return self._epoch_entry

    @property
    def epoch_exit(self):
        if self._epoch_exit is None:
            return self._epoch_periapsis
        else:
            return self._epoch_exit

    @property
    def soi_periapsis_magnitude(self):
        return self._soi_periapsis_magnitude

    @property
    def soi_entry_magnitude(self):
        if self.soi_entry_magnitude is None:
            return self._soi_periapsis_magnitude
        else:
            return self._soi_entry_magnitude

    @property
    def soi_exit_magnitude(self):
        if self.soi_exit_magnitude is None:
            return self._soi_periapsis_magnitude
        else:
            return self._soi_exit_magnitude

    @property
    def orbit_object_periapsis_epoch(self):
        return self._orbit_object_periapsis_epoch

    @property
    def orbit_object_entry_epoch(self):
        if self._orbit_object_entry_epoch is None:
            return self._orbit_object_periapsis_epoch
        else:
            return self._orbit_object_entry_epoch

    @property
    def v_planet_i(self):
        if self._orbit_object_entry_epoch:
            return self._orbit_object_entry_epoch.state.v
        else:
            return self._orbit_object_periapsis_epoch.state.v

    @property
    def v_planet_f(self):
        if self._orbit_object_entry_epoch:
            return self._orbit_object_exit_epoch.state.v
        else:
            return self._orbit_object_periapsis_epoch.state.v

    @property
    def orbit_object_exit_epoch(self):
        if self._orbit_object_exit_epoch is None:
            return self._orbit_object_periapsis_epoch
        else:
            return self._orbit_object_exit_epoch

    @property
    def soi_exit_position_heliocentric(self):
        return self.orbit_object_exit_epoch.r + self.soi_exit_position_body_ecliptic

    @property
    def soi_entry_position_heliocentric(self):
        return self.orbit_object_entry_epoch.r + self.soi_entry_position_body_ecliptic

    @property
    def soi_exit_position_body_ecliptic(self):
        if self._soi_exit_position_body_ecliptic is None:
            return np.array([0, 0, 0]) * u.km
        return self._soi_exit_position_body_ecliptic

    @property
    def soi_entry_position_body_ecliptic(self):
        if self._soi_entry_position_body_ecliptic is None:
            return np.array([0, 0, 0]) * u.km
        return self._soi_entry_position_body_ecliptic

    @body.setter
    def body(self, arg: Body):
        self._body = arg
        self._periapsis_minimum = body_d_domain[body_string_lower(self.body)]['lower'] * u.km

    @soi_entry_position_body_ecliptic.setter
    def soi_entry_position_body_ecliptic(self, arg:np.ndarray):
        self._soi_entry_position_body_ecliptic = arg

    @soi_exit_position_body_ecliptic.setter
    def soi_exit_position_body_ecliptic(self, arg:np.ndarray):
        self._soi_exit_position_body_ecliptic = arg

    @epoch_periapsis.setter
    def epoch_periapsis(self, arg: datetime):
        self.soi_periapsis_magnitude = arg
        self._epoch_periapsis = arg
        self._orbit_object_periapsis_epoch = Orbit.from_body_ephem(self.body, time.Time(arg, scale='tdb'))

    @epoch_entry.setter
    def epoch_entry(self, arg: datetime):
        self.soi_entry_magnitude = arg
        self._epoch_entry = arg
        self._orbit_object_entry_epoch = Orbit.from_body_ephem(self.body, time.Time(arg, scale='tdb'))

    @epoch_exit.setter
    def epoch_exit(self, arg: datetime):
        self.soi_exit_magnitude = arg
        self._epoch_exit = arg
        self._orbit_object_exit_epoch = Orbit.from_body_ephem(self.body, time.Time(arg, scale='tdb'))

    @soi_entry_magnitude.setter
    def soi_entry_magnitude(self, arg: datetime):
        self._soi_entry_magnitude = soi(self._body, arg)

    @soi_exit_magnitude.setter
    def soi_exit_magnitude(self, arg: datetime):
        self._soi_exit_magnitude = soi(self._body, arg)

    @soi_periapsis_magnitude.setter
    def soi_periapsis_magnitude(self, arg: datetime):
        self._soi_periapsis_magnitude = soi(self._body, arg)

    @v_entry.setter
    def v_entry(self, arg):
        self._v_entry = arg

    @v_exit.setter
    def v_exit(self, arg):
        self._v_exit = arg

