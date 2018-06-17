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
        self._epoch_entry = None
        self._soi_entry_magnitude = None
        self._epoch_exit = None
        self._soi_exit_magnitude = None

        # Refined attributes
        self._soi_entry_position = None
        self._soi_exit_position = None
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
    def epoch_periapsis(self):
        return self._epoch_periapsis

    @property
    def epoch_entry(self):
        if self.epoch_entry is None:
            return self._epoch_periapsis
        else:
            return self._epoch_entry

    @property
    def epoch_exit(self):
        if self.epoch_exit is None:
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
    def orbit_object_exit_epoch(self):
        if self._orbit_object_exit_epoch is None:
            return self._orbit_object_periapsis_epoch
        else:
            return self._orbit_object_exit_epoch

    @body.setter
    def body(self, arg: Body):
        self._body = arg
        self._periapsis_minimum = body_d_domain[body_string_lower(self.body)]['lower'] * u.km

    @epoch_periapsis.setter
    def epoch_periapsis(self, arg: datetime):
        self._epoch_periapsis = arg

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
        self._soi_entry_magnitude = soi(self._body, arg)

