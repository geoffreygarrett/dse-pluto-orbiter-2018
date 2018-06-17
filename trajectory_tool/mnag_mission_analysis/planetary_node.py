from poliastro.bodies import Body
from datetime import datetime
from .tools.orbital_mechanics import *
from .tools.misc import *
from .config import *
import astropy.units as u


class PlanetaryNode(object):
    """
    Defines a planetary node that is used to defined a interplanetary manoeuvre.
    """
    def __init__(self):
        self._order_of_states = ['basic', 'rough', 'refined']
        self.__state = None

        # Basic attributes
        self._body = None
        self._epoch_periapsis = None
        self._periapsis_minimum = None
        self._soi_periapsis_magnitude = None

        # Rough attributes
        self._epoch_entry = None
        self._soi_entry_magnitude = None
        self._epoch_exit = None
        self._soi_exit_magnitude = None

        # Refined attributes
        self._soi_entry_position = None
        self._soi_exit_position = None

        # States
        self._state_basic = False
        self._state_rough = False
        self._state_refined = False

    @property
    def state(self):
        return self._state

    @property
    def _state(self):
        return self.__state[-1]

    @property
    def _state_basic(self):
        return self._state_basic

    @property
    def _state_rough(self):
        return self._state_rough

    @property
    def _state_refined(self):
        return self._state_refined

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
        return self._epoch_entry

    @property
    def epoch_exit(self):
        return self._epoch_exit

    @property
    def soi_periapsis_magnitude(self):
        return self._soi_periapsis_magnitude

    @property
    def soi_entry_magnitude(self):
        return self._soi_entry_magnitude

    @property
    def soi_exit_magnitude(self):
        return self._soi_exit_magnitude

    @body.setter
    def body(self, arg: Body):
        self._periapsis_minimum = body_d_domain[body_string_lower(self.body)] * u.km
        if self.epoch_periapsis:
            self._state_basic = True
        self._body = arg

    @epoch_periapsis.setter
    def epoch_periapsis(self, arg: datetime):
        if self.body:
            self._state_basic = True
        self._epoch_periapsis = arg

    @epoch_entry.setter
    def epoch_entry(self, arg: datetime):
        self.soi_entry_magnitude = arg
        if self.epoch_exit:
            self._state_rough = True
        self._epoch_entry = arg

    @epoch_exit.setter
    def epoch_exit(self, arg: datetime):
        self.soi_exit_magnitude = arg
        if self.epoch_entry:
            self._state_rough = True
        self._epoch_exit = arg

    @soi_entry_magnitude.setter
    def soi_entry_magnitude(self, arg: datetime):
        self._soi_entry_magnitude = soi(self._body, arg)

    @soi_exit_magnitude.setter
    def soi_exit_magnitude(self, arg: datetime):
        self._soi_entry_magnitude = soi(self._body, arg)

    @_state_basic.setter
    def _state_basic(self, arg: True):
        self._state('basic')
        self._state_basic = arg

    @_state_rough.setter
    def _state_rough(self, arg: True):
        self._state('rough')
        self._state_rough = arg

    @_state_refined.setter
    def _state_refined(self, arg: True):
        self._state('refined')
        self._state_refined = arg

    @_state.setter
    def _state(self, arg):
        _state_error = "Maintain the order of states: \n {}".format(self._order_of_states)
        if self._state == arg:
            pass
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


