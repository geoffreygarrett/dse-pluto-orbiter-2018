from poliastro.bodies import Sun, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto, Body
from datetime import datetime
from .tools.orbital_mechanics import *


class PlanetaryNode(object):
    """
    Defines a planetary node that is used to defined a interplanetary manoeuvre.
    """
    def __init__(self):
        self._body = None
        self._epoch_rp = None

    @property
    def body(self):
        return self._body

    @property
    def epoch(self):
        return self._epoch_rp

    @body.setter
    def body(self, arg: Body):
        self._body = arg

    @epoch.setter
    def epoch(self, arg: datetime):
        self._epoch_rp = arg

    @property
    def soi_rp(self):
        return soi(self._body, self._epoch_rp)
