from .planetary_node import PlanetaryNode


class PlanetaryDeparture(object):
    def __init__(self, planetary_node:PlanetaryNode):
        self.planetary_node = planetary_node
        self._body = planetary_node.body
        self._epoch_periapsis = planetary_node.epoch_periapsis
