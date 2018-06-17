from .planetary_node import PlanetaryNode


class PlanetaryRendezvous(object):
    def __init__(self, planetary_node: PlanetaryNode):
        self._body = planetary_node.body
        self._epoch_periapsis = planetary_node.epoch_periapsis
