from .planetary_node import PlanetaryNode
from trajectory_tool.mnag_mission_analysis.departure_helper import *


class PlanetaryDeparture(object):
    def __init__(self, planetary_node:PlanetaryNode):
        self.planetary_node = planetary_node
        self._body = planetary_node.body
        self._epoch_periapsis = planetary_node.epoch_periapsis

    def calculate_refined(self, v_exit, r_p):
        return



