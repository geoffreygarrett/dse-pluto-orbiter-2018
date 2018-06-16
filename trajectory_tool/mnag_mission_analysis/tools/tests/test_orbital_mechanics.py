from poliastro.bodies import Earth
from trajectory_tool.mnag_mission_analysis.tools import orbital_mechanics
from datetime import datetime


def test_soi():
    """
    Curtis, H. (2010). Orbital Mechanics for Engineering Students (2nd ed.) (2nd ed., p. 358). Elsevier Science Limited.
    :return:
    """
    test = 0.924e6
    body_earth = Earth
    calc = orbital_mechanics.soi(body_earth, datetime.today()).value
    assert abs(calc - test) <= 0.05 * test
