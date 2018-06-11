
from trajectory_tool.ale_solar_system.solar_system import SolarSystem
import datetime
import time
import numpy as np


ss = SolarSystem()

epoch = 2471151.500000


# use this http://aa.usno.navy.mil/data/docs/JulianDate.php
this = ss.coordinates_of('PLUTO', epoch)
distance = np.linalg.norm(this)
print(distance, 'km')