#======================================================================================================================================
# Creates variables to input to the GMAT script
#======================================================================================================================================
import numpy as np
import core
import matplotlib.pyplot as plt
from matplotlib import cm

from poliastro.plotting import OrbitPlotter
from poliastro.bodies import Earth, Mars, Sun, Pluto, Venus, Jupiter
from poliastro.twobody import Orbit
from poliastro import iod
from astropy.coordinates import solar_system_ephemeris
from astropy import time
import math
from astropy import units as u
from core import TrajectoryTool


import datetime
from collections import namedtuple
from grav_ass import Hyperbolic_calculator_resources_v0p3


def state_to_vector(ss):
    r, v = ss.rv()
    x, y, z = r.to(u.km).value
    vx, vy, vz = v.to(u.km / u.s).value
    return np.array([x, y, z]), np.array([ vx, vy, vz])

G = 6.67408*10**-11*10**-9 #km^3/kg/s^2
AU = 1.496E8                  #km

#Earth
Me = 5.972E24                 #kg
mu_e = G*Me                   #km^3...
a_e = AU
R_e = 6378.                   #km

#Sun
Ms = 1.989E30                  #kg
mu_s = G*Ms                   #km^3...

#Jupiter
Mj = 1.898E27                 #kg
mu_j = G*Mj                   #km^3...
#a_j = 5.05*AU                 #this is fake news, the real one is below
a_j = 5.2044*AU
R_j = 317.7*R_e               #km

#Venus
Mv = 4.867E24                 #kg
mu_v =  G*Mv                  #km^3...
a_v = 0.723332*AU

# Pluto
Mp = 1.309E22                 #kg
mu_p = G*Mp                   #km^3...
a_p = 39.48*AU                #km
e = 0.2488



#%%
_test = TrajectoryTool()

# TEST EJP -------------------------------------------------------------------------------------------------
__raw_itinerary1 = ['earth', 'jupiter', 'pluto']
__raw_itinerary2 = {'id': 0,
                    'launch_date': datetime.datetime(2028, 12, 20, 0, 0),
                    'durations': [1.67397, 22.96712]
                    }
# ----------------------------------------------------------------------------------------------------------

processed = _test.process_itinerary(__raw_itinerary2, __raw_itinerary1, _mode='plot', _grav_ass=True)

jup_part = processed[1]
epoch = (jup_part['l']).epoch1
print('Epoch of arrival = ' + str(epoch))
#%%

v_s_i = jup_part['v']['a']
v_s_f = jup_part['v']['d']
v_p = v_s_f = jup_part['v']['p']
body = jup_part['b']
epoch = jup_part['d']
tt = TrajectoryTool()

#opass = TrajectoryTool.optimise_gravity_assist
#opass(TrajectoryTool, v_s_i, v_s_f, v_p, body, epoch)
tt.optimise_gravity_assist(v_s_i, v_s_f, v_p, body, epoch)
