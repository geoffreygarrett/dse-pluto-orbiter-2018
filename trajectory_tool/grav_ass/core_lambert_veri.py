#======================================================================================================================================
# TEST CASE
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

_test = core.TrajectoryTool()
__raw_itinerary1 = ['earth', 'jupiter', 'pluto']
__raw_itinerary2 = {'id': 5695,
                    'launch_date': datetime.datetime(2028, 12, 20, 0, 0),
                    'durations': [2.104, 17.29]
                    }
#%%
processed = _test.process_itinerary(__raw_itinerary2, __raw_itinerary1, _mode='plot')

# #initial parameters
# pars_init = processed[0]
#
# vs_init = pars_init['v']
# v_d_init = vs_init['d']
# v_d_init = np.linalg.norm(v_d_init)
#
# #lam_init = pars_init['l']
#
# epoch_init = pars_init['d']
#
# #Final parameters
# pars_fin = processed[1]
#
# vs_fin = pars_fin['v']
# v_a_fin = vs_fin['a']
# v_a_fin = np.linalg.norm(v_a_fin)
#
# epoch_fin = pars_fin['d']
#
# mu=mu_j
# d = 72008.182  #d in km

#%%
lam_init = (processed[1])['l']
lam_fin = (processed[2])['l']
sst = lam_init.sst
v_earth = ((processed[0])['v'])['p']
v_init = lam_init.v0
v_fin = lam_init.v1
v_fin_out = lam_fin.v0
posvec_init, velvec_init = state_to_vector(sst)
epoch_init = lam_init.epoch0
epoch_fin= lam_init.epoch1

sst2 = lam_fin.sst
posvec_fin, velvec_fin = state_to_vector(sst2)




#%%
import numpy as np
print('launch epoch (UTC) = ' + str(epoch_init))
print('launch epoch (MJD) = ' + str(epoch_init.mjd))
print('Earth Velocity = ' + str(v_earth))
print('departure velocity (heliocentric) = ' + str(v_init))
print('departure coordinates (heliocentric) = ' + str(posvec_init))
print('-                                                    -')
print('arrival epoch (UTC) = ' + str(epoch_fin))
print('arrival epoch (MJD) = ' + str(epoch_fin.mjd))
print('arrival velocity (heliocentric) = ' + str(v_fin))
print('arrival coordinates (heliocentric) = ' + str(posvec_fin))
print('departure velocity (heliocentric) = ' + str(v_fin_out) )


# # Gravity assist portion
#
# G = 6.67408*10**-11*10**-9 #km^3/kg/s^2
# Mj = 1.898E27                 #kg
# mu_j = G*Mj                   #km^3...
#
# ass = Hyperbolic_calculator_resources_v0p3.grav_ass
#
# v_venus = ((processed[1])['v'])['p']
# d=72008.
# pars = ass(v_fin, v_venus, v_fin_out, d, mu_j, 'full')
# #%%
# newlam = (processed[3])['l']
# sst_jup = newlam.sst
# posvec_jup, another = state_to_vector(sst_jup)
#
#
# v_rel_in = pars[-2]
# v_rel_out = pars[1]

#===========================================================================================================================================
# Gravity assist portion

G = 6.67408*10**-11*10**-9 #km^3/kg/s^2
Mj = 1.898E27                 #kg
mu_j = G*Mj                   #km^3...

ass = Hyperbolic_calculator_resources_v0p3.grav_ass

v_venus = ((processed[1])['v'])['p']



v_in = v_fin
v_out_req = v_fin_out
d=72008
mu_j = G*Mj                   #km^3...

pars = ass(v_fin, v_venus, v_fin_out, d, mu_j, 'full')

v_out = pars[0]
