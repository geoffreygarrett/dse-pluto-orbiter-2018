#======================================================================================================================================
# Creates variables to input to the GMAT script
#======================================================================================================================================
import numpy as np
import trajectory_tool.core
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
from trajectory_tool.core import TrajectoryTool
from trajectory_tool.genetic_algorithim_analysis.genetic_algorithim import Chromosome
import os


import datetime
from collections import namedtuple
#from grav_ass import Hyperbolic_calculator_resources_v0p3

def save_first_cases(chromosome, writedata, data_titles, filename, bigfolder):
    base_path = os.path.dirname(os.path.realpath(__file__))
    folder = "chromosome_{}".format(chromosome)
    filepath = os.path.join(base_path,bigfolder, folder, filename)

    # CREATES DIRECTORY IF IT DOESN NOT EXIST
    directory = os.path.join(base_path,bigfolder, folder)
    if not os.path.exists(directory):
        os.makedirs(directory)
    print("Saving {}".format(filepath))
    print(filepath)
    with open(filepath, 'w') as f:
        #f.write('Chromosome: ' + chromosome)
        f.write(data_titles)
        f.write(writedata)

def save_other_cases(chromosome, writedata, data_titles, filename, bigfolder):
    base_path = os.path.dirname(os.path.realpath(__file__))
    folder = "chromosome_{}".format(chromosome)
    filepath = os.path.join(base_path,bigfolder, folder, filename)

    print("Saving {}".format(filepath))
    print(filepath)
    with open(filepath, 'a') as f:
        f.write(writedata)


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
chromosome = '1057 50756 00000 7998'
__raw_itinerary2, __raw_itinerary1 = Chromosome.mapper(chromosome)
# TEST EJP -------------------------------------------------------------------------------------------------
# __raw_itinerary1 = ['earth', 'jupiter', 'pluto']
# __raw_itinerary2 = {'id': 0,
#                     'launch_date': datetime.datetime(2028, 12, 20, 0, 0),
#                     'durations': [1.8136986301369864, 21.915068493150685]
#                     }
# ----------------------------------------------------------------------------------------------------------

processed = _test.process_itinerary(__raw_itinerary2, __raw_itinerary1, _mode='plot2D', _grav_ass=True)
tt = TrajectoryTool()
#%%
#opass = TrajectoryTool.optimise_gravity_assist
#opass(TrajectoryTool, v_s_i, v_s_f, v_p, body, epoch)
v_s_i = processed[1]['v']['a']
v_s_f_initial = processed[1]['v']['d']
v_planet_initial = processed[1]['v']['p']
body_assisting = processed[1]['b']
body_next = processed[1+1]['b']
epoch_entry = processed[1]['d']
epoch_next_body = processed[1+1]['d']

NotEarthData=[]
for i in range(len(processed)):
    bodydata = []
    print('-                                             -')
    if i == 0:
        EarthData = []
        ss0 = processed[1]['l'].ss0
        ss1 = processed[1]['l'].ss1
        posvec, velvec = state_to_vector(ss0)
        body = str(processed[0]['b'])
        #bodydata['Body'] = str(body)
        print('Body = ', body)
        EarthData.append(body[:-4])

        launch_epoch_str = str(ss0.epoch)
        #bodydata['Launch_Epoch'] = launch_epoch_str
        print('Epoch at Launch = ' + launch_epoch_str)
        EarthData.append(launch_epoch_str)
        print('Position Vector At Launch = '+str(posvec))
        EarthData.append(str(posvec))
        print('Velocity Vector At Launch = ' + str(processed[i]['v']['d']))
        EarthData.append(str(processed[i]['v']['d']))
        #lambert_pseudo_data['Body_launch'] = bodydata
        print('Planetary velocity vector = '+str(processed[i]['v']['p']))
        EarthData.append(str(processed[i]['v']['p']))

    if i != 0 and i!= len(processed)-1:
        bodydata = []
        print('Body = ' + str(processed[i]['b']))
        bodydata.append((str(processed[i]['b']))[:-4])
        print('Epoch of arrival = ' + str(processed[i]['d']))
        bodydata.append(str(processed[i]['d']))
        ss = processed[i]['l'].ss1
        print('Position vector = ', state_to_vector(ss)[0])
        bodydata.append(state_to_vector(ss)[0])
        print('Incoming velocity vector = '+str(processed[i]['v']['a']))
        bodydata.append(str(processed[i]['v']['a']))
        print('Outgoing velocity vecotr = ', processed[i]['v']['d'])
        bodydata.append(processed[i]['v']['d'])
        print('Planetary velocity vector = '+str(processed[i]['v']['p']))
        bodydata.append(str(processed[i]['v']['p']))
        NotEarthData.append(bodydata)
    if i == len(processed)-1:
        print('Body = ' + str(processed[i]['b']))
        bodydata.append((str(processed[i]['b']))[:-4])
        print('Epoch of arrival = ' + str(processed[i]['d']))
        bodydata.append(str(processed[i]['d']))
        ss = processed[i]['l'].ss1
        print('Position vector = ', state_to_vector(ss)[0])
        bodydata.append(state_to_vector(ss)[0])
        print('Incoming velocity vector = '+str(processed[i]['v']['a']))
        bodydata.append(str(processed[i]['v']['a']))
        #print('Epoch of exit = ', processed[i]['v']['d'])
        print('Planetary velocity vector = '+str(processed[i]['v']['p']))
        bodydata.append(str(processed[i]['v']['p']))
        NotEarthData.append(bodydata)

refined_ass = tt.refined_gravity_assist(v_s_i, v_s_f_initial, v_planet_initial, body_assisting, body_next, epoch_entry, epoch_next_body,mode='plot3D', verification=True)
#less_refined_ass = tt.optimise_gravity_assist(v_s_i, v_s_f_initial, v_planet_initial, body_assisting, epoch_entry, mode='plot3D', verification=True)

#SAVE STUFF:
#define save folder

#%%
bigfolder = 'verification_data_dump'
#Save Lambert Stuff
lambertdata=''
lambertdata += "{:<50} \t {:<50} \t {:<50}\n".format(str(v_s_i),
                              str(v_s_f_initial),
                               'bellend')

data_titles = "{:<50} \t {:<50} \t {:<50} \t {:<50} \t {:<50} \t {:<50}\n".format('Body',
                                                  'Arrival Epoch',
                                                  'Position Vector (km)',
                                                    'Incoming Velocity Vector (km/s)',
                                                    'Outgoing Velocity Vector (km/s)',
                                                    'Planetary Velocity Vector (km/s)')

#SaveEarth Cases
EarthState =  "{:<50} \t {:<50} \t {:<50} \t {:<50} \t {:<50} \t {:<50}\n".format(str(EarthData[0]),
                              str(EarthData[1]),
                              str(EarthData[2]),
                                                    '--',
                                                    str(EarthData[3]),
                                                    str(EarthData[4]))

save_first_cases(chromosome, EarthState, data_titles, 'somedata.txt', bigfolder)

#Save Jupiter Cases
JupiterState =  "{:<50} \t {:<50} \t {:<50} \t {:<50} \t {:<50} \t {:<50}\n".format(str(NotEarthData[0][0]),
                              str(NotEarthData[0][1]),
                              str(NotEarthData[0][2]),
                                                    NotEarthData[0][3],
                                                    str(NotEarthData[0][4]),
                                                    str(NotEarthData[0][5]))

save_other_cases(chromosome, JupiterState, data_titles, 'somedata.txt', bigfolder)

#Save Pluto Cases
PlutoState =  "{:<50} \t {:<50} \t {:<50} \t {:<50} \t {:<50} \t {:<50}\n".format(str(NotEarthData[1][0]),
                              str(NotEarthData[1][1]),
                              str(NotEarthData[1][2]),
                                                    NotEarthData[1][3],
                                                    '--',
                                                    str(NotEarthData[1][4]))

save_other_cases(chromosome, PlutoState, data_titles, 'somedata.txt', bigfolder)

#%%
a=6992863713.300394
b=1336628199.534705
c=-1449298454.956082

#print(str(np.sqrt(a**2+b**2+c**2)))
