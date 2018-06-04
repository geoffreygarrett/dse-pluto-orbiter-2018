# -*- coding: utf-8 -*-
"""
CHANGELOG

- Removed uneccessary function definitions and imports
- Changed grav_ass to have a new input (mode) allowing short results or a variable dump
- Cleaned up and commented out the test case

@author: matth
"""

# Function to calculate angle between 2 vectors
def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            angle_between((1, 0, 0), (1, 0, 0))
            0.0
            angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    angle = np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
    return angle, (2*np.pi-angle)

# find vector rotation matrix
def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2.0)
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])


# Function to get the turning angle/ theta
def angle_check(vec1, vec2):

      # calculate all possible angles and list them

      start_angles = angle_between(vec1, vec2)
      angles = [start_angles[0], start_angles[1], -start_angles[0], -start_angles[1]]

      # crossproduct of vectors gives axis of rotation
      crossprod = np.cross(vec1, vec2)

      # loop to list feasible angles
      feasible_angles = []
      for angle in angles:

            # create rotation matrix
            rot_matrix = rotation_matrix(crossprod, angle)

            # perform rotation
            vec1_rotated = np.dot(rot_matrix, vec1)
            vec1_rotated_unit = vec1_rotated / np.linalg.norm(vec1_rotated)
            vec2_unit = vec2 / np.linalg.norm(vec2)

            # check if parallel, if true then add to feasible angles
            parallel_check = np.dot(vec1_rotated_unit, vec2_unit)

            if 0.99 < parallel_check < 1.01:
                  feasible_angles.append(abs(angle))

      # return the smallest of the feasible angles
      return min(feasible_angles)

#%%
#======================================================================================================================================
# Function to do the gravity assist
#======================================================================================================================================

#For mode, enter 'short' for just v_out and delta. Enter 'full' for a dump of all pertinent parameters
def grav_ass(v_in, v_bod, v_out_req, d, mu, mode):

    #define required turning angle
    delta_req = angle_check(v_in, v_out_req)

    #define incoming theta
    theta_i = angle_check(v_in, v_bod)

    #Calculate achieved delta
    v_rel_in = v_in - v_bod                  # relative velocity
    v_inf = np.linalg.norm(v_rel_in)         # define v_inf in same units as v_rel_in

    a = -mu/v_inf**2                                      #claculate SMA
    b = d*np.sin(theta_i)                                    #calculate b vector magnitude
    e = np.sqrt(1 + b**2/a**2)                            #calculate eccentricity
    Rpe = a*(1-e)                                           #calculate peri radius
    #Define turningn angle delta
    delta = 2*np.arcsin(1/e)

    if d<0:
        delta = -delta
        theta_i = -theta_i

    #Define outgoing theta
    theta_f = theta_i - delta

    #Define outgoing relative velocity
    v_rel_out_req = v_out_req - v_bod
    rot_axis = np.cross(v_rel_in, v_rel_out_req)
    rot_matrix = rotation_matrix(rot_axis, theta_f)
    v_rel_out = np.dot(rot_matrix, v_rel_in)

    #Find heliocentric velocity out
    v_out = v_rel_out + v_bod

    if mode == 'short': return v_out, delta

    if mode == 'full': return v_out, v_rel_out, rot_matrix, rot_axis, v_rel_out_req, theta_f, theta_i, delta, e, b, a, Rpe, v_inf, v_rel_in, delta_req


# #======================================================================================================================================
# # TEST CASE
# #======================================================================================================================================
# import numpy as np
# import core
# import matplotlib.pyplot as plt
# from matplotlib import cm
#
# from poliastro.plotting import OrbitPlotter
# from poliastro.bodies import Earth, Mars, Sun, Pluto, Venus, Jupiter
# from poliastro.twobody import Orbit
# from poliastro import iod
# from astropy.coordinates import solar_system_ephemeris
# from astropy import time
# import math
# from astropy import units as u
#
# import datetime
# from collections import namedtuple
#

# G = 6.67408*10**-11*10**-9 #km^3/kg/s^2
# AU = 1.496E8                  #km
#
# #Earth
# Me = 5.972E24                 #kg
# mu_e = G*Me                   #km^3...
# a_e = AU
# R_e = 6378.                   #km
#
# #Sun
# Ms = 1.989E30                  #kg
# mu_s = G*Ms                   #km^3...
#
# #Jupiter
# Mj = 1.898E27                 #kg
# mu_j = G*Mj                   #km^3...
# #a_j = 5.05*AU                 #this is fake news, the real one is below
# a_j = 5.2044*AU
# R_j = 317.7*R_e               #km
#
# #Venus
# Mv = 4.867E24                 #kg
# mu_v =  G*Mv                  #km^3...
# a_v = 0.723332*AU
#
# # Pluto
# Mp = 1.309E22                 #kg
# mu_p = G*Mp                   #km^3...
# a_p = 39.48*AU                #km
# e = 0.2488
#
# _test = core.TrajectoryTool()
# __raw_itinerary1 = ['earth', 'venus', 'jupiter', 'pluto']
# __raw_itinerary2 = {'id': 4332,
#                     'launch_date': datetime.datetime(2024, 1, 1, 0, 0),
#                     'durations': [0.2102, 1.0114, 11.7784]
#                     }
#
# processed = _test.process_itinerary(__raw_itinerary2, __raw_itinerary1, _mode='fast')
#
# jupars = processed[2]
# vs_ju = jupars['v']
# v_in = vs_ju['a']
# v_out_req = vs_ju['d']
# v_bod = vs_ju['p']
#
#
# mu=mu_j
# d = 72008.182  #d in km
#
# #Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto in km
# min_rads = [2539.7, 6427.67, 6540.7568, 3552.2672, 72008.182, 61405.511, 26088.564, 25152.092, 1287]
# mus = [3.3E23*G, mu_v, mu_e, 6.39E23*G, mu_j, 5.68E26*G, 8.68E25*G, 1.02E26*G, mu_p]
# #%%
# for rad in min_rads:
#     fullpars = grav_ass(v_in, v_bod, v_out_req, d, mu, 'full')
#
#     a = fullpars[-4]
#     e = fullpars[-6]
#
#     Rpe = a*(1-e)
