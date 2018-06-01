# -*- coding: utf-8 -*-
"""
Created on Thu May 24 16:21:40 2018

@author: matth
"""

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




#Given knowledge of initial: 
      # s/c heliocentric velocity               Vsi
      # s/c heliocentric flight path angle      phi_si
      # planet velocity                         Vp
      # planet flight path angle                phi_p
      # miss distance along orbit (+ve ahead)   d
      # gravitational parameter of the planet   mu
      
#then calculate final:
      # s/c heliocentric velocity               Vsf
      # s/c heliocentric flight path angle      phi_sf


#======================================================================================================================================
# Define some simple functions
#======================================================================================================================================

#function to get vector magnitude based on x and y components
def get_magnitude(x,y): return np.sqrt(x**2 + y**2)

#Function to convert velocity and flight path angle to vector. input angle in degrees
def vel2vec(V, phi):
      vec = np.array([0,0], dtype=np.float)
      vec[0] = V*np.cos(np.radians(phi))
      vec[1] = V*np.sin(np.radians(phi))
      return vec

#Function to get vector magnitude
def vec2mag(vector):
      return np.sqrt(vector[0]**2 + vector[1]**2)

#Function for calculating delta V for first burn of hohmann
def get_hoh_dv1(mu,r1,r2):
      return np.sqrt(mu/r1)*(np.sqrt(2*r2/(r1+r2))-1)

#functions to convert from C3 to DV
def C32dv(C3, mu, r): return np.sqrt(C3 + 2*mu/r) - np.sqrt(mu/r)

#function to get circular velocity given mu and r
def Vcirc(mu,r): return np.sqrt(mu/r)

#Calculate (and print) all keplerian elements from a, e and r
def dokep(a, e, r, mu):
      Ap = a*(1+e)
      Pe = a*(1-e)
      T = 2*np.pi*np.sqrt(a**3/mu)
      print('SMA = ' + str(a/AU) + ' AU')
      print('ECC = ' + str(e))
      print('RAD = ' + str(r/AU) + ' AU')
      print('APO = ' + str(Ap/AU) + ' AU')
      print('PER = ' + str(Pe/AU) + ' AU')
      print('PERIOD = ' + str(T/60/60/24/365.25) + ' years')
      
      return a, e, r, Ap, Pe, T

#Calculate orbital velocity given a r and mu
def get_orb_V(a,r,mu): return np.sqrt(2*mu/r - mu/a)


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

#======================================================================================================================================
# Define some simple constants
#======================================================================================================================================
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
#======================================================================================================================================
# Function to do the gravity assist
#======================================================================================================================================

#TEST CASE PARAMETERS
# _test = core.TrajectoryTool()
# __raw_itinerary1 = ['earth', 'venus', 'jupiter', 'pluto']
# __raw_itinerary2 = {'id': 4332,
#                     'launch_date': datetime.datetime(2024, 1, 1, 0, 0),
#                     'durations': [0.2102, 1.0114, 11.7784]
#                     }
#
# processed = _test.process_itinerary(__raw_itinerary2, __raw_itinerary1, _mode='fast')
#
# jupars = processed['2']
# vs_ju = jupars['v']
# v_in = vs_ju['a']
# v_out_req = vs_ju['d']
# v_bod = vs_ju['p']

#
# mu=mu_j
# d = 72008.182  #d in km


# Function takes incoming and outgoing velocity vector defined by lambert==============
# and with d and mu computes the actual v_out==========================================
# if d +ve, increase dv.===============================================================

def grav_ass(v_in, v_bod, v_out_req, d, mu):

    #define required turning angle
    #delta_req = angle_check(v_in, v_out_req)

    #define incoming theta
    theta_i = angle_check(v_in, v_bod)

    #Calculate achieved delta
    v_rel_in = v_in - v_bod                  # relative velocity
    v_inf = np.linalg.norm(v_rel_in)         # define v_inf in same units as v_rel_in

    a = -mu/v_inf**2                                      #claculate SMA
    b = d*np.sin(theta_i)                                    #calculate b vector magnitude
    e = np.sqrt(1 + b**2/a**2)                            #calculate eccentricity

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

    return v_out, delta


#======================================================================================================================================
# Function to do escape burn
#======================================================================================================================================

#r1 = a_e                #helio radius 1
#r2 = a_j                #helio radius 2
#mu1=mu_s                # solar mu
#mu2 = mu_e              # planet mu
#ri = 6378 + 250               # sc orbital radius around planet

def inter_inject(r1,r2,mu1,mu2,ri):       #r1 = orbital radius of start planet, r2 = orbital radius of target planet, mu1 = mu of sun, mu2 = mu of exit planet, ri = initial orbital radius

      DV1 = get_hoh_dv1(mu1,r1,r2)
      C3 = DV1**2
      DVb = C32dv(C3,mu2,ri)
      
      return DVb

#print(inter_inject(r1,r2,mu1,mu2,ri))
#======================================================================================================================================
# Function to calculate velocity magnitude of a body and it's flight path angle given keplerian elements
#======================================================================================================================================

#mu = mu_s
#r_pe = a_e
#r_ap = a_j
#r = 1*AU

def kep2V_phi(mu,r_pe,r_ap,r):       # mu = mu of central body, r = radius at the time
      #calc SMA from Pe and Ap
      a = (r_pe + r_ap)/2
      #Calc e from P and SMA
      e = 1 - r_pe/a
      #Calc V at radius r from mu and a
      V = np.sqrt(2*mu/r - mu/a)
      #Calc h (specific relative angular momentum) from a e and mmu
      h = np.sqrt(a*(1-e**2)*mu)
      #Calc phi from h r and V
      phi = np.degrees(np.arccos(h/(r*V)))
      
      return V, phi

#======================================================================================================================================
# Function to convert a velocity and flight path angle into keplerian elements
#======================================================================================================================================

#mu = mu_s
#phi = np.radians(42.626)
#V = 17.479
#r = 3*AU

#mu = central body mu (ie sun)
def V_phi2kep(V, phi, r, mu):
      a = (2/r - V**2/mu)**-1
      e = np.sqrt(1-(np.cos(phi)*r*V)**2/(a*mu))
      return a, e















