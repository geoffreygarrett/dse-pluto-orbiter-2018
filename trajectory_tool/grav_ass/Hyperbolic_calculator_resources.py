# -*- coding: utf-8 -*-
"""
Created on Thu May 24 16:21:40 2018

@author: matth
"""

import numpy as np

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

#Vp = 12740./1000.00       #km/s
#phi_p = 2.4       #deg
#Vsi = 9.470       #km/s
#phi_si = 39.2     #deg
#d = -2500000.     #km
#mu = 1.26686*10**17*10**-9 #km^3/sihdfsi

#function takes some inputs, then outputs the final spacecraft velocity vector and flight path angle
def gravass(Vp, phi_p, Vsi, phi_si, d, mu_p):
      
      if d > 0:
            phi_si = -phi_si
      
      #Define planet and sc velocity vectors
      Vp_vec = vel2vec(Vp,phi_p)
      Vsi_vec = vel2vec(Vsi, phi_si)
      
      #Define initial velocity vector of planet wrt sc
      Vspi_x = Vsi_vec[0] - Vp_vec[0]
      Vspi_y = Vsi_vec[1] - Vp_vec[1]
      Vspi_vec = np.array([Vspi_x, Vspi_y])
      
      #Define magnitude of Vspi_vec
      Vspi = vec2mag(Vspi_vec)
      #Define excess velocity (roughly equal to Vspi)
      Vinf = Vspi
      
      #Define theta_i (incoming asymptote angle)
      theta_i = np.arctan2(Vspi_y,Vspi_x)
      theta_i_deg = np.degrees(theta_i)
      
      #define b (impact parameter), a (SMA), e (eecentricity) and delta (turning angle)
      b = d*np.sin(theta_i)
      a = -mu_p/Vinf**2
      e = np.sqrt(1+b**2/a**2)
      delta = 2*np.arcsin(1/e)
      
      #define theta_f (outgoing asymptote angle)
      theta_f = theta_i-delta
      
      #Define final velocity vector of planet wrt sc
      Vspf = vel2vec(Vspi, np.degrees(theta_f))
      
      #Define final velocity vector of spacecraft (and its magnitude)
      Vsf_x = Vspf[0] + Vp_vec[0]
      Vsf_y = Vspf[1] + Vp_vec[1]
      Vsf_vec = np.array([Vsf_x,Vsf_y])
      
      Vsf = vec2mag(Vsf_vec)
      
      #Define Final flight path angle of sc
      phi_sf = np.arctan2(Vsf_y,Vsf_x)
      
      return Vsf_vec, Vsf, np.degrees(phi_sf)



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















