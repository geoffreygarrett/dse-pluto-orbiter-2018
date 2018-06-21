import numpy as np
import matplotlib.pyplot as plt

#Function to get periapsis velocity on incoming hyperbolic
def get_Vp_hyp(v_inf, mu2, rp): return np.sqrt(v_inf**2 + (2*mu2)/rp)

#Function to get periapsis velocity on required capture orbit
def get_Vp_cap(mu2, e_cap, rp): return np.sqrt((mu2*(1+e_cap))/rp)

#Velocity to get Dv insertion after a given an electric impulse
def getDvInsertion(Va_initial, V_pluto, DV_electric, mu2, rp, e_cap):
    V_inf_initial = Va_initial - V_pluto
    DV_electric_vec = DV_electric * (V_inf_initial/np.linalg.norm(V_inf_initial))
    V_inf_vec = V_inf_initial - DV_electric_vec
    V_inf = np.linalg.norm(V_inf_vec)
    #print(V_inf)

    Vp_hyp = get_Vp_hyp(V_inf, mu2, rp)
    Vp_cap = get_Vp_cap(mu2, e_cap, rp)

    DV_insertion = Vp_hyp - Vp_cap
    #print(DV_insertion)
    return DV_insertion, V_inf, np.linalg.norm(V_inf_initial)

def straight_line(x, m, c): return m*x + c

###########################################################################
#GET THE INSERTION DV

#Define some things

e_cap = 0.
rp = 1588.                      #km
mu2 = 0.87E3#8.719*10**11*10**-9      #km^3/2^2
DV_electric = 5.736              #km/s
Va_initial = np.array([6.814488756803494,               0.6566676716964039,              -1.123373303125703])
V_pluto = np.array([2.294204405818791,         3.790566817555419,        0.5042927954002003])


Va_geoff = np.array([7.02,0.4986,-2.11])
Vp_geoff = np.array([2.243,3.83,0.518])
Va_geoff_norm = Va_geoff/np.linalg.norm(Va_geoff)
#Va_final_geoff = Va_geoff - Va_geoff_norm*DV_electric
V_inf_geoff = -Vp_geoff + Va_geoff

DV_insertion1 = getDvInsertion(Va_initial, V_pluto, DV_electric, mu2, rp, e_cap)
DV_insertion_geoff = getDvInsertion(Va_geoff, Vp_geoff, DV_electric, mu2, rp, e_cap)

eccs = np.arange(0,1, 0.001)
DV_insertions = []
straight_dvs = []
for ecc in eccs:
    DV_insertion = getDvInsertion(Va_initial, V_pluto, DV_electric, mu2, rp, ecc)
    DV_insertions.append(DV_insertion[0])

m = -DV_insertions[0]
c = DV_insertions[0]
for ecc in eccs:
    straight_dvs.append(straight_line(ecc, m, c))

plt.plot(eccs, DV_insertions)
plt.plot(eccs, straight_dvs)
plt.xlabel('ecc')
plt.ylabel('DV_insertion')
plt.grid()
plt.show()
