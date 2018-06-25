import numpy as np
import os
from trajectory_tool.Insertion_Rough_Dv import getDvInsertion, get_Vp_cap, get_Vp_hyp, Va_initial, V_pluto, DV_electric, mu2, rp
import json

from trajectory_tool.gmat_reading.Collection_of_Datas import charon_ecc_incs

mu_p = 8.71E+02         #mu of pluto
mu_c = 1.06E+02         #mu of charon
a_c = 1.7536E+04          #SMA of charon around pluto
r_sc_c_p = 1.00E+03     #orbital radius of periapsis of sc around charon

def get_escape_dv(mu_p, mu_c, a_c, r_sc_c_p, e_c ):
    SMA = r_sc_c_p/(1-e_c)  #get charon SMA

    V_p_c = np.sqrt(2*(mu_c/r_sc_c_p) - mu_c/(SMA))   #get peri V
    #V_esc_c = np.sqrt(2*mu_c/r_sc_c_p)                  #get escape v around charon

    V_c = np.sqrt(2*(mu_p/a_c) - mu_p/(a_c))          #Charon v around pluto
    V_esc_p = np.sqrt(2*mu_p/a_c)                       #escape velocity of P-C at Charo alt
    DV_esc_p = V_esc_p-V_c                              # DV to escape from P-C at Charon Alt

    a_esc_inf = -mu_c/DV_esc_p**2                       # SMA at charon to achieve DV_esc_p

    V_esc_inf = np.sqrt(2*mu_c/SMA - mu_c/a_esc_inf)    #Velocity at peri to escape

    DV_sustain = V_esc_inf - V_p_c                      #DV needed for sustainability burn
    #print(SMA, V_p_c, V_c, V_esc_p, DV_esc_p, a_esc_inf, V_esc_inf)
    return DV_sustain*1000
charon_ecc_incs_escapes = charon_ecc_incs
for case in charon_ecc_incs_escapes:
    ecc = charon_ecc_incs_escapes[case][0]
    DV_escape = get_escape_dv(mu_p, mu_c, a_c, r_sc_c_p, ecc)
    charon_ecc_incs_escapes[case].append(DV_escape)

newdirec = 'C:\\Users\\matth\\Documents\\realdocs_work\\delft\\schoolwork\\Year 3\\Q4 (DSE)\\final\\Maintenance files\\Escape'
fname_save = os.path.join(newdirec, 'escape_summary.txt')
with open(fname_save, 'w') as fp:
    fp.write('case: ecc \t inc \t escape dv(m/s) \n \n')
    fp.write(str(charon_ecc_incs_escapes))
