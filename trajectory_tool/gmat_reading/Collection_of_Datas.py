import numpy as np
import os
from trajectory_tool.Insertion_Rough_Dv import getDvInsertion, get_Vp_cap, get_Vp_hyp, Va_initial, V_pluto, DV_electric, mu2, rp
import json

#PLUTO MAINTENANCE OPTIONS
pluto_direc = 'C:\\Users\\matth\\AppData\\Local\\GMAT\\R2017a\\output\\Maintenance files\\Pluto_Only\\just_maintenances'
fnames = os.listdir(pluto_direc)

data={}
for fname in fnames:
    some_data = np.loadtxt(os.path.join(pluto_direc, fname))
    DV_per_month = some_data[-3]
    data[fname[-6:-4]] = DV_per_month

PS = data
print('PLuto shizzle: ', data)

# np.savetxt(fname_ps, PS)
newdirec = 'C:\\Users\\matth\\AppData\\Local\\GMAT\\R2017a\\output\\Maintenance files\\Pluto_Only'
fname_save = os.path.join(newdirec, 'pluto_maintain')
with open(fname_save, 'w') as fp:
    json.dump(PS, fp)


#CHARON MAINTENANCE OPTIONS
charon_direc = 'C:\\Users\\matth\\AppData\\Local\\GMAT\\R2017a\\output\\Maintenance files\\Charon_Only\\just_maintenances'
fnames = os.listdir(charon_direc)

data={}
for fname in fnames:
    some_data = np.loadtxt(os.path.join(charon_direc, fname))
    DV_per_month = some_data[-3]
    data[fname[-6:-4]] = DV_per_month

CS = data
print('Charon shizzle: ', data)
newdirec = 'C:\\Users\\matth\\AppData\\Local\\GMAT\\R2017a\\output\\Maintenance files\\Charon_Only'
fname_save = os.path.join(newdirec, 'charon_maintain')
with open(fname_save, 'w') as fp:
    json.dump(CS, fp)

#TRANSFER ORBIT OPTIONS

transfer_direc = 'C:\\Users\\matth\\AppData\\Local\\GMAT\\R2017a\\output\\Maintenance files\\Transfer\\JustDVs'
fnames = os.listdir(transfer_direc)

data={}
for fname in fnames:
    filepath = os.path.join(transfer_direc, fname)
    all_data = np.genfromtxt(filepath, skip_footer=2)
    mag_data = np.genfromtxt(filepath, skip_header=7)
    TOI_DV, COI_DV = mag_data[0], mag_data[1]
    DV_tot = TOI_DV + COI_DV
    DVs = [TOI_DV, COI_DV, DV_tot]
    data[fname[-6:-4]] = DVs

transfer_dvs = data

#Left = transfer case, right = Pluto-Charon Combo of cases
transfer_case_dic = {'V2':'4-2', 'V3':'7-2', 'V4':'8-2', 'V5':'4-6', 'V6':'7-6', 'V7':'8-6'}

pluto_months = 6
charon_months = 6

pluto_ecc_incs = {'V2': [0.16, 60], 'V3': [0.16, 90], 'V4': [0.18, 77], 'V5': [0.09, 79], 'V6': [0.09, 90], 'V7': [0.04, 85], 'V8': [0.31, 45], 'V9': [0.25, 60]}
charon_ecc_incs = {'V2': [0.12, 65], 'V3': [0.12, 75], 'V4': [0.1, 75], 'V5': [0.1, 90], 'V6': [0.05, 65], 'V7': [0.05, 90]}


allDVs = {}
for case in transfer_case_dic:
    #orbit DVs
    combo = transfer_case_dic[case]
    pluto_orbit = 'V'+combo[0]
    charon_orbit = 'V'+combo[-1]

    #maintenance DVs
    pluto_maintenance = pluto_months* PS[pluto_orbit]
    charon_maintenance = charon_months* CS[charon_orbit]

    #Insertion DV
    DV_insertion = getDvInsertion(Va_initial, V_pluto, DV_electric, mu2, rp, pluto_ecc_incs[pluto_orbit][0])[0]

    totalDV=transfer_dvs[case][-1]*1000 + pluto_maintenance + charon_maintenance + DV_insertion*1000
    caseDV = {'Transfer DVs [TOI, COI, total]': np.array(transfer_dvs[case])*1000, 'Pluto maintenance': pluto_maintenance, 'Charon maintenance': charon_maintenance, 'Insertion': DV_insertion*1000, 'Total': totalDV}
    # has [total transfer DV, TOI DV, COI DV, pluto maintenance DV, charon maintenance DV, orbit insertion DV, total of all DV]
    allDVs[transfer_case_dic[case]] = caseDV

best_case = allDVs['4-6']
print(best_case)

newdirec = 'C:\\Users\\matth\\AppData\\Local\\GMAT\\R2017a\\output\\Maintenance files\\Transfer'
fname_save = os.path.join(newdirec, 'transfer_summary')
with open(fname_save, 'w') as fp:
    fp.write(str(allDVs))

