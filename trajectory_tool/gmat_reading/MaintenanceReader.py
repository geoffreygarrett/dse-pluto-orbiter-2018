import numpy as np
import os
location = 'C:\\Users\\matth\\AppData\\Local\\GMAT\\R2017a\\output'
body = 'Charon'
version = '9'
endnameread = 'OrbitMaintenance'
endnamewrite = 'OrbitMaintenance_DVs_V'+version
ext='.txt'
fname = os.path.join(location, body + endnameread + ext)
alldata = np.loadtxt(fname, comments='O')

periods=[]
apos=[]
peris=[]
Cmaintain_ones=[]
incs=[]
Cmaintain_twos=[]
julians=[]

number=-1
for i in alldata:
    number+=1
    periods.append(i[0])
    apos.append(i[1])
    peris.append(i[2])
    Cmaintain_ones.append(i[3])
    incs.append(i[4])
    Cmaintain_twos.append(i[5])
    julians.append(i[6])

time_taken = julians[-1] - julians[0]
av_period_days = np.mean(periods)/60/60/24

Cmaintain_ones_abs = [abs(x) for x in Cmaintain_ones]
Cmaintain_twos_abs = [abs(x) for x in Cmaintain_twos]

av_Cmaintain_ones = np.mean(Cmaintain_ones_abs)
av_Cmaintain_twos = np.mean(Cmaintain_twos_abs)
av_Cmaintains = np.mean([av_Cmaintain_ones, av_Cmaintain_twos])
dv_over_time = (sum(Cmaintain_ones_abs) + sum(Cmaintain_twos_abs))*1000
dv_per_month = (30/time_taken)*dv_over_time
charon_time = 6
total_dv = dv_per_month*charon_time
all_Cmaintains = Cmaintain_ones + Cmaintain_twos
max_man = max([abs(x) for x in all_Cmaintains])
min_man = min([abs(x) for x in all_Cmaintains])


print('THINGS AT CHARON')
print('Over the time (days): ', time_taken)
print('Average period (days): ', av_period_days)
print('Average maintenance per orbit for ones (m/s): ', av_Cmaintain_ones*1000)
print('Average maintenance per orbit for twos (m/s): ', av_Cmaintain_twos*1000)
print('Maximum maneuver size (m/s): ', max_man*1000)
print('Minimum maneuver size (m/s): ', min_man*1000)
print('Average maintenance per orbit total (m/s): ', av_Cmaintains*1000)
print('Total maintenance per month: ', dv_per_month)
print('Total maintenance over ', charon_time, ' months (m/s): ', total_dv)
print('-                                             -')

txt_header = '[Simulation time, \naverage orbital period(days), \naverage maneuver magnitude (m/s), \nmax maneuever mag (m/s), \nmin maneuever mag (m/s), \nDV per month (m/s), \nmaintenance time for total DV (months), \ntotal DV (m/s)]\n'
fname_dvs = os.path.join(location, body + endnamewrite + ext)
DVs_save = np.array([time_taken,
                     av_period_days,
                     av_Cmaintains,
                     max_man*1000,
                     min_man*1000,
                     dv_per_month,
                     charon_time,
                     total_dv])

np.savetxt(os.path.join(location, fname_dvs), DVs_save, header=txt_header)

#########################################################################################################################
#########################################################################################################################
#########################################################################################################################

body = 'Pluto'
fname = os.path.join(location, body + endnameread + ext)
alldata = np.loadtxt(fname, comments='O')

periods=[]
apos=[]
peris=[]
Cmaintain_ones=[]
incs=[]
Cmaintain_twos=[]
julians=[]

number=-1
for i in alldata:
    number+=1
    periods.append(i[0])
    apos.append(i[1])
    peris.append(i[2])
    Cmaintain_ones.append(i[3])
    incs.append(i[4])
    Cmaintain_twos.append(i[5])
    julians.append(i[6])

time_taken = julians[-1] - julians[0]
av_period_days = np.mean(periods)/60/60/24

Cmaintain_ones_abs = [abs(x) for x in Cmaintain_ones]
Cmaintain_twos_abs = [abs(x) for x in Cmaintain_twos]

av_Cmaintain_ones = np.mean(Cmaintain_ones_abs)
av_Cmaintain_twos = np.mean(Cmaintain_twos_abs)
av_Cmaintains = np.mean([av_Cmaintain_ones, av_Cmaintain_twos])
dv_over_time = (sum(Cmaintain_ones_abs) + sum(Cmaintain_twos_abs))*1000
dv_per_month = (30/time_taken)*dv_over_time
charon_time = 6
total_dv = dv_per_month*charon_time
all_Cmaintains = Cmaintain_ones + Cmaintain_twos
max_man = max([abs(x) for x in all_Cmaintains])
min_man = min([abs(x) for x in all_Cmaintains])


print('THINGS AT PLUTO')
print('Over the time (days): ', time_taken)
print('Average period (days): ', av_period_days)
print('Average maintenance per orbit for ones (m/s): ', av_Cmaintain_ones*1000)
print('Average maintenance per orbit for twos (m/s): ', av_Cmaintain_twos*1000)
print('Maximum maneuver size (m/s): ', max_man*1000)
print('Minimum maneuver size (m/s): ', min_man*1000)
print('Average maintenance per orbit total (m/s): ', av_Cmaintains*1000)
print('Total maintenance per month: ', dv_per_month)
print('Total maintenance over ', charon_time, ' months (m/s): ', total_dv)
print('-                                                                             -')

fname_dvs = os.path.join(location, body + endnamewrite + ext)
DVs_save = np.array([time_taken,
                     av_period_days,
                     av_Cmaintains,
                     max_man*1000,
                     min_man*1000,
                     dv_per_month,
                     charon_time,
                     total_dv])

np.savetxt(os.path.join(location, fname_dvs), DVs_save, header=txt_header)





