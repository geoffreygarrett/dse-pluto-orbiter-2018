import numpy as np
import os
from shutil import copyfile, copy
import json
import matplotlib.pyplot as plt

def domaintenance(base_path_source, base_path_tosave, body, delta):
    location=base_path_source
#location = 'C:\\Users\\matth\\AppData\\Local\\GMAT\\R2017a\\output'
    #body = 'Charon'
    endnameread = 'OrbitMaintenance'
    endnamewrite = 'OrbitMaintenance_DVs_'+delta
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


    print('THINGS AT ' + body)
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
    fname_dvs = os.path.join(base_path_tosave, body + endnamewrite + ext)
    DVs_save = np.array([time_taken,
                     av_period_days,
                     av_Cmaintains,
                     max_man*1000,
                     min_man*1000,
                     dv_per_month,
                     charon_time,
                     total_dv])
    print(os.path.join(base_path_tosave, fname_dvs))
    np.savetxt(os.path.join(base_path_tosave, fname_dvs), DVs_save, header=txt_header)

    fname_raw = 'OrbitMaintenance_raw_' + delta
    savepath_raw = os.path.join(base_path_tosave, body+ fname_raw+ ext)
    copyfile(fname, savepath_raw)

    fname_kep_src = 'KeplerianIintial' + body + ext
    fname_kep_drc = body +'KeplerianIintial'  +'_'+ delta + ext
    kep_src = os.path.join(location, fname_kep_src )
    kep_drc = os.path.join(base_path_tosave, fname_kep_drc)
    copyfile(kep_src, kep_drc)


    #np.savetxt(savepath_raw, alldata_raw)
    return DVs_save

body = 'Pluto'
ecc_inc = 'inc'
source = 'C:\\Users\\matth\\AppData\\Local\\GMAT\\R2017a\\output'
saveto = os.path.join('C:\\Users\\matth\\Documents\\realdocs_work\\delft\\schoolwork\\Year 3\\Q4 (DSE)\\final\\Sensitivity\\Maintenance', body, ecc_inc)

#uncomment below to save things
#this = domaintenance(source, saveto, body, 'iM20')


##############################################################################################################################################################
#plotting area

def sortdata(body, ecc_inc):

    #body = 'Pluto'
    #ecc_inc = 'ecc'
    direc = os.path.join('C:\\Users\\matth\\Documents\\realdocs_work\\delft\\schoolwork\\Year 3\\Q4 (DSE)\\final\\Sensitivity\\Maintenance', body, ecc_inc)

    fnames = os.listdir(direc)
    if 'summary.txt' in fnames: fnames.remove('summary.txt')

    dv_fnames=[]
    for fname in fnames:
        if 'OrbitMaintenance_DVs' in fname:
            dv_fnames.append(fname)


    MaintenanceDatas = []
    MaintenanceCases = []
    toappend =[]
    for file in dv_fnames:
        if file[-5]=='5' and file[-6] !='1': case = file[-7:-4]
        else: case = file[-8:-4]
        if file[0] == 'P': case_full = 'P-'+case
        else: _full = 'C-'+case
        #
        # DV_data = np.loadtxt(os.path.join(direc, file), comments='#')
        # #Orbit_data = np.loadtxt(os.path.join(direc, ))
        # DV_per_month = DV_data[-3]
        #MaintenanceDatas.append(data)
        MaintenanceCases.append(case)

    #creates nested list where each element is a case, each subelement is it's file components
    splitfiles =[]
    for case in MaintenanceCases:
        list = []
        for file in fnames:
            if case in file: list.append(file)
        splitfiles.append(list)

    #saves to a text file
    MaintenanceHeaders = ['Case', 'e', 'i', 'DV_Per_Month']
    alldata=[]
    with open(os.path.join(direc, 'summary.txt'), 'w') as fp:
        fp.write(str(MaintenanceHeaders).replace('[','').replace(']','').replace("'",''))
        for element in splitfiles:
            case = body[0] + '-' + element[0].split('_')[-1].split('.')[0]
            subdata=[case]
            for file in element:
                if body + 'KeplerianIintial' in file:
                    kep = np.loadtxt(os.path.join(direc,file), dtype=str)
                    ecc, inc = kep[1][2],kep[1][3]
                    if body == 'Charon': inc = kep[1][-1]
                    subdata.append(float(ecc))
                    subdata.append(float(inc))
                elif body + 'OrbitMaintenance_DVs' in file:
                    DV_data = np.loadtxt(os.path.join(direc, file))
                    dv_per_month = DV_data[-3]
                    subdata.append(dv_per_month)
            alldata.append(subdata)

            case_params = [case, ecc, inc, dv_per_month]

            towrite = str(case_params).replace('[','').replace(']','').replace("'",'')
            fp.write('\n' + towrite)
    alldata = np.array(alldata)
    return alldata

bodies = ['Charon', 'Pluto']
ecc_incs = ['ecc', 'inc']


for body in bodies:
    for ecc_inc in ecc_incs:
        a_data = sortdata(body, ecc_inc)
        if body == bodies[0] and ecc_inc == ecc_incs[0]: C_e_data=a_data
        if body == bodies[0] and ecc_inc == ecc_incs[1]: C_i_data=a_data
        if body == bodies[1] and ecc_inc == ecc_incs[0]: P_e_data=a_data
        if body == bodies[1] and ecc_inc == ecc_incs[1]: P_i_data=a_data


def plotit(data_sum, variable):
    #data_sum = P_e_data
    if variable == 'ecc': col = 1
    else: col = 2
    data_sum = data_sum[data_sum[:,col].argsort()]
    #data_sum = np.ndarray.sort(data_sum, axis=-1)
    variable_no = data_sum[:,col]
    dvs = data_sum[:,-1]

    plt.figure(data_sum[0][0][0:3])
    if variable=='ecc': plt.xlabel('Eccentricity')
    else: plt.xlabel('Inclination [' + (u'\N{DEGREE SIGN}') + ']')
    # start = float(variable_no[0])
    # stop = float(variable_no[-1])
    # step = (stop-start)/1000
    # plt.xticks(np.arange(start, stop, step=step))
    variable_new=[]
    for i in variable_no: variable_new.append(float(i))
    dvs_new=[]
    for i in dvs: dvs_new.append(float(i))
    plt.ylabel(u'\N{GREEK CAPITAL LETTER DELTA}'+'V per month [m/s]')
    plt.plot(variable_new, dvs_new)

# plotit(P_i_data, 'inc')
# plotit(P_e_data, 'ecc')
# plotit(C_i_data, 'inc')
# plotit(C_e_data, 'ecc')

###################################################################################################
# Transfer section                                                                                #
###################################################################################################


def dotransfer(base_path_source, base_path_tosave, body, delta):
    location=base_path_source
#location = 'C:\\Users\\matth\\AppData\\Local\\GMAT\\R2017a\\output'
    #body = 'Charon'
    delta = '_'+delta
    endnameread = 'TransferDV'
    endnamewrite = 'TransferDV'+delta
    ext='.txt'
    fname = os.path.join(location, endnameread + ext)
    DVs = np.loadtxt(fname, skiprows=7)
    print('TOI DV = ', DVs[0])
    print('COI DV = ', DVs[1])

    savepath = os.path.join(base_path_tosave, endnamewrite+ext)
    copyfile(fname, savepath)

    endnameread_pluto = 'PlutoKeplerianInitial'
    readpath_pluto = os.path.join(base_path_source, endnameread_pluto+ext)
    savepath_pluto = os.path.join(base_path_tosave, endnameread_pluto + delta + ext)
    copyfile(readpath_pluto, savepath_pluto)

    endnameread_charon = 'CharonKeplerianFinal'
    readpath_charon = os.path.join(base_path_source, endnameread_charon + ext)
    savepath_charon = os.path.join(base_path_tosave, endnameread_charon + delta + ext)
    copyfile(readpath_charon, savepath_charon)



    #np.savetxt(savepath_raw, alldata_raw)
    #return DVs_save

body = 'Pluto'
ecc_inc = 'ecc'
source = 'C:\\Users\\matth\\AppData\\Local\\GMAT\\R2017a\\output'
saveto = os.path.join('C:\\Users\\matth\\Documents\\realdocs_work\\delft\\schoolwork\\Year 3\\Q4 (DSE)\\final\\Sensitivity\\Transfer', body, ecc_inc)

#dotransfer(source, saveto, body, 'eP15')


########################################################################################
#Insertion sensitivity                                                                 #
########################################################################################

