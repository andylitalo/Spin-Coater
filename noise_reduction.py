# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 16:52:43 2017

@author: Andy
"""

## -*- coding: utf-8 -*-

import glob
import matplotlib.pyplot as plt
import cPickle as pkl
import numpy as np
from scipy import signal
import csv
import numpy as np
from scipy.signal import butter, lfilter, freqz
import matplotlib.pyplot as plt
import copy
import Functions as Fun
import SignalProcessingFunctions as SPF


# User parameters
# Data parameters
folder = 'Code from John 10-14\\'
saveFolder = 'noise_reduction\\'
name = 'all_data.pkl'
# Condition parameters
QList = np.array([500,1000,1500,2000,2500,3000,3250])
RPMList =  np.array([0,30,120,250,500,750,1000])
conditionList = ['Dry'] #['Dry','Water', '3mM SDS', '6.5mM SDS', '10mM SDS', 'Water with 10mM SDS'] 
# Plot parameters
plotRDiff = True
# Save Parameters
savePlots = True
saveRCritData = False
saveRCritName = 'r_crit.pkl'
# EP and Disk Parameters
diskFrac = 0.9
ePNoise = 0.2 # noise in eP value, approximately - used to determine if rivulets were formed
frac_eP = 0.1 # assign r_crit when this fraction (above 1) of max EP reached

# Derived Quantities and Constants
diskRad = 15 # cm
radiusNoisy = diskRad / 4
wetTarget = diskFrac * diskRad

# Cleanup existing windows
plt.close('all')
    
# Get list of data files to process
fileList = glob.glob(folder+name)
N1 = len(fileList)

# Initialize r_crit data structure
r_crit_data = {}
r_crit_mat = np.zeros([len(QList), len(RPMList)])
for i in range(len(conditionList)):
    r_crit_data[conditionList[i]] = r_crit_mat

for j in range(N1):
    
    # Load data and get list of fields
    if 'allData' not in globals():
        with open(fileList[j],'rb') as f:
            allData = pkl.load(f)
            print fileList[j]
    keyList = allData.keys()
    N = len(keyList)
    
    # Initialize data storage arrays
    plotData = {}
    for condition in conditionList:
        plotData[condition] = {}
        plotData[condition]['ePMax'] = np.zeros((len(QList),len(RPMList)))
        plotData[condition]['wetTime'] = np.zeros((len(QList),len(RPMList)))
        plotData[condition]['wetVolume'] = np.zeros((len(QList),len(RPMList)))

    
    for i in range(N):
    
        # Parse the experiment data
        data = allData[keyList[i]]
        RPM = data['RPM']
        QTarget = data['flowRate']
        condition = data['condition']
        if condition not in conditionList:
            continue
        if (QTarget not in QList) or (RPM not in RPMList):
            continue
        
        # Experimental Conditions
        print 'RPM = ' + str(RPM)
        print 'Q = ' + str(QTarget)
        print 'Condition = ' + condition
        
        # Set the indicies of the array for assignment of plotting data
        i_Q = np.argmax(QList==QTarget)
        j_RPM = np.argmax(RPMList==RPM)

        # Parse remaining experiment data
        fps = data['fps']
        t0 = data['t0']
        time = data['time']
        aMax = data['aMax']
        aMin2 = data['aMin2']
        aMean = data['aMean']
        eP = data['excessPerimeter']
        
        if len(time) < 1:
            continue
        
        # Adjust data fields
        timeRinsing = time[time > t0] - t0
        aMax = aMax[time > t0]
        aMin2 = aMin2[time > t0]
        aMean = aMean[time > t0]
        eP = eP[time > t0]
        Q = Fun.convert_flowrate(QTarget)
        
        #######################################################################
        ### REDUCE NOISE ###
        #######################################################################
            
        # Filter parameters
            ###### ADJUST WITH FPS ######
        unscaledRPM = 250
        timeWindow = 0.01 # sec
        if RPM == 0:
            continue
        nn_MF = int(timeWindow*fps*unscaledRPM/RPM) # number of nearest neighbors in filtering window
        nn_LA = int(timeWindow*fps*unscaledRPM/RPM) # number of nearest neighbors in filtering window
        
        # Filters
        eP_CircleStart = eP
        eP_CircleStart[aMean < radiusNoisy] = 1
        eP_MF = signal.medfilt(eP_CircleStart, 2*nn_MF + 1)
        eP_LA = SPF.local_average_filter(eP_MF, nn_LA)     
        
        eP_filtered = eP_LA
        
        #######################################################################
        ##################### FIND CRITICAL RADIUS ############################
        
        ePMax = max(eP_filtered)
        if ePMax < 1 + ePNoise:
            print 'No critical radius - ePMax threshold not reached'
            continue
        ind_ePMax = [i for i in range(len(eP_filtered)) if eP_filtered[i] == ePMax][0]
        if ind_ePMax == 0:
            print 'No critical radius - ePMax reached immediately'
            continue
        
        ############# Limit ePDiff to EP's <= 1/2 max EP ######################
        inds_ePHalfMax = [i for i in range(ind_ePMax) if eP_filtered[i] <= ((ePMax-1)/2 + 1)]
        ind_last_ePHalfMax = inds_ePHalfMax[len(inds_ePHalfMax) - 1]
        ePDiff = eP_filtered[1:ind_last_ePHalfMax] - eP_filtered[0:ind_last_ePHalfMax-1]    
        
        # Find last index where EP decreased before ePMax - critical radius reached at next index
        inds_neg = [i for i in range(len(ePDiff)) if ePDiff[i] < 0]
        if len(inds_neg) == 0:
            print 'No critical radius - no negative index before max EP'
            continue
        ind_last_neg = inds_neg[len(inds_neg)-1]
        
        # Find last index where EP < some fraction of max EP
        inds_frac = [i for i in range(len(ePDiff)) if eP_filtered[i] < frac_eP*(ePMax - 1) + 1]
        if len(inds_frac) == 0:
            print 'No critical radius - EP never exceeds ' + frac_eP + ' of max EP'
            continue
        ind_last_frac = inds_frac[len(inds_frac)-1]
        ind_r_crit = max(ind_last_neg + 1, ind_last_frac)
        t_r_crit = time[ind_r_crit]
        print 't_r_crit = ' + str(t_r_crit)
        r_crit = aMean[ind_r_crit + (len(time) - len(aMean))]
        print 'r_crit = ' + str(r_crit)
        
        # Save r_crit
        r_crit_data[condition][i_Q, j_RPM] = r_crit
        
        # PLOT
        plt.plot(time, eP_filtered, 'g-', linewidth=3, label='filtered data')
        plt.plot(time, eP, 'k.', markersize=3, label='raw data')
        if plotRDiff:
            plt.plot(time, aMax-aMin2, 'r.', markersize=5, label='aMax - aMin2')
        plt.plot((t_r_crit, t_r_crit), (1, ePMax), 'r--', label='critical point')
        plt.xlabel('Time [sec]')
        plt.grid()
        plt.legend(loc=0)
        pltTitle = 'Finding Critical Radius - RPM = %i, Q = %i' %(RPM, int(round(Q,-1)))
        plt.title(pltTitle)
        if savePlots:
            savePlotName = Fun.get_save_name(saveFolder, RPM, QTarget, condition, '')
            plt.savefig(savePlotName, bbox_inches='tight')
        plt.show()
        plt.close()
 
if saveRCritData:
    r_crit_data['QList'] = QList
    r_crit_data['RPMList'] = RPMList
    r_crit_data['conditionList'] = conditionList       
    with open(saveFolder + saveRCritName,'wb') as f:
        pkl.dump(r_crit_data,f)