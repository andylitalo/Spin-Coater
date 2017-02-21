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
import Functions as Fun
import SignalProcessingFunctions as SPF

# Data parameters
dataFolder = 'Data\\'
dataFile = 'all_data.pkl'
# Save Parameters
savePlots = False
saveRCritData = False
saveFolder = 'locate_r_Crit\\'
saveRCritName = 'r_crit_splash.pkl'
# Condition parameters
QList = np.array([500,1000,1500,2000,2500,3000,3250])
RPMList =  np.array([0,30,120,250,500,750,1000])
conditionList = ['Dry','Water', '3mM SDS', '6.5mM SDS', '10mM SDS', 'Water with 10mM SDS'] 
# MinSep Analysis Parameters
minSepCondition = False # require r_crit reached if aMax - aMin2 > minSep
minSep = 1.7 # [cm] critical value of aMax - aMin2 that sets upper limit on 
minSepFracCondition = False 
minSepFrac = 0.2 # fraction of perturbation to contact line indicating rCrit
# rSplash Analysis Parameters
rSplashCondition = True
rSplashFile = 'r_splash.pkl'
# Plot parameters
showPlots = True
plotSep = True # plot aMax - aMin2
legend_fs = 8 # font size of legend
# EP and Disk Parameters
diskFrac = 0.9 # fraction of disk to consider in analysis
ePNoise = 0.2 # noise in eP value, approximately - used to determine if rivulets were formed
fracEP = 0.1 # assign rCrit when this fraction (above 1) of max EP reached

# Derived Quantities and Constants
diskRad = 15.0 # cm
radiusNoisy = diskRad / 4.0
wetTarget = diskFrac * diskRad


##############################################################################


# Cleanup existing windows
plt.close('all')
    
# Get list of data files to process
fileList = glob.glob(dataFolder + dataFile)
N1 = len(fileList)

# Initialize rCrit data structure
rCritData = {}   

# Analyze data
for i in range(N1):
    
    # Load data if not already loaded
    if 'allData' not in globals():
        with open(fileList[i],'rb') as f:
            allData = pkl.load(f)
    # Get data keys
    keyList = allData.keys()
    N = len(keyList)
    
    # Load rSplash Data    
    if rSplashCondition:
        with open(dataFolder + rSplashFile, 'rb') as f:
            rSplashData = pkl.load(f)
        splashQList = rSplashData['QList']

    # Initialize data storage arrays
    plotData = {}
    for condition in conditionList:
        plotData[condition] = {}
        plotData[condition]['ePMax'] = np.zeros((len(QList),len(RPMList)))
        plotData[condition]['wetTime'] = np.zeros((len(QList),len(RPMList)))
        plotData[condition]['wetVolume'] = np.zeros((len(QList),len(RPMList)))
        rCritData[condition] = np.zeros([len(QList), len(RPMList)])
    
    for j in range(N):
    
        # Parse the experiment data
        data = allData[keyList[j]]
        RPM = data['RPM']
        QTarget = data['flowRate']
        condition = data['condition']
        if condition not in conditionList:
            continue
        if (QTarget not in QList) or (RPM not in RPMList):
            continue
        
        # Experimental Conditions
        print (RPM, QTarget, condition)
        
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
        
        # Skip if empty experiment
        if len(time) == 0:
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
        
        # Skip stationary experiments since no rCrit        
        if RPM == 0:
            continue
        
        # Filter parameters
        unscaledRPM = 250 
        timeWindow = 0.02 # sec
        nn_MF = int(timeWindow*fps*unscaledRPM/RPM) # number of nearest neighbors in filtering window
        nn_LA = int(timeWindow*fps*unscaledRPM/RPM) # number of nearest neighbors in filtering window
        
        # Filters
        eP_CircleStart = eP
        eP_CircleStart[aMean < radiusNoisy] = 1 # set initial EP to 1
        eP_MF = signal.medfilt(eP_CircleStart, 2*nn_MF + 1) # median filter
        eP_LA = SPF.local_average_filter(eP_MF, nn_LA)
        
        ePFiltered = eP_LA
        
        #######################################################################
        ##################### FIND CRITICAL RADIUS ############################
        
        # Locate max EP
        ePMax = max(ePFiltered)
        if ePMax < 1 + ePNoise:
            print 'No critical radius - ePMax threshold not reached'
            continue
        indEPMax = [k for k in range(len(ePFiltered)) if ePFiltered[k] == ePMax][0]
        if indEPMax == 0:
            print 'No critical radius - ePMax reached immediately'
            continue
        
        minInd = 0
        # Locate index of splash radius
        if rSplashCondition:
            indQ = [k for k in range(len(splashQList)) if splashQList[k] == QTarget]
            rSplash = rSplashData[condition][indQ]
            minInd = [k for k in range(len(aMean)-1) if (aMean[k+1]>=rSplash and aMean[k] <= rSplash)]
            if not minInd:
                minInd = 0
            else:
                minInd = minInd[0]
            tRSplash = time[minInd]
        
        # Calculate local difference in EP - limit to EP's <= 1/2 max EP
        indsEPHalfMax = [k for k in range(indEPMax) if ePFiltered[k] <= ((ePMax-1)/2 + 1)]
        indLastEPHalfMax = indsEPHalfMax[len(indsEPHalfMax) - 1]
        ePDiff = ePFiltered[1:indLastEPHalfMax] - ePFiltered[0:indLastEPHalfMax-1]    
        
        # Find last index where EP decreased before ePMax - critical radius reached at next index
        indsNeg = [k for k in range(minInd, len(ePDiff)) if ePDiff[k] < 0]
        if len(indsNeg) == 0:
            print 'No critical radius - no negative index before max EP'
            continue
        indLastNeg = indsNeg[len(indsNeg)-1]
        
        # Find last index where EP < some fraction of max EP
        indsFrac = [k for k in range(len(ePDiff)) if ePFiltered[k] < fracEP*(ePMax - 1) + 1]
        if len(indsFrac) == 0:
            print 'No critical radius - EP never exceeds ' + str(fracEP) + ' of max EP'
            continue
        indLastFrac = indsFrac[len(indsFrac)-1]
        
        # rCrit at greater of rCrit found through last negative and fractional approaches
        indRCrit = max(indLastNeg + 1, indLastFrac)
        
        ################# TRIAL ###############################################
        # Ensure separation between aMax and aMin2 isn't too big before rCrit
        # Adjust rCrit accordingly
        if minSepCondition:
            sep = aMax - aMin2
            indsMinSep = [k for k in range(minInd, len(sep[:indEPMax])) if sep[k] <= minSep]
            indLastMinSep = indsMinSep[len(indsMinSep) - 1]
            indRCrit = min(indLastMinSep, indRCrit)
            if indRCrit == indLastMinSep:
                print 'rCrit updated with minSep heuristic'
                
        if minSepFracCondition:
            sep = aMax - aMin2
            sep[aMean < radiusNoisy] = 0
            sepMax = max(sep)
            indSepMax = [k for k in range(len(sep)) if sep[k] == sepMax][0]
            indsMinSepFrac = [k for k in range(minInd, len(sep[:indSepMax])) if sep[k] <= minSepFrac*aMean[k]]
            if len(indsMinSepFrac) != 0:            
                indLastMinSepFrac = indsMinSepFrac[len(indsMinSepFrac) - 1]
                indRCrit = min(indLastMinSepFrac, indRCrit)
                if indRCrit == indLastMinSepFrac:
                    print 'rCrit updated with minSepFrac heuristic'
        #######################################################################
        # Print results
        tRCrit = time[indRCrit]
        print 'tRCrit = ' + str(tRCrit)
        rCrit = aMean[indRCrit + (len(time) - len(aMean))]
        print 'rCrit = ' + str(rCrit)
        
        # Save rCrit
        rCritData[condition][i_Q, j_RPM] = rCrit
        
        # Plot
        if showPlots:
            plt.plot(time, ePFiltered, 'g-', linewidth=3, label='filtered data')
            plt.plot(time, eP, 'k.', markersize=3, label='raw data')
            if plotSep:
                plt.plot(time, aMax - aMin2, 'r.', markersize=5, label='aMax - aMin2')
            # Show critical point
            axes = plt.gca()
            yLim = axes.get_ylim()
            plt.plot((tRCrit, tRCrit), yLim, 'r--', label='critical point')
            if rSplashCondition:
                plt.plot((tRSplash, tRSplash), yLim, 'b--', label='splash radius')
            # Formatting
            plt.xlabel('Time [sec]')
            plt.ylabel('EP / aMax - aMin2 [cm]')
            plt.grid()
            plt.legend(loc=0, fontsize=legend_fs)
            pltTitle = 'Finding Critical Radius - RPM = %i, Q = %i, %s' %(RPM, int(round(Q,-1)), condition)
            plt.title(pltTitle)
            if savePlots:
                savePlotName = Fun.get_save_name(saveFolder, RPM, QTarget, condition, '')
                plt.savefig(savePlotName, bbox_inches='tight')
            plt.show()
            plt.close()
 
# Saving rCrit
if saveRCritData:
    rCritData['QList'] = QList
    rCritData['RPMList'] = RPMList
    rCritData['conditionList'] = conditionList       
    with open(dataFolder + saveRCritName,'wb') as f:
        pkl.dump(rCritData,f)