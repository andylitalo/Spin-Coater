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
import Functions as Fun

# Data parameters
dataFolder = 'Data\\'
dataFile = 'all_data.pkl'
# Save Parameters
savePlots = False
saveData = False
saveName = 'r_splash.pkl'
# Analysis parameters
timeStep = 0.01 # [s], 0.02 s is same time overwhich data is smoothed in locate_r_crit
diskRad = 15.0
radiusNoisy = diskRad / 4.0 # [cm] radius within which image-processing is poor
radiusMax = 10 # [cm] radius outside which stagnation will mean water reached wafer edge
# Condition parameters
QList = np.array([500,1000,1500,2000,2500,3000,3250])
RPMList =  np.array([0])
conditionList = ['Dry'] #['Dry','Water', '3mM SDS', '6.5mM SDS', '10mM SDS', 'Water with 10mM SDS'] 
# Plot parameters
showPlots = True
legend_fs = 12 # font size of legend

##############################################################################


# Cleanup existing windows
plt.close('all')
    
# Get list of data files to process
fileList = glob.glob(dataFolder + dataFile)
N1 = len(fileList)

# Initialize rCrit data structure
rSplashData = {}

# Analyze data
for i in range(N1):
    
    # Load data if not already loaded
    if 'allData' not in globals():
        with open(fileList[i],'rb') as f:
            allData = pkl.load(f)
    # Get data keys
    keyList = allData.keys()
    N = len(keyList)

    # Initialize data storage arrays
    for condition in conditionList:
        rSplashData[condition] = np.zeros([len(QList)])
    
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
        
        # Set the indices of the array for assignment of plotting data
        i_Q = np.argmax(QList==QTarget)

        # Parse remaining experiment data
        fps = data['fps']
        t0 = data['t0']
        time = data['time']
        aMax = data['aMax']
        aMean = data['aMean']
        
        # Skip if empty experiment
        if len(time) == 0:
            continue
        
        # Adjust data fields
        aMax = aMax[time > t0]
        aMean = aMean[time > t0]
        Q = Fun.convert_flowrate(QTarget)
        
        # Calculate number of frames to step when checking for stagnation
        step = int(round(fps*timeStep))
        
        # Record mean radius when aMax stagnates - this is what is helpful for identifying peaks in EP
        aDiff = aMax[step:len(aMax)] - aMax[0:len(aMax)-step]        
        rSplashInd = [i for i in range(len(aDiff)) if aDiff[i] == 0]
        if not rSplashInd:
            rSplash = radiusNoisy
            rSplashInd = 0 # TODO remove
        else:
            rSplashInd = rSplashInd[0]
            rSplash = aMean[rSplashInd]
            if rSplash > radiusMax:
                rSplash = radiusNoisy
         
        tRCrit = [time[k] for k in range(len(aMean)-1) if aMean[k+1] >= rSplash and aMean[k] < rSplash][0]
        print "Splash Radius is " + str(rSplash)
        print "Step is " + str(step) + " frames"
        print "aMax at splash index is " + str(aMax[rSplashInd]) # TODO remove
        print "aMax after step is " + str(aMax[rSplashInd + step]) # TODO remove
        rSplashData[condition][i_Q] = rSplash
       
        # Plot
        if showPlots:
            plt.figure()
            plt.plot(time[0:len(time)-step], 100*abs(aDiff), 'g.', linewidth=3, label='100aDiff')
            plt.plot(time[0:len(time)], aMax, 'k.', linewidth=3, label='aMax')
#            # Show critical point
            axes = plt.gca()
            yLim = axes.get_ylim()
            plt.plot((tRCrit, tRCrit), yLim, 'r--', label='critical point')
            plt.plot((tRCrit+timeStep, tRCrit+timeStep), yLim, 'r--', label='cpt end')
            # Formatting
            plt.xlabel('Time [sec]')
            plt.ylabel('aMax [cm]')
            plt.grid()
            plt.legend(loc=1, fontsize=legend_fs)
            pltTitle = 'Finding Splash Radius - Q = %i, Condition - %s' %(int(round(Q,-1)), condition)
            plt.title(pltTitle)
            plt.ylim([-1, diskRad])
#            if savePlots:
#                savePlotName = Fun.get_save_name(saveFolder, RPM, QTarget, condition, '')
#                plt.savefig(savePlotName, bbox_inches='tight')
            plt.show()
            plt.close()
 
# Saving rCrit
if saveData:
    rSplashData['QList'] = QList
    rSplashData['conditionList'] = conditionList       
    with open(dataFolder + saveName,'wb') as f:
        pkl.dump(rSplashData,f)