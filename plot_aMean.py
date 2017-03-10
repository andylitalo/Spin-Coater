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
import Functions as Fun\

# Data parameters
dataFolder = 'Data\\'
dataFile = 'all_data.pkl'
rCritFile = 'r_crit.pkl'
# Save Parameters
savePlots = True
saveData = False
saveFolder = 'plot_aMean\\'
# Analysis parameters
diskFrac = 0.9 # fraction of disk considered in analysis
diskRad = 15.0
#QList = np.array([500,1000,1500,2000,2500,3000,3250])
#RPMList =  np.array([0, 30, 120, 250, 500, 750, 1000])
#conditionList = ['Dry','Water', '3mM SDS', '6.5mM SDS', '10mM SDS', 'Water with 10mM SDS'] 
# Plot parameters
showPlots = True
legend_fs = 9 # font size of legend

##############################################################################


# Loop through datasets

# Cleanup existing windows
plt.close('all')
    
# Get list of data files to process
fileList = glob.glob(dataFolder + dataFile)
N1 = len(fileList)

# Initialize rTrans data structure
rTransData = {}

# Analyze data
for i in range(N1):
    
    # Load data if not already loaded
    if 'allData' not in globals():
        with open(fileList[i],'rb') as f:
            allData = pkl.load(f)
    # Get data keys
    keyList = allData.keys()
    N = len(keyList)
    
    with open(dataFolder + rCritFile, 'rb') as f:
        rCritData = pkl.load(f)
        
    # Extract parameter lists
    QList = rCritData['QList']  
    RPMList = rCritData['RPMList'] 
    conditionList = rCritData['conditionList']

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
        saveKey = (RPM, QTarget, condition)
        print saveKey
        

        # Load aMean, time data

        # Parse remaining experiment data
        fps = data['fps']
        t0 = data['t0']
        time = data['time']
        aMean = data['aMean']        
        
        # Skip if empty experiment
        if len(time) == 0:
            continue
        
        # Adjust data fields
        aMean = aMean[time > t0]
        Q = Fun.convert_flowrate(QTarget)
        
              
        # Load rCrit
        rCritMat = rCritData[condition] # matrix of r_crit where rows=Q and cols=RPM
        i_Q = np.argwhere(QList == QTarget)
        j_RPM = np.argwhere(RPMList == RPM)
        rCrit = rCritMat[i_Q, j_RPM][0]
        
        # Plot aMean data on semilog plot
        
         # Plot
        if showPlots:
            plt.figure()
#            plt.plot(time[0:len(time)-step], 100*abs(aDiff), 'g.', linewidth=3, label='100aDiff')
            plt.semilogy(time[0:len(aMean)], aMean, 'k.', linewidth=3, label='aMean')
#            # Show critical point
            axes = plt.gca()
            xLim = axes.get_xlim()
            plt.plot(xLim, (rCrit, rCrit), 'r--', label='critical point')
#            plt.plot((tRCrit+timeStep, tRCrit+timeStep), yLim, 'r--', label='cpt end')
            # Formatting
            plt.xlabel('Time [sec]')
            plt.ylabel('aMean [cm]')
            plt.grid()
            plt.legend(loc=0, fontsize=legend_fs)
            pltTitle = 'Mean Radius: %i RPM, Q = %i, Condition - %s' %(RPM, int(round(Q,-1)), condition)
            plt.title(pltTitle)
            if savePlots:
                savePlotName = Fun.get_save_name(saveFolder, RPM, QTarget, condition, '')
                plt.savefig(savePlotName, bbox_inches='tight')
            plt.show()
            plt.close()
 

# Locate sharp differences in slope

## Save radii of sharp differences
#       
#       
## Saving rCrit
#if saveData:
#    rSplashData['QList'] = QList
#    rSplashData['conditionList'] = conditionList       
#    with open(dataFolder + saveName,'wb') as f:
#        pkl.dump(rSplashData,f)