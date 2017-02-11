# -*- coding: utf-8 -*-
"""
Created on Sun Jan 29 07:57:09 2017

@author: Andy
"""

import matplotlib.pyplot as plt
import cPickle as pkl
import numpy as np
import Functions as Fun

# Data parameters
dataFolder = 'Data\\'
dataFile = 'r_crit_minSepFrac_0-2.pkl'
# Save parameters
savePlot = True
saveFolder = 'plot_r_crit\\'
saveName = 'r_crit_vs_RPM_minSepFrac_0-2_log'
# Disk Parameters
diskRad = 15.0 # [cm]
radiusNoisy = diskRad / 4.0
# Plot parameters
logScale = True
legend_fs = 8
xLim = [100, 1000]
yLim = [diskRad/10, diskRad]


###############################################################################


# Load data
with open(dataFolder + dataFile,'rb') as f:
    rCritData = pkl.load(f)

# Extract parameter lists
QList = rCritData['QList']  
RPMList = rCritData['RPMList'] 
conditionList = rCritData['conditionList']

# Plot on one figure
fig = plt.figure()
for i in range(len(conditionList)):
    condition = conditionList[i]
    rCritMat = rCritData[condition] # matrix of r_crit where rows=Q and cols=RPM
    # Plot one trendline per Q
    for r in range(len(QList)):
        QTarget = QList[r]
        Q = Fun.convert_flowrate(QTarget)
        rCrit = rCritMat[r, :]
        if np.size(np.nonzero(rCrit)) == 0:
            continue # skip flow rates with all rCrit = 0
            print Q
    
        # Clean out zeros
        cleanInds = [j for j in range(len(RPMList)) if rCrit[j] != 0]
        RPMListClean = RPMList[cleanInds]
        rCritClean = rCrit[cleanInds]
        # Plot trendline
        label = 'Q = %i mL/min' % (int(round(Q, -1)))
        plt.plot(RPMListClean, rCritClean, linewidth=2, marker='o', label=label)
    
    # Format plot
    plt.xlabel('Rotation Rate [RPM]')
    plt.ylabel('Critical Radius [cm]')
    if logScale:
        plt.xscale('log')
        plt.yscale('log')
    plt.xlim(xLim)
    plt.ylim(yLim)
    plt.grid()
    plt.legend(loc=0, fontsize=legend_fs)
    pltTitle = 'Critical Radius vs. RPM, %s' %(condition)
    plt.title(pltTitle)
    if savePlot:
        condition = condition.replace('.','-') # can't save with decimals in name
        savePlotName = saveFolder + saveName + '_cond_' + condition
        plt.savefig(savePlotName, bbox_inches='tight')
    plt.show()
    plt.close()
