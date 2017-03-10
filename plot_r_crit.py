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
dataFile = 'r_crit.pkl'
# Save parameters
savePlot = True
saveFolder = 'plot_r_crit\\'
saveName = 'r_crit_log_aristoff'
# Disk Parameters
diskRad = 15.0 # [cm]
radiusNoisy = diskRad / 4.0
aJet = 0.2 # jet radius [cm]
nu = 0.01 # kinematic viscosity of water in cgs [Poise/(g/cm3)]
# Plot parameters
logScale = True
aristoffCond = True
legend_fs = 8
xLim = [100, 1000]
yLim = [1, diskRad]
# Conversions
s_per_min = 1./60


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
cmap = plt.get_cmap('jet')
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
        color = cmap(float(r)/len(QList))
        label = 'Q = %i mL/min' % (int(round(Q, -1)))
        plt.plot(RPMListClean, rCritClean, linewidth=2, marker='o', label=label,c=color)
        if aristoffCond:
            Qcm3s = Q*s_per_min
            Re = Qcm3s/(aJet*nu)
            rAristoff = 0.315*aJet*Re**(1./3)
            plt.plot(RPMListClean, rAristoff*np.ones_like(RPMListClean), linewidth=2,c=color)
    
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
