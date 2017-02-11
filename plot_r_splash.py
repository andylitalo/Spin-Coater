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
dataFile = 'r_splash.pkl'
# Save parameters
savePlot = False
saveFolder = 'plot_r_crit\\'
saveName = 'r_splash_vs_Q'
# Disk Parameters
diskRad = 15.0 # [cm]
radiusNoisy = diskRad / 4.0
# Plot parameters
logScale = False
legend_fs = 8
xLim = [300, 3000]
yLim = [diskRad/100, diskRad/2]


###############################################################################


# Load data
with open(dataFolder + dataFile,'rb') as f:
    rSplashData = pkl.load(f)

# Extract parameter lists
QList = rSplashData['QList']  
conditionList = rSplashData['conditionList']

# Plot on one figure
fig = plt.figure()
for i in range(len(conditionList)):
    condition = conditionList[i]
    rSplash = rSplashData[condition] # matrix of r_crit where rows=Q and cols=RPM
    QActual = []
    # Plot one trendline per Q
    for r in range(len(QList)):
        QActual += [Fun.convert_flowrate(QList[r])]
        
    # Plot trendline
    label = condition
    plt.plot(QActual, rSplash, linewidth=2, marker='o', label=label)
    
    # Format plot
    plt.xlabel('Jet Flow Rate [RPM]')
    plt.ylabel('Splash Radius [cm]')
    if logScale:
        plt.xscale('log')
        plt.yscale('log')
    plt.xlim(xLim)
    plt.ylim(yLim)
    plt.grid()
    plt.legend(loc=0, fontsize=legend_fs)
    pltTitle = 'Splash Radius vs. Q, %s' %(condition)
    plt.title(pltTitle)
    if savePlot:
        condition = condition.replace('.','-') # can't save with decimals in name
        savePlotName = saveFolder + saveName + '_cond_' + condition
        plt.savefig(savePlotName, bbox_inches='tight')
    plt.show()
    plt.close()
