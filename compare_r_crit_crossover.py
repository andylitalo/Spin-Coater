# -*- coding: utf-8 -*-
"""
Created on Sat Feb 04 01:03:35 2017

@author: Andy
"""

import cPickle as pkl
import matplotlib.pyplot as plt
import numpy as np
import Functions as Fun

# Data parameters
dataFolder = 'Data\\'
rCritFile = 'r_crit.pkl'
rCrossoverFile = 'r_diffusion.pkl'
# Plot parameters
showRCrit = False
showRCrossover = False
# Save parameters
saveDiff = True
saveFolder = 'compare_r_crit_crossover\\'
saveName = 'crit_minus_diffusion'

# Derived Quantities and Constants
diskRad = 15.0 # cm
radiusNoisy = diskRad / 4


###############################################################################


# Open r_crit.pkl
# Load data
with open(dataFolder + rCritFile,'rb') as f:
    rCritData = pkl.load(f)

with open(dataFolder + rCrossoverFile, 'rb') as f:
    rCrossoverData = pkl.load(f)
        
# Extract parameter lists
QList = rCritData['QList']  
RPMList = rCritData['RPMList'] 
conditionList = rCritData['conditionList']

for condition in conditionList:
    rCritMat = rCritData[condition]
    rCrossoverMat = rCrossoverData[condition]
    
    if showRCrit:
        plt.contourf(RPMList, QList, rCritMat, cmap='hot', interpolation='nearest') # why this interpolation?
        plt.colorbar()
        # Formatting
        plt.xlabel('RPM')
        plt.ylabel('Q Target [mL/min]')
        pltTitle = 'Critical Radius, ' + condition
        plt.title(pltTitle)
        plt.show()
        plt.close()
    
    if showRCrossover:
        plt.pcolor(RPMList, QList, rCrossoverMat, cmap='hot') # why this interpolation?
        plt.colorbar()
        # Formatting
        plt.xlabel('RPM')
        plt.ylabel('Q Target [mL/min]')
        pltTitle = 'Crossover Radius, ' + condition
        plt.title(pltTitle)
        plt.show()
        plt.close()
       
    # Clean out zeros
    zeros = (rCritMat == np.zeros_like(rCritMat))
    rCrossoverMat[zeros] = 0
    diffMat = rCritMat - rCrossoverMat
    maxDiff = np.amax(diffMat)
    diffMat[zeros] = 15 
    # Difference Plot
    plt.pcolor(RPMList, QList, diffMat, cmap='hot', vmax=maxDiff+1) # why this interpolation?
    plt.colorbar()
    # Formatting
    plt.xlabel('RPM')
    plt.ylabel('Q Target [mL/min]')
    pltTitle = 'Difference in critical and crossover radii, ' + condition
    plt.title(pltTitle)

    if saveDiff: # this must come before plt.show() for proper saving
        savePlotName = Fun.get_save_name(saveFolder + saveName, 0, 0, condition, "")
        plt.savefig(savePlotName, bbox_inches='tight')
    plt.show()  
    plt.close()

            