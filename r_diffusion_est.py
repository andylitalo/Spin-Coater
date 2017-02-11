# -*- coding: utf-8 -*-
"""
Created on Fri Feb 03 14:17:50 2017

@author: Andy
"""

import cPickle as pkl
import numpy as np
import Functions as Fun

# Data paramters
dataFolder = 'Data\\'
rCritFile = 'r_crit.pkl'
allDataFile = 'all_data.pkl'
# Save parameters
saveData = True
saveFolder = 'Data\\'
saveName = 'r_diffusion.pkl'

# Derived Quantities and Constants
diskRad = 15.0 # cm
radiusNoisy = diskRad / 4
rho_kgm3 = 1000
mu_kgms = 1E-3
m3PermL = 1E-6
minPerSec = 1.0/60.0
radPerRot = 2*np.pi
mPercm = 1.0/100.0


###############################################################################


 # Load data if not already loaded
if 'allData' not in globals():
    with open(dataFolder + allDataFile,'rb') as f:
        allData = pkl.load(f)
    
# Open r_crit.pkl
# Load data
with open(dataFolder + rCritFile,'rb') as f:
    rCritData = pkl.load(f)

# Extract parameter lists
QList = rCritData['QList']  
RPMList = rCritData['RPMList'] 
conditionList = rCritData['conditionList']

rDiffusion = {}
for condition in conditionList:
    rCritMat = rCritData[condition] # matrix of r_crit where rows=Q and cols=RPM
    rDiffusionMat = np.zeros_like(rCritMat)
    for r in range(len(QList)):
        QTarget = QList[r]
        Q = Fun.convert_flowrate(QTarget)
        for c in range(len(RPMList)):
            RPM = RPMList[c]
            # Change zeros to heuristic for radius with low error from rivulets/image-processing noise
            rCrit = rCritMat[r, c]
            if rCrit == 0:
                rCrit = radiusNoisy + 1 # 1 cm added to be far from noise radius
            
            key = (condition, QTarget, RPM)
            if key not in allData.keys():
               rDiffusionMat[r, c] = 0
            elif RPM == 0:
                rDiffusionMat[r, c] = diskRad
            else:
                data = allData[key]
    
                # Parse remaining experiment data
                fps = data['fps']
                t0 = data['t0']
                time = data['time']
                aMax = data['aMax']
                aMin2 = data['aMin2']
                aMean = data['aMean']
                eP = data['excessPerimeter']
                
                indCrit = [i for i in range(len(aMean)-1) if (aMean[i+1] >= rCrit and aMean[i] <= rCrit)][0]
                tCrit = time[indCrit]
                tFlow = tCrit - t0
                Q_m3s = m3PermL * minPerSec * Q # m^3/s
                volDispensed_m3 = Q_m3s * tFlow # m^3
                wettedArea_m2 = np.pi * (mPercm * aMean[indCrit])**2 # m^2
                d_m = volDispensed_m3 / wettedArea_m2 # m
                tD = rho_kgm3*d_m**2/mu_kgms # s
                if tD < time[0] - t0:
                    indDiffusion = 0
                else:
                    indDiffusion = [i for i in range(len(time)-1) if (time[i+1]-t0 >= tD and time[i]-t0 <= tD)][0]
                rDiffusionMat[r, c] = aMean[indDiffusion] # cm
          
    print rDiffusionMat
    # Save estimated film thickness in an excel spreadsheet
    rDiffusion[condition] = rDiffusionMat
    if saveData:
        with open(saveFolder + saveName, 'wb') as f:
            pkl.dump(rDiffusion, f)