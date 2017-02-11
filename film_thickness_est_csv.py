# -*- coding: utf-8 -*-
"""
Created on Fri Feb 03 14:17:50 2017

@author: Andy
"""

import cPickle as pkl
import numpy as np
import Functions as Fun
import csv

# Data paramters
dataFolder = 'Data\\'
rCritFile = 'r_crit.pkl'
rinseDataFile = 'all_data.pkl'
# Save parameters
saveFolder = 'Data\\'

# Derived Quantities and Constants
diskRad = 15 # cm
radiusNoisy = diskRad / 4.0
m3PermL = 1E-6
minPerSec = 1.0/60.0
mPercm = 1.0/100.0


###############################################################################

 # Load data if not already loaded
if 'rinseData' not in globals():
    with open(dataFolder + rinseDataFile,'rb') as f:
        rinseData = pkl.load(f)
    
# Open r_crit.pkl
# Load data
if 'rCritData' not in globals():
    with open(dataFolder + rCritFile,'rb') as f:
        rCritData = pkl.load(f)

# Extract parameter lists
QList = rCritData['QList']  
RPMList = rCritData['RPMList'] 
conditionList = rCritData['conditionList']

for condition in conditionList:
    rCritMat = rCritData[condition] # matrix of r_crit where rows=Q and cols=RPM
    filmThicknessMat = np.zeros_like(rCritMat)
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
            if key not in rinseData.keys():
               filmThicknessMat[r, c] = 0
            else:
                data = rinseData[key]
    
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
                filmThicknessEst_m = volDispensed_m3 / wettedArea_m2 # m
                filmThicknessMat[r, c] = filmThicknessEst_m / mPercm # cm
        
                print (condition, RPM, QTarget, rCrit)
                print (volDispensed_m3 / m3PermL, wettedArea_m2 / (mPercm)**2, filmThicknessEst_m / mPercm)
                
    # Save estimated film thickness in an excel spreadsheet
    ## SAVE .CSV FILE ##
    # Create save name
    saveName = saveFolder + "cond_" # rounds Q to nearest 10s
    if condition.find('SDS') < 0:
        saveName += condition
    elif condition.find('Water') < 0:
        saveName += 'SDS_'
        concentration = float(condition[0:condition.find('mM')])
        saveName += "%d" % int(10*concentration)
    else:
        saveName += 'SDS_100_on_Water'
    saveName += '.csv'
    
    # Save file
    with open (saveName, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        filmThicknessRPM = np.vstack((RPMList, filmThicknessMat))
        colQList = np.hstack((np.zeros(1), QList))
        filmThicknessRPMQ = np.vstack((colQList, np.transpose(filmThicknessRPM)))
        writer.writerows(np.transpose(filmThicknessRPMQ))