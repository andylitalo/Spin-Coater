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
allDataFile = 'all_data.pkl'
# Save parameters
saveFolder = 'Data\\'
savePrefix = 'trials'
# Condition parameters
QList = np.array([500,1000,1500,2000,2500,3000,3250])
RPMList =  np.array([0,30,120,250,500,750,1000])
conditionList = ['Dry','Water', '3mM SDS', '6.5mM SDS', '10mM SDS', 'Water with 10mM SDS'] 


###############################################################################

 # Load data if not already loaded
if 'allData' not in globals():
    with open(dataFolder + allDataFile,'rb') as f:
        allData = pkl.load(f)

for condition in conditionList:
    completedTrials = np.zeros([len(QList), len(RPMList)])
    for r in range(len(QList)):
        Q = QList[r]
        for c in range(len(RPMList)):
            RPM = RPMList[c]
            key = (condition, Q, RPM)
            if key in allData.keys():
               completedTrials[r, c] = 1
                
    # Save estimated film thickness in an excel spreadsheet
    ## SAVE .CSV FILE ##
    # Create save name
    saveName = Fun.get_save_name(saveFolder + savePrefix, 0, 0, condition, '.csv')
    
    # Save file
    with open (saveName, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        completedTrialsRPM = np.vstack((RPMList, completedTrials))
        colQList = np.hstack((np.zeros(1), QList))
        completedTrialsRPMQ = np.vstack((colQList, np.transpose(completedTrialsRPM)))
        writer.writerows(np.transpose(completedTrialsRPMQ))