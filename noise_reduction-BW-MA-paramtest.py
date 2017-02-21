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
import csv
import numpy as np
from scipy.signal import butter, lfilter, freqz
import matplotlib.pyplot as plt
import copy
import Functions as Fun
import SignalProcessingFunctions as SPF

# User parameters
folder = 'Code from John 10-14\\'
# saveFolder = '..\\..\\..\\Data\\t_aMax_aMin2_aMean_EP\\'
name = 'all_data.pkl'
QList = np.array([500]) #np.array([500,1000,1500,2000,2500,3000,3250])
RPMList = np.array([120]) #np.array([0,30,120,250,500,750,1000])
conditionList = ['Dry'] #['Dry','Water', '3mM SDS', '6.5mM SDS', '10mM SDS', 'Water with 10mM SDS'] 
diskFrac = 0.9
viewTimeSignal = False

# Derived Quantities and Constants
diskRad = 15 # cm
wetTarget = diskFrac * diskRad

# Cleanup existing windows
plt.close('all')
    
# Get list of data files to process
fileList = glob.glob(folder+name)
N1 = len(fileList)

for j in range(N1):
    
    # Load data and get list of fields
    if 'allData' not in globals():
        with open(fileList[j],'rb') as f:
            allData = pkl.load(f)
            print fileList[j]
    keyList = allData.keys()
    N = len(keyList)
    
    # Initialize data storage arrays
    plotData = {}
    for condition in conditionList:
        plotData[condition] = {}
        plotData[condition]['ePMax'] = np.zeros((len(QList),len(RPMList)))
        plotData[condition]['wetTime'] = np.zeros((len(QList),len(RPMList)))
        plotData[condition]['wetVolume'] = np.zeros((len(QList),len(RPMList)))

    
    for i in range(N):
    
        # Parse the experiment data
        data = allData[keyList[i]]
        RPM = data['RPM']
        QTarget = data['flowRate']
        condition = data['condition']
        if condition not in conditionList:
            continue
        if (QTarget not in QList) or (RPM not in RPMList):
            continue
        # Set the indicies of the array for assignment of plotting data
        ind0 = np.argmax(QList==QTarget)
        ind1 = np.argmax(RPMList==RPM)

        # Parse remaining experiment data
        fps = data['fps']
        t0 = data['t0']
        time = data['time']
        aMax = data['aMax']
        aMin2 = data['aMin2']
        aMean = data['aMean']
        eP = data['excessPerimeter']
        
        if len(time) < 1:
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
        
        
        ################### BUTTERWORTH FILTER ################################ 
        print 'Butterworth Filter'
        for cutoff in [10,20,30]:
            # Filter requirements.
            order = 6
            fs = fps       # sample rate, Hz
    #        cutoff = 20  # desired cutoff frequency of the filter, Hz
            
            # Get the filter coefficients so we can check its frequency response.
            b, a = SPF.butter_lowpass(cutoff, fs, order)
            
            # Plot the frequency response.
            w, h = freqz(b, a, worN=8000)
            plt.subplot(2, 1, 1)
            plt.plot(0.5*fs*w/np.pi, np.abs(h), 'b')
            plt.plot(cutoff, 0.5*np.sqrt(2), 'ko')
            plt.axvline(cutoff, color='k')
            plt.xlim(0, 0.5*fs)
            plt.title("Lowpass Filter Frequency Response")
            plt.xlabel('Frequency [Hz]')
            plt.grid()
    
            # Filter the data, and plot both the original and filtered signals.
            ePFiltered = SPF.butter_lowpass_filter(eP, cutoff, fs, order)
            
            plt.subplot(2, 1, 2)
            plt.plot(time, eP, 'b-', label='data')
            plt.plot(time, ePFiltered, 'g-', linewidth=2, label='filtered data')
            plt.xlabel('Time [sec]')
            plt.grid()
            plt.legend()
            
            plt.subplots_adjust(hspace=0.35)
            plt.show()
            
        #######################################################################
        ###################### LOCAL AVERAGE #################################
        print 'Local Average Filter'
        # Filter
        for n in range(0,101,20):
            eP_LA = SPF.local_average_filter(eP, n)
            
            # Plot
#            plt.plot(time, eP, 'b-', label='data')
            plt.plot(time, eP_LA, 'g-', linewidth=2, label='filtered data')
            plt.xlabel('Time [sec]')
            plt.grid()
            plt.legend()
            pltTitle = 'Local Average Filtering with n = ' + str(n)
            plt.title(pltTitle)
            plt.show()
            
        #######################################################################
        ######################### MEDIAN FILTER ###############################
        print 'Median Filter'
        # Filter
        for n in range (1, 156, 20):
            eP_MF = signal.medfilt(eP, n) # note: n must be odd
            
            # Plot
#            plt.plot(time, eP, 'b-', label='data')
            plt.plot(time, eP_MF, 'g-', linewidth=2, label='filtered data')
            plt.xlabel('Time [sec]')
            plt.grid()
            plt.legend()
            pltTitle = 'Median Filtering with n = ' + str(n)
            plt.title(pltTitle)
            plt.show()
            
        #######################################################################
        ######################## KILL SPIKES (cf. Tom O'Haver, UMD, 2016) #####
        ### !!! TODO seems to do nothing...
        print 'Kill Spikes'
        # Parameters
        thresh = 0.2
        width = 150
        # Filter
        eP_KS = SPF.kill_spikes_filter(time, eP, width, thresh)
        
        # Plot
#            plt.plot(time, eP, 'b-', label='data')
        plt.plot(time, eP_KS, 'g-', linewidth=2, label='filtered data')
        plt.xlabel('Time [sec]')
        plt.grid()
        plt.legend()
        pltTitle = 'Kill-spikes filtering with width = ' + str(width) + ' and thresh = ' + str(thresh)
        plt.title(pltTitle)
        plt.show()
        
        #######################################################################
        ##################      
#            