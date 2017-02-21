# -*- coding: utf-8 -*-
"""
Created on Sat Feb 18 20:00:18 2017

@author: Andy
"""


import scipy
import cPickle as pkl
import scipy.interpolate
import numpy as np
from numpy.fft import fft
import matplotlib.pyplot as plt
import SignalProcessingFunctions as SPF
from scipy import signal


# User parameters
folder = 'Data\\Fourier Decompositions\\'
fileName = 'Dry_0120RPM_0500mLmin_0500FPS_contour.pkl'
numInfoFields = 7 # number of fields in data dictionary for info on experiment
selectedK = 15 # fourier component to plot--for some reason >=20 isn't included in kVals?
frame2Sample = 800


# Load data
with open(folder + fileName, 'rb') as f:
    allData = pkl.load(f)
   
data = allData['Dry_0120RPM_0500mLmin_0500FPS_c']
condition = data['condition']
flowRate = data['flowRate']
RPM = data['RPM']
fps = data['fps']
t0 = data['t0']
time = data['time']
tMaxEP = data['tMaxEP']
keys = data.keys()

kMag = np.zeros([len(time),1])
sSample = []
aSample = []
for i in range(len(keys)-numInfoFields):
    
    contour = data[keys[i]]

    s = contour[0,:]
    a = contour[1,:]
    N = len(s)
    
    # fit cubic spline
    cs = scipy.interpolate.interp1d(s, a) # find better interpolation function
    
    # resample
    x = np.linspace(0, max(s), N)
    y = cs(x)
    
    # Smooth s vs. a
    # Filter parameters
    unscaledRPM = 250 
    timeWindow = 0.02 # sec
    nn_MF = int(timeWindow*fps*unscaledRPM/RPM) # number of nearest neighbors in filtering window
    nn_LA = int(timeWindow*fps*unscaledRPM/RPM) # number of nearest neighbors in filtering window
    
    # Filters
    y_MF = signal.medfilt(y, 2*nn_MF + 1) # median filter
    y_LA = SPF.local_average_filter(y_MF, nn_LA)
    
    yFilt = y_LA
    
    if keys[i] == frame2Sample:
        sSample = x
        aSample = yFilt
        
    # Frame-by-frame Fourier Analysis - record dominant mode
    k = fft(yFilt)
    kVals = np.concatenate((np.arange(0,N/2-1),np.arange(-N/2,-1)))
    kMag[i] = k[(kVals == selectedK)]
    
# plot magnitude of a few dominant modes over time, as in Fraysse and Homsy
plt.plot(time[time<=tMaxEP], abs(kMag[time<=tMaxEP]))
plt.xlabel('Time (s)')
plt.ylabel('Magnitude of 3rd Fourier Component')

plt.figure()
plt.plot(x, y)
plt.title('Interpolated s vs. a')
plt.xlabel('s')
plt.ylabel('a')

plt.figure()
plt.plot(x, yFilt)
plt.title('Interpolated s vs. a - Smoothed')
plt.xlabel('s')
plt.ylabel('a')

plt.figure()
plt.plot(sSample, aSample)
plt.title('Interpolated s vs. a - Smoothed - Sampled from frame #' + str(frame2Sample))
plt.xlabel('s')
plt.ylabel('a')