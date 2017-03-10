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
rCritFilePath = 'Data\\r_crit.pkl'
numInfoFields = 7 # number of fields in data dictionary for info on experiment--7 default, 8 if stride included
selectedK = 15
frame2Sample = 800
diskRad = 15
diskFrac = 0.9


# Load data
if 'allData' not in globals():
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

# Load rCrit Data
with open(rCritFilePath, 'rb') as f:
    rCritData = pkl.load(f)

QList = rCritData['QList']
RPMList = rCritData['RPMList']
i_Q = np.argmax(QList==flowRate)
j_RPM = np.argmax(RPMList==RPM)
rCrit = rCritData[condition][i_Q, j_RPM]

NMax = 0
nFrames = len(keys)-numInfoFields
iStart = 0
iEnd = 0
aMeanPrev = 0
for i in range(nFrames):
    contour = data[keys[i]]
    # Update max number of Fourier modes, NMax
    s = contour[0,:]
    N = len(s)
    NMax = max(NMax, N)
    # Check if unconstrained rivulet growth stage reached (iStart:iEnd)
    a = contour[1,:]
    aMean = np.mean(a)
    aMax = np.max(a)
    if aMeanPrev < rCrit and aMean >= rCrit:
        iStart = i
    aMeanPrev = aMean
    if aMax >= diskFrac*diskRad:
        iEnd = i
        break
 
NMax = NMax - (NMax%2) # round NMax down to nearest even number
kMat = np.zeros([NMax, iEnd - iStart])
rZ = NMax/2 # row for zero mode (k = 0)

for i in range(iStart, iEnd):    
    contour = data[keys[i]]

    s = contour[0,:]
    a = contour[1,:]
    N = len(s)
    
    # fit cubic spline
    cs = scipy.interpolate.interp1d(s, a, kind='cubic')
    
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
    kMat[rZ:(rZ+N/2),i] = k[:N/2] # k organized from n = 0:N/2-1, -N/2:-1 modes
    kMat[(rZ-N/2):rZ,i] = k[(N+1)/2:] # +1 skips largest mode for cases with odd N
    
plt.figure()
plt.plot(sSample, aSample)
plt.title('Interpolated s vs. a - Smoothed - Sampled from frame #' + str(frame2Sample))
plt.xlabel('s')
plt.ylabel('a')