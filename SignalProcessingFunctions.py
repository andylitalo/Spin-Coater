# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 10:49:32 2017

@author: Andy
"""

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

def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y
    
def local_average_filter(data, n):
    '''
    Returns locally averaged data where average is taken over n nearest-
    neighbors.
    '''
    data_LA = copy.copy(data)
    nPts = len(data)
    for i in range(1,n+1):
        data_LA += np.concatenate((data[0:i], data[0:nPts-i]))
        data_LA += np.concatenate((data[i:], data[nPts-i:]))
    data_LA /= 2*n + 1
    return data_LA
    
def kill_spikes_filter(x_data, y_data, width, thresh):
    data_KS = copy.copy(y_data)
    for i in range(len(data_KS)-width):
        if abs(data_KS[i+1] - data_KS[i]) > thresh:
            for j in range(width):
                data_KS[i+j] = np.interp(x_data[i+j], [x_data[i],x_data[i+width]], [data_KS[i], data_KS[i+width]])
    return data_KS