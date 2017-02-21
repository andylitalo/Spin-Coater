# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 15:34:40 2015

@author: Fuller
"""

import numpy as np
import cv2

class Stats:

    def __init__(self, name=None, fullName=None, ylabel=None, style='b.', defaultVal=0, color=(0,0,1)):
        self.name = name
        self.fullName = fullName
        self.ylabel = ylabel
        self.style = style
        self.defaultVal = defaultVal
        
    def set_name(self, name):
        self.name = name
        
    def get_name(self):
        return self.name

    def set_full_name(self, fullName):
        self.fullName = fullName
        
    def get_full_name(self):
        return self.fullName

    def set_ylabel(self, ylabel):
        self.ylabel = ylabel
        
    def get_ylabel(self):
        return self.ylabel 
        
    def set_style(self, style):
        self.style = style
        
    def get_style(self):
        return self.style

    def set_defaultVal(self, defaultVal):
        self.defaultVal = defaultVal
        
    def get_defaultVal(self):
        return self.defaultVal
        
    def set_color(self, color):
        self.color = color
        
    def get_color(self):
        return self.color
        
    def generate_circle(self,R,center,N=100):
        """
        Generate an array of x and y values that lie evenly spaced on a circle 
        with the specified center and radius.
        """
        theta = np.linspace(0.0,2.0*np.pi,N)
        y = R*np.sin(theta) + center[1]
        x = R*np.cos(theta) + center[0]
        return x,y
        
    def generate_periodic_circle(self,R, var, xCenter, yCenter, freq, fn):
        """
        Generates points that map out a circle with a given periodic function
        overlayed.
        
        fn can be the following:
        Fun.sine -> overlays sine wave with amplitude "var"
        Fun.sawtooth -> overlays the increasing sawtooth that begins at (R-var) and
        increases up to (R+var) in each period
        Fun.triangle -> overlays the increasing-decreasing triangle, which begins at
        (R-var), reaches (R+var) halfway through the period, and then decreases
        back to (R-var) by the end of the period.
        """
        # The perimeter for a triangle wave is about sqrt((circum)^2+(2*var*freq)^2),
        # A little bit is added to get the formula used below, just in case.
        nPts = int(2*np.sqrt((2*np.pi*R)**2.0 + (2*var*freq)**2.0)) # ensure sufficient points to fill perimeter
        theta = np.linspace(0,2*np.pi, nPts)
        # Initalize xy, array of tuples of x- and y-values
        xy = np.array([[0,0]])
        if freq == 0:
            period = 0
        else:
            period = 2*np.pi/freq
        for i in range(len(theta)):
            pert = fn(theta[i], period, var)
            r = R + pert
            # Convert from polar to Cartesian coordinates
            x = r*np.cos(theta[i]) + xCenter
            y = r*np.sin(theta[i]) + yCenter
            xy = np.concatenate((xy,np.array([[x,y]])))
        return xy[1:len(xy)]

    def sawtooth(self,x, period, var):
        """
        Calculates y-value for an increasing sawtooth with amplitude "var" above
        and below y = 0, i.e. it begins at (0,-var), passes through (period/2,0),
        and ends at (period,+var) in each period. 
        """
        if period == 0:
            return x
        else:
            x = x % period
            frac = x/period
            y = (2*frac-1)*var
            return y
        
    def sine(self,x, period, var):
        """
        Calculates y-value for a sine wave starting at the origin.
        """
        return var*np.sin((2*np.pi/period)*x)
        
    def triangle(self,x, period, var):
        """
        Calculates y-value for an increasing-decreasing triangle function with 
        amplitude "var" above and below y = 0, i.e. it begins at (0,-var), 
        increases to (period/2,+var), and decreases to (period,-var).
        """
        if period == 0:
            return x
        else:
            x = x % period
            frac = 2*abs(0.5 - (x/period))
            y = (1-2*frac)*var
            return y
            
    def get_corrected_arclength(self,pts,closed=False):
        """
        Smooths the digital curve defined by the row-column tuples in the numpy
        array "pts" using a 5-point average, i.e., it replaces each tuple with an 
        average of the two previous points, the two succeeding points, and the 
        point itself. The arc-length is then calculated by scaling up the image
        to 3-decimal-place precision, applying the OpenCV arcLength function, and
        scaling back down.
        """
        
        l = len(pts)
        ptsDown2 = np.concatenate((pts[2:l],pts[0:2]))
        ptsDown1 = np.concatenate((pts[1:l],np.array([(pts[0][0],pts[0][1])])))
        ptsUp1 = np.concatenate((np.array([(pts[l-1][0],pts[l-1][1])]),pts[0:l-1]))
        ptsUp2 = np.concatenate((pts[l-2:l],pts[0:l-2]))
        summedPts = ptsDown2 + ptsDown1 + pts + pts + ptsUp1 + ptsUp2
        avePts = summedPts/5.0
        zoomAvePts = 1000.0*avePts
        arcLength = cv2.arcLength(zoomAvePts.astype(int),closed)/1000.0
        
        return arcLength
        
    def get_value(self, statName, aveRad, diskRad, 
                    diskCenter,r,rMin, rMax, amp,adjustedArcLength, wettedArea):
                      
        if statName == "aveRad":
            nonDimRad = float(aveRad)/diskRad
            
            return nonDimRad

                
        if statName == "rMin":
            nonDimRMin = float(rMin)/diskRad
            
            return nonDimRMin

            
        if statName == "rMax":
            nonDimRMax = float(rMax)/diskRad
            
            return nonDimRMax

            
        #TODO Debug calculation of excPerim so diskMaskIm returns 1.0
        if statName == "excPerim":
            center = (2*aveRad,2*aveRad)
            X,Y = self.generate_circle(aveRad,center)
            xy = np.array([(X[l],Y[l]) for l in range(len(X))])
            perfectCirc = float(self.get_corrected_arclength(xy.astype(int),closed=True))
            excPerim = adjustedArcLength/perfectCirc #(2*np.pi*aveRad)/1.09#(2*np.pi*diskRad/diskCirc)            
            
            return excPerim

            
        if statName == "fracWetArea":
            return wettedArea/(np.pi*diskRad**2.0)

            
        if statName == "sineFreq":
            freq = 0
            while True:
                sinusoidalCnt = self.generate_periodic_circle(aveRad,amp,diskCenter,freq,self.sine)
                testArcLength = cv2.arcLength(sinusoidalCnt.astype(int),True)
                if testArcLength >= adjustedArcLength:
                    break
                else:
                    freq += 1

            return freq
            
        if statName == "sawtoothFreq":
            freq = 0
            while True:
                sawtoothCnt = self.generate_periodic_circle(aveRad,amp,diskCenter,freq,self.sawtooth)
                testArcLength = cv2.arcLength(sawtoothCnt.astype(int),True)
                if testArcLength >= adjustedArcLength:
                    break
                else:
                    freq += 1 

            return freq
            
        if statName == "triangleFreq":
            freq = 0
            while True:
                triangleCnt = self.generate_periodic_circle(aveRad,amp,diskCenter,freq,self.triangle)
                testArcLength = cv2.arcLength(triangleCnt.astype(int),True)
                if testArcLength >= adjustedArcLength:
                    break
                else:
                    freq += 1 

            return freq
        
        if statName == "quarterAmplFreq":

            # Calculate frequency                
            rQ = rMin + 1/float(4)*(rMax-rMin); freq = 0
            for l in range(1,len(r)):
                rQOutPrev = r[l-1] > rQ
                rQOut = r[l] > rQ
                if rQOutPrev and not rQOut:
                    freq += 1
                    
            return freq
            
        if statName == "halfAmplFreq":
            # Calculate frequency                
            rHalf = (rMax+rMin)/2; freq = 0
            for l in range(1,len(r)):
                rHalfOutPrev = r[l-1] > rHalf
                rHalfOut = r[l] > rHalf
                if rHalfOutPrev and not rHalfOut:
                    freq += 1
                    
            return freq
            
            #TODO Fix this calculation
        if statName == "fftFreq":
            fft = np.fft.rfft(r)
            fft = fft[1:]
            maxVal = fft.max()
            maxMode = fft.searchsorted(maxVal)+1 # +1 corrects for removal of 0th element
        
            return maxMode
            
        if statName == "stDev":
            stDev = np.std(r)
            scaledStDev = stDev/aveRad
            
            return scaledStDev
            
    def main(self, statsList, aveRad, diskRad, diskCenter,r,rMin, rMax, amp,
                  adjustedArcLength, wettedArea):

        statsVec = np.zeros(len(statsList))
        for i in range(len(statsList)):
            statsVec[i] = self.get_value(statsList[i], aveRad, diskRad, 
                    diskCenter,r,rMin, rMax, amp,adjustedArcLength, wettedArea)
                  
        return statsVec
        
if __name__ == '__main__':
    print