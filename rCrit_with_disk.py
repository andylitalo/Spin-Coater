## -*- coding: utf-8 -*-

import glob
import matplotlib.pyplot as plt
import cPickle as pkl
import numpy as np
from scipy import signal

import VideoFunctions as VF
import ImageProcessingFunctions as IPF
import Functions as Fun

# User parameters
folder = '..\\Data\\New Disk\\'
name = 'all_data_SDS.pkl'
rCritFile = 'r_crit_splash.pkl'
RPMList = [750] # [0,30,120,250,500,750,1000]
QTargetList = [500] # [500,1000,1500,2000,2500,3000,3250]
conditionTargetList = ['3mM SDS'] # ['Dry','Water','3mM SDS','6.5mM SDS','10mM SDS']
wetTarget = 0.9*15
diskScale = 15
step = 1
ePcutoff = 5
figscale = (10,5)
DPI = 150
removeMaskPoints = False
noRotate = True
shiftCenter = True
trueCenter = [509,547]


# Plot properties
FS1 = 18 # axes labels
LW = 2 # Line width for plots
   
# Cleanup existing windows
plt.close('all')
pi = np.pi
invert = True
    
# Load data
temp = glob.glob(folder+name)[0]
with open(temp,'rb') as f:
    allData = pkl.load(f)
    
with open(rCritFile, 'rb') as f:
    rCritData = pkl.load(f)
    
# Compute plotting properties
keyList = allData.keys()
N = len(keyList)
nData = len(RPMList)
dPhi = 2.0*np.pi/float(nData)
offset = -(pi/2 - dPhi)


for conditionTarget in conditionTargetList:
    
    for QTarget in QTargetList:
        # Initialize figure variables
        plt.figure(1,figsize=figscale)

        for i in range(0,N):

            # Parse the filename to get video info
            keyName = keyList[i]
            data = allData[keyName]
            RPM = data['RPM']
            Q = data['flowRate']
            condition = data['condition']
            dataFile = folder+condition+'\\'+keyName+'_data.pkl'
            videoPath = folder + condition + '\\' + keyName + '.mp4'
            
            # Skip the conditions to required for plot
            if (condition != conditionTarget):
                continue
            if (Q != QTarget) or (RPM not in RPMList):
                continue
            
            # Parse remaining data needed for plots    
            fps = data['fps']
            time = data['time']
            t0 = data['t0']
            eP = signal.medfilt(data['excessPerimeter'],3)
            aMax = data['aMax']
            aMean = data['aMean']
            eP[aMax<ePcutoff] = 1.0
            ePMax = 1.05*np.max(eP)
            tMax = 1.05*np.max(time-t0)
            contourMask = data['dilatedMask']
            # Choose frame numbers for plotting
            nFrames = len(time)
            frameNum = [int(round(time[z]*fps)) for z in range(nFrames)]
            
            # Load original data file and video file
            temp = glob.glob(dataFile)[0]
            with open(temp,'rb') as f:
                container = pkl.load(f)
            hMatrix = container['hMatrix']
            maskData = container['maskData']
            temp = maskData['diskMask']
            maskData['mask'] = temp
            center = maskData['diskCenter']
            R = maskData['diskRadius']
            pix_per_cm = container['maskData']['diskRadius']/diskScale
            cropRect = [center[0]-R,center[0]+R,center[1]-R,center[1]+R]
            vid = VF.get_video_object(videoPath)
            
            # Load rCrit and determine corresponding frame number
            rCritQList = rCritData['QList']
            rCritRPMList = rCritData['RPMList']
            i_Q = [k for k in range(len(rCritQList)) if QTarget == rCritQList[k]]
            j_RPM = [k for k in range(len(rCritRPMList)) if RPM == rCritRPMList[k]]
            rCrit = rCritData[condition][i_Q, j_RPM]
            num = [k for k in range(len(aMean)) if aMean[k+1] >= rCrit and aMean[k] <= rCrit]
            
            # Get contour and video frame at indicated time
            contour1 = container['data'][num]
            theta = container['theta'][num]
            # Adjust the rotation of the mask to match the stored data
            if RPM > 0:
                mask = IPF.rotate_image(contourMask,theta,
                                        size=np.shape(contourMask))
                mask = (mask == 255)
            # Remove contour points that are covered by the mask
            if removeMaskPoints:
                temp = mask[np.int16(contour1[:,1]),np.int16(contour1[:,0])]
                ind0 = np.argmax(temp == False)
                ind1 = len(temp) - np.argmax(temp[::-1] == False) - 1
                # Reorder temp if mask is covering beginning or end of temp
                if ((ind0 == 0) or (ind1 == len(temp)-1)) or (ind1-ind0>0.75*len(temp)):
                    halfInd = len(temp)/2
                    temp = np.concatenate((temp[halfInd:],temp[:halfInd]))
                    ind0 = np.argmax(temp == False)
                    ind1 = len(temp) - np.argmax(temp[::-1] == False) - 1
                    temp[ind0:ind1] = False
                    midInd = len(temp)-halfInd
                    temp = np.concatenate((temp[midInd:],temp[:midInd]))
                else:
                    temp[ind0:ind1] = False
                contour = contour1[temp]
            else:
                contour = contour1
            
            x = contour[:,0]
            y = contour[:,1]
            
            if aMax[i] < ePcutoff:
                R = aMax[i]*pix_per_cm
                C = 100
                if shiftCenter:
                    cen1 = trueCenter
                else:
                    cen1 = center
                x,y = Fun.generate_circle(R,cen1,N=C)
            
            frame = VF.extract_frame(vid,num,hMatrix,maskData)
            # Rotate video frame to match orientation of contour
            if noRotate:
                if aMax[i] >= ePcutoff:
                    x,y = Fun.rotate_points(x,y,theta,center,units='degrees')
            else:
                frame = IPF.rotate_image(frame,theta,size=np.shape(frame))
            
            
            # Make Plots
            plt.figure(1)
            plt.cla()
            plt.imshow(frame,cmap='gray')
            plt.plot(x,y,'r')
            plt.axis(cropRect)
            ax = plt.gca()
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)
            plt.tight_layout()
            plt.title('Critical Radius Reached - RPM = ' + str(RPM) + '; Q = ' + str(Q))
            
            saveName = 'image_%06i.png'%(num)
            plt.savefig(saveName,dpi=DPI)