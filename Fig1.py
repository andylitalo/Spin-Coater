## -*- coding: utf-8 -*-

import glob
import matplotlib.pyplot as plt
import cPickle as pkl
import numpy as np
from scipy import signal

import Functions as Fun
import VideoFunctions as VF
import ImageProcessingFunctions as IPF

# User parameters
folder = '..\\Data\\New Disk\\'
name = 'all_data.pkl'
# [500,1000,1500,2000,2500,3000,3250]
# [0,30,120,250,500,750,1000]
# ['Dry','Water','3mM SDS','6.5mM SDS','10mM SDS','Water with 10mM SDS']
Q = 500
RPM = 500
C = 'Dry'
tagList = [(C,Q,0),(C,Q,120),(C,Q,250),(C,Q,500)]
wetTarget = 0.9*15

# Plot properties
FS1 = 12 # axes labels
LW = 2 # Line width for plots
figscale = 5 # Set size of figure larger will produce relatively smaller font
colorList = ['b','g','r','c','m','k']
i = 0

   
# Cleanup existing windows
plt.close('all')
pi = np.pi
invert = True
    
# Load data
try:
  allData
except NameError:
    temp = glob.glob(folder+name)[0]
    with open(temp,'rb') as f:
        allData = pkl.load(f)
else:
    pass
    
# Compute plotting properties
nData = len(tagList)
dPhi = 2.0*np.pi/float(nData)
offset = -(pi/2 - dPhi)


# Initialize figure variables
plt.figure(1,figsize=(figscale,figscale))
plt.cla()
plt.figure(2,figsize=(figscale,figscale))
plt.cla()
handles1 = []
legList = []
legVal = []
image = [0]

for tag in tagList:

    # Parse the filename to get video info
    try:
        data = allData[tag]
    except:
        print tag, ' not found in stored data'
        continue

    # Parse the filename to get video info
    RPM = data['RPM']
    Q = data['flowRate']
    condition = data['condition']
    dataFile,videoPath = Fun.build_file_name(folder,condition,tag,data['fps'])
    
    # Parse remaining data needed for plots    
    fps = data['fps']
    time = data['time']
    aMax2 = data['aMax2']
    eP = signal.medfilt(data['excessPerimeter'],3)
    # Find frame number for plotting
    indWet = np.argmax(aMax2 > wetTarget)
    frameNum = int(time[indWet]*fps)
    
    # Load original data file and get contour at indicated time
    temp = glob.glob(dataFile)[0]
    with open(temp,'rb') as f:
        container = pkl.load(f)
    contour = container['data'][frameNum]
    theta = container['theta'][frameNum]%(-360)/180.*pi
    hMatrix = container['hMatrix']
    maskData = container['maskData']
    x = contour[:,0]
    y = contour[:,1]
    center = maskData['diskCenter']
    R = maskData['diskRadius']
    
    # Compute needed angles of rotation for processing data
    thetaDeg = theta/pi*180
    phiStart = (2*pi - dPhi - theta)%(2*pi)
    phiStop = (2*pi - theta)%(2*pi)
    shiftFactor = tagList.index(tag)
    angle = shiftFactor*dPhi + theta + offset
    
    
    # Center contour data and compute angle at each point
    x1 = x - center[0]
    y1 = y - center[1]
    phi = np.arctan2(y1,x1)
    phi[phi<0] += 2*pi
    # Cut to angles of interest
    if phiStart > phiStop:
        phiCrop = phiStart - 2*pi
        indPlot = (phi>=phiStart)+(phi<=phiStop)
    else:
        phiCrop = phiStart
        indPlot = (phi>=phiStart)*(phi<=phiStop)
    x = x[indPlot]
    y = y[indPlot]
    # Reorder so that the angle is continually increasing
    dx = x[1:] - x[:-1]
    dy = y[1:] - y[:-1]
    ds = np.sqrt(dx**2 + dy**2)
    ind = np.argmax(dx>dPhi/3.0*R)
    if ind > 0:
        ind += 1
    x = np.concatenate((x[ind:],x[:ind]))
    y = np.concatenate((y[ind:],y[:ind]))
    # Rotate points to desired position for plotting
    x,y = Fun.rotate_points(x,y,angle,center)
    
    ePMax = eP[indWet]
    if condition == 'Water with 10mM SDS':
        ePMax = 1.0
    
    legList += [r'%s; $E_p$ = %.2f' %(condition,ePMax)]

    plt.figure(1)
    if invert:
        plt.gca().invert_yaxis()
        invert = False
    handles1 += plt.plot(x,y,colorList[i]+'-',linewidth=LW)
    plt.axis('scaled')
    plt.pause(0.001)
    i += 1
    
    # Load the video frame from the original video file
    videoPath = glob.glob(videoPath)[0]
    vid = VF.get_video_object(videoPath)
    frame = VF.extract_frame(vid,frameNum,hMatrix,maskData)
    # Rotate video frame to match orientation of contour
    frame = IPF.rotate_image(frame,thetaDeg,size=np.shape(frame))
    if len(image) == 1:
        image = np.zeros_like(frame)
    xm,ym = Fun.generate_circle(R,center,
                                N=1000,t0=phiCrop,t1=dPhi)
    xm = np.concatenate(([center[0]],xm,[center[0]]))
    ym = np.concatenate(([center[1]],ym,[center[1]]))
    points = [(xm[z],ym[z]) for z in range(len(xm))]
    mask = IPF.create_polygon_mask(image,points)[0]
    tempFrame = IPF.rotate_image(frame*mask,-angle/np.pi*180,
                                 size=np.shape(frame))
    image = image + tempFrame
    
    plt.figure(2)
    plt.cla()
    plt.imshow(image,cmap='gray')
    plt.pause(0.001)                                

plt.figure(1)
tempR = maskData['maskRadius']
x,y = Fun.generate_circle(tempR,center)
handles1 += plt.plot(x,y,'k--',linewidth=LW)
plt.axis('scaled')
plt.pause(0.001)


legList += [r'Edge; $E_p$ = %.2f' %(1.0)]

plt.figure(1)
tempR *= 1.02
plt.axis([center[0]-tempR,center[0]+tempR,center[1]+tempR,
          center[1]-tempR])
ax = plt.gca()
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
ax.tick_params(axis='both',pad=-FS1)
ax.legend(legList,loc='center',numpoints=2,fontsize=FS1,frameon=False)
#ax.legend(tuple([handles1[i] for i in temp]), tuple([legList[i] for i in temp]),
#          loc='center left',bbox_to_anchor=(1, 0.5),numpoints=2)
plt.tight_layout(pad=0)
plt.axis('off')
plt.savefig('../Figures/Fig1a.pdf')
plt.pause(0.001)

plt.figure(2)
ax = plt.gca()
plt.axis([center[0]-tempR,center[0]+tempR,center[1]+tempR,
          center[1]-tempR])
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
ax.tick_params(axis='both',pad=-FS1)
plt.tight_layout(pad=0)
plt.savefig('../Figures/Fig1b.pdf')
plt.pause(0.001)