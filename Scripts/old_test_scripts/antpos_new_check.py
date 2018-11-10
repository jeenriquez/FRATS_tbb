'''Script to check that the new relative antenna position function is correct.
'''

import pycrtools as cr
from pycrtools import metadata as md
import numpy as np
import pdb;# pdb.set_trace()
import glob
import Beam_Tools as bt
import os

LOFARDATA=os.environ["LOFARDATA"].rstrip('/')+'/'
filenames=glob.glob(LOFARDATA+ '*.h5')

filter = 'HBA0'
diff = 0
Two_d = 1
if Two_d:
    filenames = [filenames[0]]

xpos1=cr.hArray(float,(len(filenames),96))
ypos1=cr.hArray(float,(len(filenames),96))
zpos1=cr.hArray(float,(len(filenames),96))

stations=[]


for i,file in enumerate(filenames):
    f=cr.open(file)
    f['ANTENNA_SET']=filter

    if diff:
        antpos = f.getRelativeAntennaPositions() - md.getRelativeAntennaPositionsNew(f['STATION_NAME'][0],filter)
    else:
        antpos = md.getRelativeAntennaPositionsNew(f['STATION_NAME'][0],filter)

    #X-axis positions
    xpos=cr.hArray(float,[96])
    xpos[...].copy(antpos[...,0])
    xpos1[i]=xpos

    #Y-axis positions
    ypos=cr.hArray(float,[96])
    ypos[...].copy(antpos[...,1])
    ypos1[i]=ypos

    #Z-axis positions
    zpos=cr.hArray(float,[96])
    zpos[...].copy(antpos[...,2])
    zpos1[i]=zpos

    stations.append(f['STATION_NAME'][0])

ant_num=range(96)
cr.plt.ion()

if diff:
    cr.plt.rcParams['font.size']=20
    cr.plt.rcParams['axes.titlesize']=30
    cr.plt.suptitle('Relative Antenna Possition Comparison. ANTENNA_SET = ' + filter)

    cr.plt.subplot(2,2,1)
    cr.plt.plot(ant_num,xpos1.Transpose().toNumpy())
    cr.plt.xlabel("Antenna Number ")
    cr.plt.ylabel('Difference [m]')
    cr.plt.title("(X-Axis)")

    cr.plt.subplot(2,2,2)
    cr.plt.plot(ant_num,ypos1.Transpose().toNumpy())
    cr.plt.xlabel("Antenna Number ")
    cr.plt.ylabel('Difference [m]')
    cr.plt.title("(Y-Axis)")

    cr.plt.subplot(2,2,3)
    cr.plt.plot(ant_num,zpos1.Transpose().toNumpy())
    cr.plt.xlabel("Antenna Number ")
    cr.plt.ylabel('Difference [m]')
    cr.plt.title("(Z-Axis)")

    cr.plt.subplot(2,2,0)
    cr.plt.plot(ant_num,zpos1.Transpose().toNumpy()*0)
    cr.plt.legend(stations,loc=3)

if Two_d:
    cr.plt.plot(xpos1[0].vec(),ypos1[0].vec(),'o')
    inner=[17,19,29,31]
    inner=cr.hArray(int,len(inner),inner)
    cr.plt.plot(xpos1[0,list(inner)].vec(),ypos1[0,list(inner)].vec(),'o')
    cr.plt.plot(xpos1[0,list(inner+48)].vec(),ypos1[0,list(inner+48)].vec(),'o')



