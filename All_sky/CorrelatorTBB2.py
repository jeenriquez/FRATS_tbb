#! /usr/bin/env python

##########################################################################
#
# Software to correlate LOFAR TBB data (h5 format) from a Single Station
#
# individual channelwidth can be choosen or can be calculated for differ.
# RM resolution
#
##########################################################################


from pycrtools import *
from pycrtools import metadata as md
import numpy as np
import pycrtools as cr
import sys
import os
import os.path
import re
import optparse
import numpy as np
import matplotlib.pyplot as plt
import pyfits
from pylab import *
from numpy import *
from math import *
from scipy import *
from numpy.fft import fft
import pickle
from scipy.special import iv
import time
from datetime import datetime
import threading
from tempfile import mkdtemp
import Queue


start = time.time()

getFFTData = False      # get FFTData instead of timeseries-data
var = False
divide = False

if getFFTData:
    divide = False

BSbig = 200         # for divide = True

#-------------get in files and choose options--------------------

usage = "usage: %prog [options] file0.h5 file1.h5 ..."
parser = optparse.OptionParser(usage=usage)

parser.add_option('-s', '--startblock',type="int",    help = 'startblock to read in the data', default = 0 )
parser.add_option('-b', '--blocksize', type="int",    help = 'blocksize as a factor of 512 (also f.e. 0.1 for 51 channels is possible )', default = 1)
parser.add_option('-n', '--nblocks',   type="int",    help = 'number of blocks to read in; 0 means all', default = 2)
parser.add_option('-f', '--startchan', type="int",    help = 'startfrequency in MHz ', default = None)
parser.add_option('-e', '--endchan',   type="int",    help = 'endfrequency in MHz ', default = None)
parser.add_option('-o', '--outfile',   type="str",    help = 'outfile for the cm', default = "correlation")
parser.add_option('-r', '--rcumode',   type="int",    help = 'RCU mode; 3: LBA, 5-7: HBA ', default = 5)

(options, args) = parser.parse_args()

startblock = int(options.startblock)
wide = options.blocksize
nblocks = int(options.nblocks)
startfreq = float(options.startchan)
endfreq = float(options.endchan)
outfile = str(options.outfile)
rcumode = int(options.rcumode)

if not startfreq  or not endfreq :
    raise ValueError('Need startchan and endchan within the appropriate band range (eg. 100-200 MHz for HBA - rcumode=5)')

if wide == "var":
    var = True
#else:
#   wide = int(wide)

#-------------calculating needed frequencies for given RM--------

if var:
    print "calculation of the frequencies \n"
    RM = 100            # max RM that is expected
    deltaphi = 180          # max deltaphi (depol.) within one SB
    fre = startfre
    freliste = []
    BSliste = []

    const = 0.5*deltaphi/(3*3*10000*RM)*2*pi/360
    deltav0 = const*fre**3
    wide = int(round(0.2/deltav0))  # equal to BS0

    while math.floor(fre) <= endfre:
        deltav = const*fre**3
        BS = 0.2/deltav
        freliste.append(fre)
        fre += deltav
        BSliste.append(round(wide/BS))

    lenBSliste = len(BSliste)
    print " - ready -\n"

#----------------------------------------------------------------

#startfre = int(startfreq/200*1024*wide)
#endfre = int(endfreq/200*1024*wide)

share = 1
if divide:
    if wide <= BSbig:
        share = int(round(BSbig/wide))
        wide *= share

blocksize = int(wide * 1024)

files = args
count = 0
for n in files:
    count += 1

data = [array]*count
shiftmax = 0
f = [array]*count

#-------------read in RCUs---------------------------------------

for i in range(0,count,1):
    data[i] = cr.open(files[i],blocksize)
    f[i] = cr.TBBData(files[i])

    if i == 0:
        rcun = data[i]["DIPOLE_NAMES"]
        nant = len(rcun)
        shifts = zeros((count,nant),dtype = int)

    if i != 0:
        rcun2 = data[i]["DIPOLE_NAMES"]
        rcun = hstack((rcun, rcun2))
        del rcun2

#-------------read in antenna positions -------------------------

for i in range(0,count,1):
    data[i] = cr.open(files[i],blocksize)
    f[i] = cr.TBBData(files[i])

    if i == 0:
        antposp = data[i]["ITRF_ANTENNA_POSITION"].toNumpy()
        nant = len(antposp)

    if i != 0:
        antposp2 = data[i]["ITRF_ANTENNA_POSITION"].toNumpy()
        antposp = vstack((antposp, antposp2))
        del antposp2

#-------------shifts between differnt RCUs-----------------------

    shifts[i] = data[i]["SAMPLE_NUMBER"]

    localmax = max(shifts[i])
    if localmax > shiftmax:
            shiftmax = localmax

for i in range(0,count,1):
    shifts[i] = ((shiftmax-shifts[i]))

    if getFFTData:
        shifts[i] /= blocksize

rculiste = [array]*len(rcun)
for rcu in range(0,len(rcun),1):
    rculiste[rcu] = (int(rcun[rcu])%1000)

rcun = rculiste
#antnr = max(rcun)+1
antnr = 0
for d in data:
    antnr += d["NOF_SELECTED_DATASETS"]


#----------------------------------------------------------------

if rcumode == 3: nyquistZone = 1
if rcumode == 5: nyquistZone = 2
if rcumode == 6: nyquistZone = 3
if rcumode == 7: nyquistZone = 3

lenfre = blocksize/2 + 1

if nblocks == 0:
    nblocks = int((data[0]["DATA_LENGTH"][0]-shiftmax/1024-startblock-10000)/1024/wide)

if not var:
#   freliste = (cr.open(files[0], blocksize/share)).getFrequencies().toNumpy()[startfre:endfre]
    freliste = (cr.open(files[i],blocksize/share)).empty("FREQUENCY_DATA")
    cr.hFFTFrequencies(freliste, data[i]["SAMPLE_FREQUENCY"][0], nyquistZone)
    freliste = freliste.toNumpy()
    freliste /= 1000000
    if (startfreq != 0) and (endfreq != 0):
    #   if rcumode == 5:
    #       endfre = (abs(freliste-startfreq)).argmin()
    #       startfre = (abs(freliste-endfreq)).argmin()
    #   else:
        startfre = (abs(freliste-startfreq)).argmin()
        endfre = (abs(freliste-endfreq)).argmin()
        freliste = freliste[startfre:endfre]

stationname = data[0]['STATION_NAME'][0]
timestamp = (data[0])["TIME"][0]
timedate = datetime.fromtimestamp(timestamp)        # real time
timeutc = datetime.utcfromtimestamp(timestamp)      # UTC time
timetuple = time.gmtime(timestamp)          # time tuple

del shiftmax, files, rculiste

#-------------coeffizients for PolyphaseFilter-------------------
# (S.Faint, W.Read, 2003)

if not getFFTData:
    #print "\nStart calculating filter coefficients\n"

    a = 16.0/pi #for Kaiser window

    if divide:
        blocksize /= share
    h = cr.hArray(complex, dimensions = (1,1,blocksize))

    for n in range(0,blocksize,1):
    # Kaiser-Window:
    #   wn = iv(0,pi*a*sqrt(1-(2.0*n/blocksize-1)*(2.0*n/blocksize-1)))/(iv(0,pi*a))
    # Hamming-window:
    #   wn = 0.54+0.46*math.cos(2*n/blocksize*pi)
    # Blackman-Window:
        wn = 0.42 - 0.5*math.cos(2.0*n/blocksize*pi) + 0.08*cos(4.0*n/blocksize*pi)

        h[0,0,n] = wn * sinc((n-float(blocksize)/2+0.5)/blocksize)
        del wn

    if divide:
        blocksize *= share

    h = h.toNumpy()

    #print " - ready -\n"

#-------------reading in the data--------------------------------

class readin(threading.Thread):
    def __init__(self, i, block, offset):
        self.i = i
        self.block = block
        self.offset = offset
        threading.Thread.__init__(self)
    def run(self):
        if getFFTData:
            self.offset = cr.hArray(self.block + shifts[self.i])
            data[self.i].getFFTData(fftdata[self.i*nant:(self.i+1)*nant], self.offset)
        else:
            self.offset = cr.hArray(self.block + shifts[self.i])
            self.offset = cr.hArray((self.block)*blocksize + shifts[self.i]+startblock)
#                        data[self.i]["BLOCKSIZE"] = blocksize
#                        print "BLA", fftdata.shape(), self.i*nant, (self.i+1)*nant, self.offset
#                        data[self.i].getTimeseriesData(fftdata[self.i*nant:(self.i+1)*nant], self.offset)
            cr.hReadTimeseriesData(fftdata[self.i*nant:(self.i+1)*nant], self.offset, blocksize, f[self.i])

#-------------correlate the data---------------------------------

class correlation(threading.Thread):
    def __init__(self, fre, fftdatan):
        self.fre = fre
        self.fftdatan = fftdatan
        threading.Thread.__init__(self)
    def run(self):
        acm[:,:,self.fre] += dot((self.fftdatan.transpose()).conjugate(), self.fftdatan)

#-------------main slope for reading and correlation-------------

class main(threading.Thread):
    def __init__(self):
        threading.Thread.__init__(self)
    def run(self):
        for block in range(0, nblocks,1 ):

            ready = [array]*count
            for i in range(0, count,1):
                offset = 0
                threadread = readin(i, block, offset)
                threadread.start()
                ready[i] = threadread

            for threadread in ready:
                threadread.join()
            del ready

            fftdatan = fftdata.toNumpy()

            if block != 0:
                for threadcorr in corry:
                    threadcorr.join()
                del corry

            fftdatan = hsplit(fftdata.toNumpy(),share)

            if not getFFTData:
                fftdatan *= h
                fftdatan = (fft(fftdatan))[:,:,startfre:endfre]
            else:
                fftdatan = fftdatan[:][:][startfre:endfre]

            corry = [array]*len(freliste)
            if var:
                diffBS = 0
                for fre in range(0,lenBSliste,1):
                    diffBS += BSliste[fre]
                    threadcorr = correlation(fre, fftdatan[:,:,int(diffBS)])
                    threadcorr.start()
                    corry[fre] = threadcorr
            else:
                for fre in range(0,endfre-startfre,1):
                    threadcorr = correlation(fre,fftdatan[:,:,fre])
                    threadcorr.start()
                    corry[fre] = threadcorr

            if block == (nblocks-1):
                for threadcorr in corry:
                    threadcorr.join()

            print 'block: {0}/{1} \r'.format(repr(block+1).rjust(len(str(nblocks))),nblocks),
            sys.stdout.flush()
#-------------start of main programm-----------------------------

fftdatan = zeros((share,antnr, len(freliste)),dtype = complex)

if getFFTData:
    fftdata = cr.hArray(complex, dimensions = (antnr,lenfre))
else:
    fftdata = cr.hArray(float, dimensions = (antnr,blocksize))

acm = zeros((antnr,antnr,len(freliste)),dtype = complex)

print "\n\nNumber of input files:    %i" % (count)
print "Blocksize:                %i" % (blocksize/wide)
print "Startblock:               %i " % (startblock)
print "Startfrequency:           %6.2f MHz" % (startfreq)
print "Endfrequency:             %6.2f MHz" % (endfreq)
print "------------------------------------------------"

print "\nSTART CORRELATION\n"

threadmain = main()
threadmain.start()
threadmain.join()

del fftdata, fftdatan, data, h
if var:
    del BSliste

##------------sorting of the correlation matrix------------------
#
#corrm = zeros((antnr,antnr,len(freliste)),dtype = complex)
#
#i1 = 0
#for rcu1 in rcun:
#   i2 = 0
#   for rcu2 in rcun:
#       corrm[rcu1][rcu2] = (acm[i1][i2]).conjugate()
#       i2 += 1
#   i1 += 1
#
#del acm
#
#print "\n - ready - \n"
corrm = acm.conjugate()
if rcumode == 5:
    corrm = corrm[:,:,::-1]


#-------------write correlation matrix in a file-----------------

print "Start writing \n"


of = open(outfile+".pickle",'w')
antposp = md.convertITRFToLocal(cr.hArray(antposp)).toNumpy()
pickle.dump(antposp,of,1)
pickle.dump(timedate,of,1)
pickle.dump(timeutc,of,1)
pickle.dump(timetuple,of,1)
pickle.dump(stationname, of,1)
pickle.dump(freliste,of,1)


np.save(outfile+".npy", corrm)

nant = corrm.shape[0]

print "data sucessfully written to file:", outfile + ".pickle", "and", outfile + ".npy"

#outfile = ("acm_"+str(nblocks)+"_2048.out")

end = time.time()

print "\nTime:   ", round((end-start)/60,2), "min,   ", round(end-start), " sec"


#-------------code for reading in the corrmatrix-----------------

"""

dataset = open(dataset)
timedate = pickle.load(dataset)
timeutc = pickle.load(dataset)
timetuple = pickle.load(dataset)
station = pickle.load(dataset)
freliste = pickle.load(dataset)
u = pickle.load(dataset)

data = pickle.load(dataset)
if u != 0:
    for uc in range(0,u-1,1):
        data2 = pickle.load(dataset)
        data = vstack((data,data2))
        del data2

#data as numpy array with [ # rcus, # rcus, # channels]

print "\n File informations \n"
print "Observing Time      ", timedate
print "Obs Time in UTC     ", timeutc
print "Station             ", station
print "Startfrequency      ", round(freliste[0],2),"MHz"
print "Endfrequency        ", round(freliste[-1],2),"MHz"
print "Number of channels  ", len(freliste), "\n"

"""

