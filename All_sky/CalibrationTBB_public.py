#!/usr/bin/env python
# -*- coding: utf-8 -*-

##########################################################################
#                                                                        #
#   Software to calibrate LOFAR data from a Single Station (All-Sky)     #
#   written by Jana Koehler (MPIfR, Bonn)                #
#                                    #
#   contact: jkoehler@mpifr-bonn.mpg.de or jana@mpifr-bonn.mpg.de    #
#                                    #
##########################################################################

from pylab import *
from numpy import *
from scipy import *
import time
import os
from optparse import OptionParser
import re
import sys
import pickle
import pyfits
import threading
# for plotting in galactic coordinates:
import matplotlib.pyplot
import matplotlib.backends.backend_agg

#################################################################

# give the directory the AntennaArrays.conf files:
LOFARSOFT = os.environ["LOFARSOFT"].rstrip('/') + '/'
antennafiledir = LOFARSOFT +"data/lofar/StaticMetaData/AntennaArrays/"

# choose which polarizations you want to get plottet:
# for galactic maps it is always Stokes I
# choose from I,Q,U,V,XX,YY,XY,YX
map_Stokes = ['I']

# if you want to use a beam model
#Jonesfile = "/.aux_mnt/pc20021a/Jones_512_n.npy"


showimages = True   # plot/show the image
saveimages = False  # save image(s) of skymaps or spectra in a .png file
savedata = False        # save data of skymaps... in a pickle file
plotsources = True      # plot source position in the maps
sourcelistplot = (["CasA","CygA","TauA","VirA","HerA","HydA","Sun","PSR B0329+54"])

if saveimages or savedata:
    directory = "test_galac" #os.getcwd()       # choose a directory, where you want to save the images/data
# in this case it writes everything to the directory, where you are running the code



#################################################################
###----------               OPTIONS                 ----------###
#################################################################

if len(sys.argv) < 2:
    print '\n   DATASET FILE IS NEEDED\n\nUse -h for help\n'
    exit()

filename = str(sys.argv[1])     # filename

usage = "usage: %prog corr_matrix_file.out [options] ..."
parser = OptionParser(usage=usage)

parser.add_option("-n", "--no-calibration",
                  action = "store_false", dest="calibrate", default=True,
                  help = "do not calibrate the data before imaging")
parser.add_option('-l', '--blflag', help = 's for just short baselines; l for just long baselines; 0 for all baselines', default = 0)
parser.add_option("-b", "--beam",
                  action = "store_true", dest="beam", default=False,
                  help = "use a beam model")
parser.add_option('-s', '--startchannel', help = 'Startchannel for the Calibartion; all for all; you also can write for example -s 53.2MHz, if you want a special frequency instate of a channel ', default = 251)
parser.add_option('-e', '--endchannel', help = 'Endchannel, if none given, only take startchannel for Calibration; you also can write for example -e 53.2MHz, if you want a special frequency instate of a channel', default = 0)
parser.add_option("-p", "--phasecal",
                  action = "store_false", dest="phasecal", default=True,
                  help = "no phase calibration will be done")
parser.add_option("-a", "--ampphacal",
                  action = "store_false", dest="ampphacal", default=True,
                  help = "no ampl-phase calibration will be done")
parser.add_option("-i", "--integratefrequencies",
                  action = "store_true", dest="integratefrequencies", default=False,
                  help = "integration over choosen frequencies will be done")
parser.add_option('-r', '--output', help = 'gain for saving Gains; map for making the Skymap; mapf for saving Skymap in FITS format; galactic for making Skymap in Galactic Coordinates; spectra for Autocorrelation Spectra from all RCUs', default = "map")
parser.add_option('-d', '--polarisation', help = '0 for XX;, 1 for YY; 2 for total intensity', default = 2)
parser.add_option('-f', "--find-sources",action = "store_true", dest="find_sources", default=False, help = "Search for sources in the sky, instead of using Ateam sources")
parser.add_option('-o', '--outfile',    help = 'name of the outfile', default = 'automatic')


(options, args) = parser.parse_args()

start_chan = (options.startchannel)
end_chan = (options.endchannel)
output = str(options.output)
out = str(options.outfile)
pol = int(options.polarisation)
blflag = str(options.blflag)


if filename[0] == '-':
    print '\n   DATASET FILE IS NEEDED\n\nUse -h for help\n'
    exit()

if out == 'automatic':
    if output == 'gain' or output == 'amp' or output == 'spectra':
        outfile = (str(output)+'.out')
    elif output == 'galactic':
        outfile = (str(output)+str(start_chan)+'.out')
    elif output == 'map':
        outfile = ('skymap_'+str(start_chan)+'.png')
    elif output == 'mapf':
        outfile = ('skymap_'+str(start_chan)+'.fits')
else:
    outfile = out

if (showimages == False) and (saveimages == False) and (savedata == False) and (output != 'gain') and (output != 'spectra'):
    print "With the choises you made, there will be now data/images saved or shown. Therefor (in theory) it doesn't make any sense to run the programm and it will exit here. But if you still want to run the code in this way (for what reason ever) go to line $$, remove the exit() and go on with having fun..."
    exit()

fre = start_chan

#################################################################
###----------              DATA INPUT               ----------###
#################################################################

#-------------reading in the dataset-----------------------------
# format: data[rcu, rcu, fre (in MHz)]
"""
datase = os.path.basename(dataset)
nrcu = int(datase[24:27])
date = datase[0:15]
freliste = arange(0,512,1)*100./512
dataset = open(dataset)
data = reshape(fromfile(dataset,dtype = complex128),(512,nrcu,nrcu))
data = rollaxis(data,1,0)
data = rollaxis(data,2,1)
"""

dataset = open(filename+".pickle") #open(dataset)
antposp = pickle.load(dataset)
timedate = pickle.load(dataset)
timeutc = pickle.load(dataset)
timetuple = pickle.load(dataset)
station = pickle.load(dataset)
freliste = pickle.load(dataset)     # in MHz

data = np.load(filename+".npy")
print data.shape

year = float(timetuple[0])
month = float(timetuple[1])
day = float(timetuple[2])
hour = float(timetuple[3])
minu = float(timetuple[4])
sec= float(timetuple[5])


print "\n\nFILE INFORMATIONS: \n"
print "Observing Time      ", timedate
print "Obs Time in UTC     ", timeutc
print "Station             ", station
print "Startfrequency      ", round(freliste[0],2),"MHz"
print "Endfrequency        ", round(freliste[-1],2),"MHz"
print "Number of Channels  ", len(freliste)
print "-------------------------------------------\n"

#----------------------------------------------------------------

timeutc = str(year)+'-'+str(month)+'-'+str(day)+' '+str(hour)+str(minu)+str(sec)

if start_chan == 'all':
    Calall = True
elif str(start_chan)[-3:] != 'MHz':
    start_chan = int(start_chan)
    Calall = False
    if end_chan == 0:
        end_chan = int(start_chan)+1
    else:
        end_chan += 1

if str(start_chan)[-3:] == 'MHz':
    start_chan = (abs(freliste-float(start_chan[:-3]))).argmin()
    Calall = False
    if end_chan == 0:
        end_chan = int(start_chan)+1
    else:
        end_chan = (abs(freliste-float(end_chan[:-3]))).argmin()+1


if start_chan > len(freliste):
    parser.error("Start Channel is out of range\nUse -h for help\n")
if end_chan > len(freliste):
    parser.error("End Channel is out of range\nUse -h for help\n")

#-------------reading in a Jones Matrix--------------------------
# Jones matrix [channel, 'polarisation', l, m] - # [512, 4, 200, 200]
# Calculated from a code provided by Tobia Carozzi

if options.beam:
    Jones = np.load(Jonesfile)
#   beamy = np.load("beam_Y_260_cascyg.npy")
#   beamx = rot90(beamy)


# if you want to have a look at the beams, then use the stuff below:

#jonesx = (abs(Jones[channel,0].real))**2 + (abs(Jones[channel,1].real))**2
#jonesy = (abs(Jones[channel,2].real))**2 + (abs(Jones[channel,3].real))**2
#imshow(jonesx,aspect='auto')
#plot(freliste,(abs(Jones[:,0,100,100].real))**2 + (abs(Jones[:,1,100,100].real))**2)
#plot(freliste,(abs(Jones[:,0,30,80].real))**2 + (abs(Jones[:,1,30,80].real))**2)

#colorbar()
#show()
#exit()


#################################################################
###----------         STATION POSITIONS             ----------###
#################################################################

#-------------Position of station and antennas-------------------
#station = str("SE607")
#if str(station) == "SE607":
#   antennafilename = "SE607-AntennaArrays.conf.LBA.dat"
#else:
#
antennafilename = antennafiledir+station+"-AntennaArrays.conf"

if (0 <= freliste[0]) and( freliste[-1] <= 100):
    field = 'LBA'
    rcumode = 3
elif (100 <= freliste[0]) and( freliste[-1] <= 200):
    field = 'HBA'
    rcumode = 5
    print "use HBA mode 5"
elif (160 <= freliste[0]) and( freliste[-1] <= 240):
    field = 'HBA'
    rcumode = 6
    print "use HBA mode 7"
elif (200 <= freliste[0]) and( freliste[-1] <= 250):
    field = 'HBA'
    rcumode = 7
    print "use HBA mode 7"
else:
    print "wrong frequency values in the data file"

nrcu = len(data)

position = open(antennafilename).readlines()

count = 0
for lines in position:
    sline=lines.split()
    if sline[0] == field : break
    count += 1

if str(station) == "SE607":
    longitude = 11.92
    latitude = 57.4
else:
    station = position[count+1].split()
    longitude = float(station[2])
    latitude = float(station[3])

#xpos = zeros((nrcu/2))
#ypos = zeros((nrcu/2))
#for line in range(0, nrcu/2, 1):
#   token = position[line+3+count].split()
#   xpos[line] = (float(token[0]))
#   ypos[line] = (float(token[1]))

xpos = antposp[0::2,0]
ypos = antposp[0::2,1]

baselines = zeros((nrcu,nrcu))

for rcu1 in range(0, nrcu/2, 1):
    for rcu2 in range(0, nrcu/2,1):
        baselines[2*rcu1,2*rcu2] = sqrt((xpos[rcu2]-xpos[rcu1])**2 + (ypos[rcu2]-ypos[rcu1])**2)
        baselines[2*rcu1+1,2*rcu2+1] = baselines[2*rcu1,2*rcu2]

#plot(xpos,ypos,'ko')
#show()
#exit()

xpos = reshape(xpos,(1,nrcu/2))
ypos = reshape(ypos,(1,nrcu/2))


#-------------general settings of parameters---------------------

# settings for the Calibration routine
flag_max = 4.
flag_min = 0.3

flagging = False
autoc = True
blflag = True

max_baselenght = 0.0
maxval = 0.8

sou = (["CasA","CygA"])#,"Tau","VirA","HerA" ])#, "HydA","PerA"])


#################################################################
###----------          SOURCE POSITIONS             ----------###
#################################################################

# [source, rascension(deg, min, sec), declination(deg, min, sec), a, b, c, Intensity(38MHz)]
# see Baars et al., 1977 :
# log S[Jy] = a + b * log v[MHz] + c * (log v[MHz])**2

#            source|   rascension    |    declination    |   a  |    b  |    c  |   I
sourcelis = array([ ("CasA", 23.0, 23.0, 27.9,  58.0, 48.0,  42.0, 5.625, -0.634, -0.023, 37200),
                    ("CygA", 19.0, 59.0, 28.3,  40.0, 44.0,   2.0, 4.695,  0.085, -0.178, 22000),
                    ("TauA",  5.0, 34.0, 32.0,  22.0,  0.0,  52.0, 3.915, -0.299,    0.0,  2430),
                    ("VirA", 12.0, 30.0, 49.4,  12.0, 23.0,  28.0, 5.023, -0.856,    0.0,  4000),
                    ("HerA", 16.0, 51.0, 37.7,   4.0, 58.0, 33.88, 4.963, -1.052,    0.0,  1800),
                    ("PSR B0329.54", 3.0 ,32.0 ,59.37, 54.0, 34.0, 44.9, 0.0, 0.0,    0.0,  0),
                    ("HydA",  9.0, 18.0,  5.7, -12.0, -5.0, -44.0, 4.497,  -0.91,    0.0,  1200) ])
            #       ("Jup", 4.58, 21.2,0,0,0,0,0,0,0,0) ])#,
            #       ("Jup", 4.997, 21.88,0,0,0,0,0,0,0,0) ])
                    #("PerA", 3.0, 19.0, 48.1, 41.0, 30.0, 42.0, 340),
                    #("SgrA", 17.0, 45.0, 40.0, -29.0, 0.0, -28.0, 1) ])

#----------------------------------------------------------------

if month <= 2:
    year -= 1
    month += 12

H = hour/24 + minu/1440 + sec/86400
JD = int(365.25*(year+4716))+int(30.6001*(month+1))+day+H+2+int((int(year/100))/4)-int(year/100)-1524.5
nJD = JD - 2451545.0        #time since standard equinox J2000
T = float(nJD)/36525
GWhourangle = mod(280.46061838+13185000.77*T+T*T/2577.765-T*T*T/38710000, 360)  #Greenwich hourangle
vernalequinox = mod(GWhourangle + longitude, 360)   #hourangle of vernal equinox


#-------------Calculate l and m for given Sources----------------

def source_l_m(so):

    if str(so) == "Sun":
        L = mod(280.460 + 0.9856474*nJD, 360)     #mean ecliptic longitude of the sun; aberration is included
        g = mod(357.528 + 0.9856003*nJD, 360)   #mean anormaly
        delta = mod((L + 1.915*math.sin(pi*g/180)+0.02*math.sin(2*pi*g/180)), 360)              #ecliptic longitude of the sun
        ecliptic = (23.439 - 0.0000004*nJD)*pi/180

        rascension = math.atan(math.cos(ecliptic)*math.sin(pi/180*delta)/cos(pi/180*delta))*180/pi
        if cos(pi/180*delta) < 0:
            rascension += 180
        declination = math.asin(math.sin(ecliptic)*sin(pi/180*delta))
        I = 1000

    else:
        for source in range(0,len(sourcelis),1):
            if str(so) == str(sourcelis[source][0]):
                break

        rascension = (float(sourcelis[source][1])+float(sourcelis[source][2])/60+float(sourcelis[source][3])/3600)*15
        declination = (float(sourcelis[source][4])+float(sourcelis[source][5])/60+float(sourcelis[source][6])/3600)*pi/180

        if str(so) == "Jup":
            rascension = float(sourcelis[source][1])*15
            declination = float(sourcelis[source][2])*pi/180
            I = 100000.0

    hourangle = (vernalequinox - rascension)*pi/180

    Elevation = math.asin(cos(declination)*cos(hourangle)*cos(latitude*pi/180)+sin(declination)*sin(latitude*pi/180))*180/pi
    Azimuth = math.atan(sin(hourangle)/(sin(latitude*pi/180)*cos(hourangle)-cos(latitude*pi/180)*tan(declination)))*180/pi+180
    if (cos(hourangle)*sin(latitude*pi/180)-tan(declination)*cos(latitude*pi/180))<0:
        Azimuth += 180
    if Azimuth < 0:
        Azimuth += 360
    Azimuth = mod(Azimuth, 360)

    #Conversion into direction cosines

    if Elevation >= 0:
        l = (cos((90-Azimuth)*pi/180))*(cos((Elevation)*pi/180))
        m = (sin((90-Azimuth)*pi/180))*(cos((Elevation)*pi/180))

        if str(so) != "Sun" and "Jup":
            I = (10**(float(sourcelis[source][7])+float(sourcelis[source][8])*log10(f/10**6)+float(sourcelis[source][9])*(log10(f/10**6))**2))#/2#/9216
        if str(so) == "Jup":
            I = 100000.0
        if str(so) == "CasA":
            n = (0.97 - 0.3 * log10(f/10**9))/100
            n2 = (-1.81+5.9*10**(-2)*log(f))/100

            for i in range(1966,2012,1):
                I += I*n2
    else:
        l = 'NaN'
        m = 'NaN'
        I = 'NaN'

    return l,m,I,Azimuth, Elevation



#################################################################
###----------               FLAGGING                ----------###
#################################################################

#-------------Flagging in Baselines------------------------------

# fits a curve through the amplitudes of the Visibilities as a function
# of the baseline lenght and then finds big outlyers

def FLAGGING(datfl):
    print "\n- Start Flagging in Baselines -"

    coeffs = polyfit(baselines[::2,::2].flatten(),(abs((datfl[::2,::2]).flatten())),14)
    x = polyval(coeffs,baselines[::2,::2].flatten())
    coeffs = polyfit(baselines[1::2,1::2].flatten(),(abs((datfl[1::2,1::2]).flatten())),14)
    y = polyval(coeffs,baselines[1::2,1::2].flatten())
    coeffs = polyfit(baselines[::2,::2].flatten(),(abs((datfl[::2,1::2]).flatten())),14)
    z = polyval(coeffs,baselines[::2,::2].flatten())

    (datfl[::2,::2])[abs(datfl[::2,::2]) > flag_max*x.reshape(nrcu/2,nrcu/2)] = 0.0
    (datfl[::2,::2])[abs(datfl[::2,::2]) < flag_min*x.reshape(nrcu/2,nrcu/2)] = 0.0
    (datfl[1::2,1::2])[abs(datfl[1::2,1::2]) >= (flag_max*y.reshape(nrcu/2,nrcu/2))] = 0.0
    (datfl[1::2,1::2])[abs(datfl[1::2,1::2]) <= (flag_min*y.reshape(nrcu/2,nrcu/2))] = 0.0
    (datfl[::2,1::2])[abs(datfl[::2,1::2]) > flag_max*z.reshape(nrcu/2,nrcu/2)] = 0.0
    (datfl[::2,1::2])[abs(datfl[::2,1::2]) < flag_min*z.reshape(nrcu/2,nrcu/2)] = 0.0

    (datfl[::2,::2])[abs(datfl[1::2,1::2]) == 0.0] = 0.0
    (datfl[1::2,1::2])[abs(datfl[::2,::2]) == 0.0] = 0.0
    (datfl[1::2,::2])[abs(datfl[::2,1::2]) == 0.0] = 0.0

    return datfl

#-------------Flagging by using the Autocorrelations-------------

# takes an average over the autocorrelation spectrum from all RCUs
# and finds outlyers
# centre frequency is used for the average

def flag_from_autocorr(data):
    print "\n- Start Flagging bad RCUs using Autocorrelations -"

    datas_cal = (data[:,:,len(freliste)/2]).diagonal()
    aver_XX = median(datas_cal[::2].real)
    aver_YY = median(datas_cal[1::2].real)
    rcuflag = zeros((len(datas_cal)),dtype = bool)

    rcuflag[::2] = (datas_cal[::2] < aver_XX*0.5) + (datas_cal[::2] > aver_XX*4.0) + (datas_cal[1::2] < aver_YY*0.5) + (datas_cal[1::2] > aver_YY*4.0)
    rcuflag[1::2] = rcuflag[::2]

    data[rcuflag] = 0.0
    data[:,rcuflag] = 0.0

    flagged_rcus = (arange(nrcu))[rcuflag]
    for rcu in range(0,len(flagged_rcus),2):
        print "   RCU %3i flagged, because of bad values" % (flagged_rcus[rcu])

    return data, rcuflag


#-------------Find RFI peaks in Spectrum-------------------------

# mit minus, statt mit division versuchen ....
"""
bpb1 = (abs(freliste-30.1)).argmin()    #devide Spectrum in 2 pieces to
bpe1 = (abs(freliste-56.0)).argmin()    #find the best fit for the Spectrum
bpb2 = (abs(freliste-54.0)).argmin()
bpe2 = (abs(freliste-82.0)).argmin()

def RFI_peaks(frelis, datadig):
    bpb1 = (abs(frelis-30.1)).argmin()  #devide Spectrum in 2 pieces to
    bpe1 = (abs(frelis-56.0)).argmin()  #find the best fit for the Spectrum
    bpb2 = (abs(frelis-54.0)).argmin()
    bpe2 = (abs(frelis-82.0)).argmin()

    #if not options.StationCorrelator_data:
#   datadig = dat.diagonal()    # just the autocorrelations

    datadiga = (datadig[bpb1:bpe1].mean(axis=1))
    coeffs = polyfit(frelis[bpb1:bpe1],log10(abs(datadiga)),9)
    y1 = 10**polyval(coeffs,frelis[bpb1:bpe1])  #fit for the first part

    datadiga = (datadig[bpb2:bpe2].mean(axis=1))
    coeffs = polyfit(frelis[bpb2:bpe2],log10(abs(datadiga)),9)
    y2 = 10**polyval(coeffs,frelis[bpb2:bpe2])  #fit for the second part

    bandpass = zeros((len(frelis)))
    bandpass[bpb1:bpb2] = y1[:bpb2-bpb1]        #combining everything
    bandpass[bpb2:bpe1] = (y1[bpb2-bpb1:]+y2[:bpe1-bpb2])/2
    bandpass[bpe1:bpe2] = y2[bpe1-bpb2:]
    datadiga = (datadig[bpb1:bpe2].mean(axis=1))
    meanbp = median(datadiga/bandpass[bpb1:bpe2])
    RFIpeaks = ones((len(freliste)),dtype = bool)
    RFIpeaks[bpb1:bpe2] = (datadiga/bandpass[bpb1:bpe2]<meanbp*1.05)

    return RFIpeaks

datadig = data.diagonal()
RFIpeaks = RFI_peaks(freliste,datadig)

frelisrfi = freliste[RFIpeaks]
datadigrfi = datadig[RFIpeaks]
print frelisrfi
print datadigrfi

RFIpeaks = RFI_peaks(frelisrfi,datadigrfi)

print len(RFIpeaks)

datadig = data.diagonal()
datadiga = (datadig[bpb1:bpe2].mean(axis=1))

semilogy(freliste[bpb1:bpe2],datadiga/bandpass[bpb1:bpe2],'b')
semilogy(freliste[bpb1:bpe2],1.01*meanbp *ones((len(bandpass[bpb1:bpe2]))),'r')
ad = arange(bpb1,bpe2,1)
semilogy(ad[RFIpeaks]*100.0/512,(datadiga/bandpass[bpb1:bpe2])[RFIpeaks],'k+')

#semilogy(freliste[bpb1:bpe2],datadiga,'b')
#semilogy(freliste[bpb1:bpe2],bandpass[bpb1:bpe2],'k')
#print bandpass[170]
#print bandpass[370]
show()
exit()
"""


#################################################################
###----------              CALIBRATION              ----------###
#################################################################

#-------------Calibration using Autocorrelations-----------------

def gain_from_autocorr(datas_c,gain,rcuflag):
    print "\n- Start Calculating Gains using Autocorrelations -"

    datas_c2 = datas_c.diagonal()
    datas_c2[rcuflag] = 0.0
    naverX = median(datas_c2[::2].real)
    naverY = median(datas_c2[1::2].real)

    gainn = zeros((nrcu),dtype = complex)
    gainn[::2] = sqrt(naverX/datas_c2[::2].real)
    gainn[1::2] = sqrt(naverY/datas_c2[1::2].real)
    gainn[rcuflag] = 0.0

    return gainn

#-------------Phase Calibration----------------------------------

# based on AIPS GCALC for phase calibration, but modified

def CAL_PHASE(diff):

    print "\n- Start Phase Calibration -"
    Gstart = ones((nrcu), dtype= complex)

    for i in range(0,20,1):     # number of iterations

        z = array(Gstart.conjugate())*array(diff)
        gnx = sum(z[::2,::2],axis = 1) + sum(z[::2,::2].conjugate(),axis = 0)
        gny = sum(z[1::2,1::2],axis = 1) + sum(z[1::2,1::2].conjugate(),axis = 0)

        gn = (hstack((gnx.reshape(nrcu/2,1),gny.reshape(nrcu/2,1)))).reshape(nrcu)

        G_ph = angle(1/(gn))
        G_ph[gn==0.0] = 0.0
        Gstart = cos(G_ph)+1j*sin(G_ph)

    return Gstart

#-------------Phase + Amplitude Calibration----------------------

# based on AIPS GCALC1 Amp+Phase Calibration, but modified

def CAL_AMP_PHASE(diff):
    print "\n- Start Amp-Phase Calibration -\n"

    flagom = ones((len(diff),len(diff[0])))
    flagom[diff==0.0] = 0.0
    Gstart = ones((nrcu), dtype= complex)
    eps = [.5e-3,.5e-4,.5e-5]

    for j in range(0,3,1):
        epsi = eps[j]
        for i in range(0,30,1):     # number of iterations

            Gstartm = dot(mat(Gstart).transpose(),mat(Gstart).conjugate())
            Gdiag = (Gstartm.diagonal()).transpose()
            t = array(1.0/sqrt((abs(diff-Gstartm))**2+epsi))
            t *= flagom

            Gstartx = array([array(Gstart)]*nrcu)
            Gstarty = Gstartx.transpose()

            gnx = array(t)*array(diff)*array(Gstarty)
            gny = array(t)*array(diff.conjugate())*array(Gstartx)

            gncx = sum(gnx[::2,::2],axis=0) + sum(gny[::2,::2],axis=1)
            gncy = sum(gnx[1::2,1::2],axis=0) + sum(gny[1::2,1::2],axis=1)
            gnc = (hstack((gncx.reshape(nrcu/2,1),gncy.reshape(nrcu/2,1)))).reshape(nrcu)

            QQix = (abs(array(Gstartx)))**2*array(t)
            QQiy = (abs(array(Gstarty)))**2*array(t)
            QQx = sum(QQix[::2,::2],axis=0) + sum(QQiy[::2,::2],axis=1)
            QQy = sum(QQix[1::2,1::2],axis=0) + sum(QQiy[1::2,1::2],axis=1)
            QQ = (hstack((QQx.reshape(nrcu/2,1),QQy.reshape(nrcu/2,1)))).reshape(nrcu)

            Gstart += 0.25*(gnc/(QQ) - Gstart)
            Gstart[QQ == 0.0] = 0.0

            Gstart[rcuflag] = 0.0

    return Gstart


#-------------Apply Gains to the data----------------------------

def Gain(Gain,dat):

    dat /= dot(Gain.conjugate().reshape(nrcu,1),(Gain).reshape(1,nrcu))

    dat[rcuflag] = 0.0
    dat[:,rcuflag] = 0.0
    dat[Gain==0.0] = 0.0
    dat[:,Gain==0.0] = 0.0

    return dat


#################################################################
###----------           APPLYING A MODEL            ----------###
#################################################################

#-------------Calculating model visibilities---------------------

# calculates the model Visibilities for a given Source list

def vis_modelcxy(sourcelistxy, pol):

    vis_model = zeros((nrcu,nrcu),dtype = complex)

    for source in range(0, len(sourcelistxy), 1):
        l = int(sourcelistxy[source][2])
        m = int(sourcelistxy[source][3])

        ampX = sqrt(sqrt(sourcelistxy[source][0]))
        ampY = sqrt(sqrt(sourcelistxy[source][1]))

        wxx = wx[0,:,l].reshape(nrcu/2,1)
        wxn = (hstack((ampX*wxx,ampY*wxx))).reshape(nrcu,1)
        wyx = wy[0,:,m].reshape(nrcu/2,1)
        wyn = (hstack((ampX*wyx,ampY*wyx))).reshape(nrcu,1)

        wet = (wxn.conjugate()*wyn)
        wet2 = (wet*wet.conjugate().transpose())
        vis_model += wet2

    return vis_model

def vis_modelc0(sourcelistx,sourcelisty, pol):

    vis_model = zeros((nrcu,nrcu),dtype = complex)

    for source in range(0, min(len(sourcelistx),len(sourcelisty)), 1):
        lx = int(sourcelistx[source][1])
        mx = int(sourcelistx[source][2])
        ly = int(sourcelisty[source][1])
        my = int(sourcelisty[source][2])
        ampX = sqrt(sqrt(sourcelistx[source][0]))
        ampY = sqrt(sqrt(sourcelisty[source][0]))

        wxx = wx[0,:,lx].reshape(nrcu/2,1)
        wxy = wx[0,:,ly].reshape(nrcu/2,1)
        wxn = (hstack((ampX*wxx,ampY*wxy))).reshape(nrcu,1)
        wyx = wy[0,:,mx].reshape(nrcu/2,1)
        wyy = wy[0,:,my].reshape(nrcu/2,1)
        wyn = (hstack((ampX*wyx,ampY*wyy))).reshape(nrcu,1)

        wet = (wxn.conjugate()*wyn)
        wet2 = (wet*wet.conjugate().transpose())
        vis_model += wet2

    return vis_model

#-------------Difference between model and data------------------

def diff_model_real(dat,vis_model):

    diff = zeros((nrcu,nrcu),dtype = complex)

    diff = dat/vis_model
    diff[baselines <= max_baselenght] = 0.
    diff[vis_model == 0.] = 0.

    return diff


#################################################################
###----------            SOURCE FINDING             ----------###
#################################################################

beam = 1.22* (3*10**2/freliste)/40/0.017
beam_deg = beam/90

#-------------Search for Sources in Skymap-----------------------
# searchs for a given number of strong sources in a given Skymap

def RFI_search(skymap,num):

    RFIso = []

    skymap = skymap[skymap[:,0]>0.0]    # remove values < 0
    skymaps = sorted(skymap, key=lambda mapi:mapi[0], reverse=True)

    RFIso.append(skymaps[0])

    for peak in range(1, len(skymaps), 1):
        inBeam = False
        for so in range(0, len(RFIso), 1):
            if (float(RFIso[so][1])-beam[fre] < float(skymaps[peak][1]) < float(RFIso[so][1])+beam[fre]) and (float(RFIso[so][2])-beam[fre] < float(skymaps[peak][2]) < float(RFIso[so][2])+beam[fre]):
                inBeam = True
            if inBeam:
                break
        if not inBeam:
            RFIso.append(skymaps[peak])
        if len(RFIso) >= num:
            break

    return RFIso

#----------------------------------------------------------------

def source_XY(sourcelisteX, sourcelisteY, skymapx, skymapy):

    sourcelisteXY = []

    if len(sourcelisteX) != len(sourcelisteY):
        minlen = min(len(sourcelisteX),len(sourcelisteY))
        sourcelisteX = sourcelisteX[:minlen]
        sourcelisteY = sourcelisteY[:minlen]

    ara = arange(0,len(sourcelisteX),1)
    testp = zeros((len(sourcelisteX)))


    for so in range(0, len(sourcelisteX),1):
        testxm = (abs(sourcelisteX[so,1]-sourcelisteY[:,1]) <= beam[fre]/2) & (abs(sourcelisteX[so,2]-sourcelisteY[:,2]) <= beam[fre]/2)
        testp += (testxm)*(so+1)


    for so in range(0,len(testp),1):
        if testp[so] != 0:
            sourcelisteXY.append((sourcelisteX[testp[so]-1,0], sourcelisteY[so,0], sourcelisteY[so,1],sourcelisteY[so,2]))
            ara[testp[so]-1] = 0
        else:
            sourcelisteXY.append((skymapx[int(sourcelisteY[so,2])*2*resol+int(sourcelisteY[so,1])][0] , sourcelisteY[so,0], sourcelisteY[so,1],sourcelisteY[so,2]))

    for so in range(0,len(ara),1):
        if ara[so] != 0:
            sourcelisteXY.append((sourcelisteX[so,0], skymapy[int(sourcelisteX[so,2])*2*resol+int(sourcelisteX[so,1])][0] , sourcelisteX[so,1],sourcelisteX[so,2]))

    return array(sourcelisteXY)


#-------------New Source Positions for Ateam---------------------

# finds the new positions of the Ateam sources in the Skymap

def new_pos(Ateam3,sourceliste3):

    Ateamn = Ateam3*1

    for so in range(0, len(sourceliste3), 1):
        for sou in range(0, len(Ateam3), 1):
            if (float(Ateam3[sou][1])-beam[fre]/2 < float(sourceliste3[so][1]) < float(Ateam3[sou][1])+beam[fre]/2) and (float(Ateam3[sou][2])-beam[fre]/2 < float(sourceliste3[so][2]) < float(Ateam3[sou][2])+beam[fre]/2):
                Ateamn[sou] = sourceliste3[so]

    return Ateamn

def new_pos_so(Ateam3,sourceliste3):

    Ateamn = Ateam3*1

    for so in range(0, len(sourceliste3), 1):
        for sou in range(0, len(Ateam3), 1):
            if (float(Ateam3[sou][1])-beam[fre]/2 < float(sourceliste3[so][1]) < float(Ateam3[sou][1])+beam[fre]/2) and (float(Ateam3[sou][2])-beam[fre]/2 < float(sourceliste3[so][2]) < float(Ateam3[sou][2])+beam[fre]/2):
                sourceliste3[so][1] = Ateam3[sou][1]
                sourceliste3[so][2] = Ateam3[sou][2]

    return sourceliste3


#################################################################
###----------                SKYMAPS                ----------###
#################################################################


# values of the direction cosines
resol = 100     # how much pixels for the skymaps

if output == 'galactic':
    ll = -1.0*arange(-1.0, 1.0, 1.0/resol)
elif rcumode == 5:
    ll = 1.0*arange(-1.0, 1.0, 1.0/resol)
else:
    ll = -1.0*arange(-1.0, 1.0, 1.0/resol)

mm = -1*ll
lenght_l = len(ll)
lenght_m = len(mm)

r = sqrt((array([ll]*lenght_l))**2 + (array([mm]*lenght_m).transpose())**2)

# boundaries for the Skymaps
rand = ones((2*resol,2*resol))
randz = ones((2*resol,2*resol))
rande = ones((2*resol,2*resol))
rand[r>1] = 'NaN'
randz[r>1] = 0.
rande[r>0.9] = 'NaN'

wx0 = [exp((dot((mat(xpos)).T,(mat(ll)))))]
wx0 = asarray(wx0)
wy0 = [exp((dot((mat(ypos)).T,(mat(mm)))))]
wy0 = asarray(wy0)


#-------------Calculating the skymap-----------------------------

# calculates the Skymap for given data and polarisation

# for the imaging part without beam:
def SKYMAP(dat, pol):

    if pol == 'XX': acm = dat[::2,::2]          # just X dipoles
    if pol == 'YY': acm = dat[1::2,1::2]            # just Y dipoles
    if pol == 'I': acm = dat[::2,::2] + dat[1::2,1::2]  # Stokes I - total intensity
    if pol == 'Q': acm = dat[1::2,1::2] - dat[::2,::2]  # Stokes Q - linear polarisation
    if pol == 'U': acm = dat[1::2,::2] + dat[::2,1::2]  # Stokes U - orthogonal lin. pol.
    if pol == 'V': acm = dat[::2,1::2] - dat[1::2,::2]  # Stokes V - circular polarisation
    if pol == 'XY': acm = dat[::2,1::2]         # XY Correlations
    if pol == 'YX': acm = dat[1::2,::2]         # YX Correlations

    prod1 = inner(wei.conjugate(),acm)
    prod1 *= wei

    skymap = (sum(prod1,axis =1))
    skymap = rot90(skymap.reshape(2*resol,2*resol))

    return skymap*rand      # is a complex number


# for the imaging part with beam model:
def SKYMAP_beam(dat, pol):

    jonesx = ((abs(Jones[fre,0]))**2 + (abs(Jones[fre,1]))**2).real
    jonesy = ((abs(Jones[fre,2]))**2 + (abs(Jones[fre,3]))**2).real
    jonesxy = (abs(Jones[fre,0]*Jones[fre,2].conjugate()+Jones[fre,1]*Jones[fre,3].conjugate()) + abs(Jones[fre,0]*Jones[fre,3].conjugate() - Jones[fre,1]*Jones[fre,2].conjugate())).real

#   jonesx = beamx*1
#   jonesy = beamy*1

    b = 1
    if pol == 'XX':                 # just X dipoles
        acm = dat[::2,::2]
        beam = jonesx*1
    if pol == 'YY':                 # just Y dipoles
        acm = dat[1::2,1::2]
        beam = jonesy*1
    if pol == 'I':                  # Stokes I - total intensity
        acm1 = dat[::2,::2]
        acm2 = dat[1::2,1::2]
        beam1 = jonesx*1
        beam2 = jonesy*1
        b = 2
    if pol == 'Q':                  # Stokes Q - linear polarisation
        acm2 = dat[1::2,1::2]
        acm1 = dat[::2,::2]
        beam1 = jonesx*1
        beam2 = jonesy*1
        b = 2
    if pol == 'U':                  # Stokes U - orthogonal lin. pol.
        acm = dat[1::2,::2] + dat[::2,1::2]
        beam = jonesxy*1
    if pol == 'V':                  # Stokes V - circular polarisation
        acm = dat[::2,1::2] - dat[1::2,::2]
        beam = jonesxy*1
    if pol == 'XY':
        acm = dat[::2,1::2]         # XY Correlations
        beam = jonesxy*1
    if pol == 'YX':
        acm = dat[1::2,::2]         # YX Correlations
        beam = jonesxy*1

    if b == 1:
        prod1 = inner(wei.conjugate(),acm)
        prod1 *= wei
        skymap = (sum(prod1,axis =1))
        skymap = rot90(skymap.reshape(2*resol,2*resol))
        skymap /= beam
    if b == 2:
        prod1 = inner(wei.conjugate(),acm1)
        prod1 *= wei
        skymap = (sum(prod1,axis =1))
        skymap = rot90(skymap.reshape(2*resol,2*resol))
        skymap /= beam1

        prod1 = inner(wei.conjugate(),acm2)
        prod1 *= wei
        skymap2 = (sum(prod1,axis =1))
        skymap2 = rot90(skymap2.reshape(2*resol,2*resol))
        skymap2 /= beam2

        if pol == 'I':
            skymap += skymap2
        if pol == 'Q':
            skymap -= skymap2

    return skymap*rand      # is a complex number


matl = array([arange(0,2*resol,1)]*2*resol)
matm = (matl.transpose()).flatten()
matl = matl.flatten()

# for the source search part:
def SKYMAP2(dat, pol):

    if pol == 0: acm = dat[::2,::2]
    if pol == 1: acm = dat[1::2,1::2]
    if pol == 2: acm = dat[::2,::2] + dat[1::2,1::2]

    prod1 = inner(wei.conjugate(),acm)
    prod1 *= wei

    skymap0 = (sum(prod1, axis = 1)).real
    skymap0 = (rot90(skymap0.reshape(2*resol,2*resol)))
    skymap = (vstack((skymap0.reshape(4*resol**2)*rand.flatten(),matl,matm))).transpose()

    return skymap           # is a real number + position (l,m)


#################################################################
###--------- TRANSFORMATION TO GALACTIC COORDINATES ----------###
#################################################################

def trafo_galactic_coordinates():

    MM, LL = meshgrid(mm,ll)
    A = 90 - arctan(-MM/LL)*180/pi
    A [ LL>0 ] += 180
    Elev = arccos(MM/sin((90.-A)*pi/180))*180./pi
    Elev [sin((90.0-A)*pi/180) == 0] = 0.
    Azimut = mod(A,360)
    Elev[MM==0] = 0.0#Elev[MM-1]
    Azimut[:,LL==0] = 0.
    Elev[:,LL==0] = 0.0#Elev[:,LL-1]
    Azimut *= randz
    Elev *= randz
    Azimut = (Azimut).real
    Elev = Elev.real

    L = mod(280.460 + 0.9856474*nJD, 360)     #mean ecliptic longitude of the sun; aberration is included
    g = mod(357.528 + 0.9856003*nJD, 360)   #mean anormaly
    delta = mod((L + 1.915*math.sin(pi*g/180)+0.02*math.sin(2*pi*g/180)), 360)   #ecliptic longitude of the sun
    ecliptic = (23.439 - 0.0000004*nJD)*pi/180

    decl = arcsin(-cos(Elev*pi/180)*cos(Azimut*pi/180)*cos(latitude*pi/180)+sin(Elev*pi/180)*sin(latitude*pi/180))*180/pi
    ta = cos(Azimut*pi/180)*sin(latitude*pi/180)+tan(Elev*pi/180)*cos(latitude*pi/180)
    t = arctan(sin(Azimut*pi/180)/ta)*180/pi
    t[ta==0] = 0.
    t[ta<0] += 180
    rect = mod(vernalequinox - t,360)

    gb = arcsin(sin(decl*pi/180)*sin(27.4*pi/180)+cos(decl*pi/180)*cos(27.4*pi/180)*cos((192.25-rect)*pi/180))*180/pi + 90.
    g1 = cos((192.25-rect)*pi/180)*sin(27.4*pi/180)-tan(decl*pi/180)*cos(27.4*pi/180)
    gl = 303-arctan(sin((192.25-rect)*pi/180)/g1)*180/pi
    gl[g1>0] += 180
    gl = mod(gl,360)

    return gl, gb


#################################################################
###----------            IMAGING ROUTINE            ----------###
#################################################################

def do_imaging():
    print "- Start Imaging -"

    if output != 'galactic':
        Stok_len = float(len(map_Stokes))
        s1 = round(sqrt(Stok_len))
        s2 = round(Stok_len/s1)
        if plotsources:
            Ateamplot = []
            for so in sourcelistplot:
                l,m,I,Azi,Elev = source_l_m(so)
                if l != 'NaN':
                    Ateamplot.append((so,l,m))

        for stok in range(0,int(Stok_len),1):
            if map_Stokes[stok] != 5:
                if options.beam:
                    skymap = SKYMAP_beam(datas,map_Stokes[stok]).real
                else:
                    skymap = SKYMAP(datas,map_Stokes[stok]).real
            else:
                if options.beam:
                    skymap = SKYMAP_beam(datas,map_Stokes[stok]).imag
                else:
                    skymap = SKYMAP(datas,map_Stokes[stok]).imag

            if savedata:
                ofile = 'datamap_'+str(map_Stokes[stok])+'_'+str(date[0:15])+'_ch_'+str(int(fre))+'.out'
                of = open(ofile,'w')
                pickle.dump(skymap,of)
                print "\nSkymap data sucessfully written to file: ", ofile
            else:
                subplot(int(s1), int(s2), int(stok+1))
                imshow((skymap), aspect='auto', interpolation='bilinear',extent=(1,-1,-1,1))
                title(str(map_Stokes[stok])+'     '+str(round(f/10**6,2))+' MHz',size=12)
                xlabel('E'+' '*15+'l'+' '*15+'W', size=10)
                ylabel('S'+' '*15+'m'+' '*15+'N', size=10)
                colorbar()

                if plotsources:
                    for so in range(0,len(Ateamplot),1):
                        plot(Ateamplot[so][1], Ateamplot[so][2],'k+')
                        text(Ateamplot[so][1]-0.05, Ateamplot[so][2]-0.05, str(Ateamplot[so][0]))

    if output == 'galactic':

        datad = zeros((nrcu,nrcu), dtype = complex)
        datad[eye(nrcu)==1] = datas[eye(nrcu)==1]*1
        datas[eye(nrcu)==1] = 0.0

        Stokes = ['I']
        stok = 0
        if Stokes[stok] != 5:
            if options.beam:
                skymap = SKYMAP_beam(datas,Stokes[stok]).real
            else:
                skymap = SKYMAP(datas,Stokes[stok]).real
        else:
            if options.beam:
                skymap = SKYMAP_beam(datas,Stokes[stok]).imag
            else:
                skymap = SKYMAP(datas,Stokes[stok]).imag


        r = sqrt((array([ll]*lenght_l))**2 + (array([mm]*lenght_m).transpose())**2)
    #   f = freliste[channel] * 1000000
        c = 2.9979245e8
        wl = c/f
        D = 60.
        n = r*D/wl*pi

        skymapd = (SKYMAP(datad,Stokes[stok]).real)*(sin(n*pi/180.)/n)**2

        skymap += skymapd


        skymap[isnan(skymap)] = 0.0

        if (Stokes[stok] == 'I') and (options.calibrate == True):
            quelle = 0
            ampcas = skymap[int(Ateamny[quelle][2]),int(Ateamny[quelle][1])]
            n = Ateam2[quelle][0]/ampcas
            skymap *= n

        gal_lon, gal_lat = trafo_galactic_coordinates()

        mapn = zeros((180+1,360+1))

        for l in range(0,len(ll),1):
            for m in range(0,len(mm),1):
                d = round(gal_lat[l,m],1)
                r = round(gal_lon[l,m],1)
                mapn[d,r] = skymap[m,l]

        mapn[isinf(mapn)] = 0.0
        mapn[isnan(mapn)] = 0.0

        mapn = ma.masked_where(mapn==0,mapn)
        lon = arange(-180.,180.,1)*(-pi/180)
        lat = arange(-90.,90.,1)*(pi/180)
        x,y = meshgrid(lon,lat)

        fig = figure(figsize=(10,5))
        fig.patch.set_facecolor('w')
        ax = fig.add_subplot(111,projection = 'mollweide')
        cmap=cm.jet
        cmap.set_bad('k',1)

        image = ax.pcolormesh(x,y,mapn,interpolation='bicubic',cmap=cmap)

        cb = fig.colorbar(image)
        axis('off')

    #   figure(figsize=(10,5))
    #   axes(projection = 'mollweide')
    #   imshow(mapn, aspect='auto',interpolation='bicubic',extent=(pi,-pi,pi/2,-pi/2))
    #   colorbar()

        if savedata:
            ofile = str(directory)+'/datagalactic_'+str(Stokes[stok])+'_'+str(date[0:15])+'_ch_'+str(int(fre))+'.out'
            of = open(ofile,'w')
            pickle.dump(mapn,of)
            print "\nGalactic data sucessfully written to file: ", ofile

    if saveimages:
        savefig(str(directory)+'/map_'+str(date[0:15])+'_ch_'+str(int(fre))+'.png')
    #   exit()
    elif showimages:
        show()
#        exit()




#################################################################
###----------          CALIBRATION ROUTINE          ----------###
#################################################################


if output == 'spectra':
    figure()
    for rcu in range(0,nrcu,2):

        semilogy(freliste,abs(data[rcu,rcu],),'r')
        semilogy(freliste,abs(data[rcu+1,rcu+1]),'b')
    legend(('X-Dipols','Y-Dipols'))
    xlabel('Frequency [MHz]')
    title('Autocorrelation Spectra of all RCUs')
    if saveimages == True:
        savefig(str(directory)+'/spectra_'+str(date[0:15])+'.png')
    else:
        show()
    exit()

if output == 'gain':
    Gainf = zeros((len(freliste),nrcu), dtype = complex)

if options.calibrate:
    print "\nCALIBRATION STRATEGY:"
    print "\nPhase     ", options.phasecal
    print "Amp+Phase ", options.ampphacal, "\n"
else:
    print "\nNo Calibration\n"
print "Beam      ", options.beam


if Calall:  # Calibration over all frequencies
    start_chan = 0
    end_chan = len(freliste)

if options.integratefrequencies:
    datasall = zeros((nrcu,nrcu), dtype = complex)

if not options.calibrate:
    for fre in range(start_chan, end_chan, 1):
        f = freliste[fre] * 10**6
        c = 2.9979245e8
        k = 2.0*pi / c  * f

        print "\nCHANNEL: %i  -  FREQUENCY: %.2f MHz" % (fre, f/10**6)

        wx = wx0**(-(1j)*k)
        wy = wy0**(-(1j)*k)
        wei = (array(wx).transpose(1,2,0)*array(wy).transpose(1,0,2)).reshape(nrcu/2,4*resol**2).transpose()

        datas = data[:,:,fre]*1
        datas[eye(nrcu)==1] = 0.0


        if options.blflag == 's':
            datfl = datas.flatten()
            flagbl = baselines.flatten() < 40
            flagbl = invert(flagbl)
            datfl[flagbl] = 0.0
            datas = datfl.reshape(nrcu,nrcu)
        if options.blflag == 'l':
            datfl = datas.flatten()
            flagbl = baselines.flatten() > 40
            flagbl = invert(flagbl)
            datfl[flagbl] = 0.0
            datas = datfl.reshape(nrcu,nrcu)

        if options.integratefrequencies:
            datasall += datas

    if options.integratefrequencies:
        datas = datasall/(end_chan-start_chan)

    if (showimages != False) or (saveimages != False) or (savedata != False):
        do_imaging()
    exit()


#-------------doing the Calibration------------------------------

if options.calibrate:

    data, rcuflag = flag_from_autocorr(data)

    for fre in range(start_chan, end_chan, 1):

        f = freliste[fre] * 10**6
        c = 2.9979245e8
        k = 2.0*pi / c  * f

        print "\nCHANNEL: %i  -  FREQUENCY: %.2f MHz" % (fre, f/10**6)

        wx = wx0**(-(1j)*k)
        wy = wy0**(-(1j)*k)
        wei = (array(wx).transpose(1,2,0)*array(wy).transpose(1,0,2)).reshape(nrcu/2,4*resol**2).transpose()

        Ateam2 = []
        for so in sou:      # get right source intensity for the given frequency
            l,m,I,Azi,Elev = source_l_m(so)
            if l != 'NaN':
                lp = int(resol-l*resol)
                mp = int(resol-m*resol)
                Ateam2.append((I,lp,mp))

        datas = data[:,:,fre]*1     # choose one channel for the calibration
        gain = ones((nrcu), dtype = complex)
        Gaing = gain+0.0*1j

        if autoc:
            gain = gain_from_autocorr(datas, gain, rcuflag)
            Gaing = gain+0.0*1j
            datas = Gain(Gaing,datas)

    #   datas[eye(nrcu)==1] = 0.0   # flagging the autocorrelations
    #   data[:,:,fre][eye(nrcu)==1] = 0.0

        datas = FLAGGING(datas)

        if options.blflag == 's':
            datfl = datas.flatten()
            flagbl = baselines.flatten() < 40
            flagbl = invert(flagbl)
            datfl[flagbl] = 0.0
            datas = datfl.reshape(nrcu,nrcu)
        if options.blflag == 'l':
            datfl = datas.flatten()
            flagbl = baselines.flatten() > 40
            flagbl = invert(flagbl)
            datfl[flagbl] = 0.0
            datas = datfl.reshape(nrcu,nrcu)

    #   plot(baselines,datas,'ko')
    #   show()
    #   exit()

    #-------------doing the Phase calibration------------------------

        if options.phasecal:

            if options.find_sources:
                skymapx = SKYMAP2(datas,0)
                skymapy = SKYMAP2(datas,1)
                sourcelisteX = RFI_search(skymapx,10)
                sourcelisteY = RFI_search(skymapy,10)

            #   sourcelisteX = new_pos_so(Ateam2,sourcelisteX)
            #   sourcelisteY = new_pos_so(Ateam2,sourcelisteY)

                vismodel = vis_modelc0(sourcelisteX,sourcelisteY,pol)
            else:
                vismodel = vis_modelc0(Ateam2,Ateam2, pol)

            diff = diff_model_real(datas,vismodel)

            # apply some weighting
            RMS = sqrt((abs(diff[::2,::2])**2).mean())
            RMSY = sqrt((abs(diff[1::2,1::2])**2).mean())

            diff[::2,::2] *= RMS/(abs(diff[::2,::2]) +RMS*0.05 )
            diff[1::2,1::2] *= RMSY/(abs(diff[1::2,1::2]) +RMSY*0.05 )
            diff[rcuflag] = 0.0
            diff[:,rcuflag] = 0.0

            G_phase = CAL_PHASE(diff)
            datas = Gain(G_phase,datas)
            Gaing *= G_phase

            if flagging:
                datas = FLAGGING(datas)

    #-------------doing the Amplitude-Phase Calibration--------------

        if options.ampphacal:

        #   if options.find_sources:
            skymapx = SKYMAP2(datas,0)
            skymapy = SKYMAP2(datas,1)
            sourcelisteX = RFI_search(skymapx,10)
            sourcelisteY = RFI_search(skymapy,10)

            Ateamnx = new_pos(Ateam2,sourcelisteX)
            Ateamny = new_pos(Ateam2,sourcelisteY)

            sourcelisteX = array(sourcelisteX)
            sourcelisteY = array(sourcelisteY)
            sourcelisteXY = source_XY(sourcelisteX, sourcelisteY, skymapx, skymapy)

            """
            quelle = 0
            ampcas = (skymapx[int((Ateamnx[quelle][2])*200+Ateamnx[quelle][1])][0]+skymapy[int((Ateamny[quelle][2])*200+Ateamny[quelle][1])][0])
            flagi = len(datas[datas[::2,::2]==0.0])+len(datas[datas[1::2,1::2]==0.0])
            n = Ateam2[quelle][0]/ampcas/((nrcu/2)**2-flagi/2)

            sourcelisteXY[:,0] *= n/4.
            sourcelisteXY[:,1] *= n/4.
            """

            vismodel = vis_modelcxy(sourcelisteXY,pol)

            diff = diff_model_real(datas,vismodel)

            G_amp_ph = CAL_AMP_PHASE(diff)

            datas = Gain(G_amp_ph,datas)
            Gaing *= G_amp_ph

            if flagging:
                datas = FLAGGING(datas)

        if output == 'gain':
            Gainf[fre] = Gaing

        data[:,:,fre] = datas

        if options.integratefrequencies:
            datasall += datas
        elif (showimages != False) or (saveimages != False) or (savedata != False):
            do_imaging()


if options.integratefrequencies:
    datas = datasall/(end_chan-start_chan)

    if (showimages != False) or (saveimages != False) or (savedata != False):
        do_imaging()

if output == 'gain':
    of = open(outfile,'w')
    pickle.dump(Gainf,of)
    print "\nGains sucessfully written to file: ", outfile
    exit()

#-------------save image as fits file----------------------------
# not really tested by now

if output == 'mapf':
    skymap = reshape((skymap),(1,1,lenght_l,lenght_m))

    hdu = pyfits.PrimaryHDU(skymap)

    hdulist = pyfits.HDUList([hdu])
    prihdr = hdulist[0].header

#   prihdr.update('NAXIS',3)
    prihdr.update('NAXIS1',lenght_l)
    prihdr.update('NAXIS2',lenght_m)
    prihdr.update('NAXIS3',1)
    prihdr.update('NAXIS4',1)
    prihdr.update('EQUINOX',2000.0)
    #prihdr.update('LONPOLE',0.)
    #prihdr.update('LATPOLE',90.)
    prihdr.update('OBSLON', longitude*0.017)
    prihdr.update('OBSLAT', latitude*0.017)
#   prihdr.update('OBSTIME', timetuple)
    prihdr.update('CTYPE1','RA---TAN')
    prihdr.update('CRVAL1',GWhourangle)  #value of the axis at the reference point "CRPIX"
    prihdr.update('CDELT1',-(180.0/lenght_l))
    prihdr.update('CRPIX1',float(lenght_l/2))
    #prihdr.update('CROTA1',0.)          #rotation, just leave it at 0.
    prihdr.update('CUNIT1','deg')       #the unit in which "CRVAL" and "CDELT" are given
    prihdr.update('CTYPE2','DEC--TAN')
    prihdr.update('CRVAL2',latitude)
    prihdr.update('CDELT2',180.0/lenght_l)
    #prihdr.update('CROTA2',0.)
    prihdr.update('CRPIX2',float(lenght_l)/2)
    prihdr.update('CUNIT2','deg     ')
    prihdr.update('CTYPE3','FREQ    ')
    prihdr.update('CRVAL3',f)
    prihdr.update('CDELT3',frequency)
    prihdr.update('CUNIT3','HZ      ')

    prihdr.update('CUNIT4','        ')

    if os.path.isfile(outfile):
        os.remove(outfile)

    hdulist.writeto(outfile)

    exit()

#exit()






