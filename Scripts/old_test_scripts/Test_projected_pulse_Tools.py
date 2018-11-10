'''
This script is used to characterize pulses (only works for L43784_D20120125T211154)
NOTE : The addition of stations is done incoherently!!!
There are several options to use this script:
    - add dispersed/dedispersed station beams
    - add clean/unclean station beams
    - add the two polarizations
    - plots a 3 window plot with the dynamic spectrum plus two plots with a "collapsed" axis.
'''

import pycrtools as cr
import numpy as np
import sys; import os
import pdb;# pdb.set_trace()
import Beam_Tools as bt


def super_plot(zoom_dynspec,pulse=0,zoom=1):
    '''Do some awesome ploting.
    
    '''

    #Create a frequency vector
    frequencies = zoom_dynspec.par.yvalues

    #Create a time vector
    times = zoom_dynspec.par.xvalues

#    cr.plt.subplot(1,len(filenames),i)
    cr.plt.clf()
    
    zoom=0
    if zoom:
        cr.plt.rcParams['font.size']=20
        cr.plt.rcParams['axes.titlesize']=30
        flat_dyn=bt.flatdyn(zoom_dynspec,pulse=pulse)
        cr.plt.subplot(2,2,1)
        cr.plt.title("Pulse Profile from : " + zoom_dynspec.par.station)
        cr.plt.plot(flat_dyn.par.xvalues,flat_dyn)
        cr.plt.xlabel("+/- Time [s]")
        cr.plt.autoscale(axis='x',tight=True)
        pdb.set_trace()
    
        if pulse!=3:  # Need to fix flatdyn for any kind of pulse ... or need to actually remove all the hardcoded stuff.
            flat_dyn,flat_dyn2=bt.flatdyn(zoom_dynspec,axis='y',pulse=pulse)
            cr.plt.subplot(2,2,4)
            cr.plt.plot(flat_dyn.par.xvalues,flat_dyn)
            cr.plt.plot(flat_dyn.par.xvalues,flat_dyn2,'r')
            cr.plt.autoscale(axis='x',tight=True)
            cr.plt.xlabel("Frequency [MHz]")
    
        cr.plt.subplot(2,2,3)
    #Time integration.
#    tbin=9
    if pulse:
        tbin=4
    else:
        tbin=6
    if not zoom:
        tbin=1
    if tbin>1:
        zoom_dynspec = cr.hArray_toNumpy(zoom_dynspec)
        zoom_dynspec = zoom_dynspec.reshape((zoom_dynspec.shape[0]/tbin,tbin,zoom_dynspec.shape[1]))
        zoom_dynspec = np.sum(zoom_dynspec,axis=1)
        zoom_dynspec = cr.hArray(zoom_dynspec)
#    tbin=8
    tbin=1
    zoom_dynspec = zoom_dynspec.Transpose()
    if tbin>1:
        zoom_dynspec = cr.hArray_toNumpy(zoom_dynspec)
        zoom_dynspec = zoom_dynspec.reshape((zoom_dynspec.shape[0]/tbin,tbin,zoom_dynspec.shape[1]))
        zoom_dynspec = np.sum(zoom_dynspec,axis=1)
        zoom_dynspec = cr.hArray(zoom_dynspec)
    zoom_dynspec = zoom_dynspec.Transpose()
    cr.plt.imshow(np.log10(cr.hArray_toNumpy(zoom_dynspec)),aspect='auto',origin='lower',cmap=cr.plt.cm.hot,extent=(times[0],times[-1],frequencies[0],frequencies[-1]))
#    pdb.set_trace()
#    cr.plt.imshow(np.log10(cr.hArray_toNumpy(zoom_dynspec)),aspect='auto',origin='lower',cmap=cr.plt.cm.hot,vmin=-4,vmax=-3,extent=(times[0],times[-1],frequencies[0],frequencies[-1]))
    cr.plt.xlabel("+/- Time [s]")
    cr.plt.ylabel("Frequency [MHz]")
    
#------------------------------------------------------------------------------------------------------------------------
#-------------------------------PREPARE FOR SUPER SCRIPT BELOW-----------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------

#filenames = ['L43784_D20120125T211154.887Z_CS002_R000_tbb.pol0.multi_station.beam','L43784_D20120125T211154.867Z_CS003_R000_tbb.pol0.multi_station.beam','L43784_D20120125T211154.887Z_CS004_R000_tbb.pol0.multi_station.beam','L43784_D20120125T211154.866Z_CS005_R000_tbb.pol0.multi_station.beam','L43784_D20120125T211154.871Z_CS006_R000_tbb.pol0.multi_station.beam','L43784_D20120125T211154.887Z_CS007_R000_tbb.pol0.multi_station.beam','L43784_D20120125T211154.887Z_CS002_R000_tbb.pol1.multi_station.beam','L43784_D20120125T211154.867Z_CS003_R000_tbb.pol1.multi_station.beam','L43784_D20120125T211154.887Z_CS004_R000_tbb.pol1.multi_station.beam','L43784_D20120125T211154.866Z_CS005_R000_tbb.pol1.multi_station.beam','L43784_D20120125T211154.871Z_CS006_R000_tbb.pol1.multi_station.beam','L43784_D20120125T211154.887Z_CS007_R000_tbb.pol1.multi_station.beam']
#filenames = ['L43784_D20120125T211154.887Z_CS002_R000_tbb.pol0.multi_station.beam','L43784_D20120125T211154.867Z_CS003_R000_tbb.pol0.multi_station.beam','L43784_D20120125T211154.887Z_CS004_R000_tbb.pol0.multi_station.beam','L43784_D20120125T211154.866Z_CS005_R000_tbb.pol0.multi_station.beam','L43784_D20120125T211154.871Z_CS006_R000_tbb.pol0.multi_station.beam','L43784_D20120125T211154.887Z_CS007_R000_tbb.pol0.multi_station.beam']
filenames = ['L43784_D20120125T211154.887Z_CS002_R000_tbb.pol1.multi_station.beam','L43784_D20120125T211154.867Z_CS003_R000_tbb.pol1.multi_station.beam','L43784_D20120125T211154.887Z_CS004_R000_tbb.pol1.multi_station.beam','L43784_D20120125T211154.866Z_CS005_R000_tbb.pol1.multi_station.beam','L43784_D20120125T211154.871Z_CS006_R000_tbb.pol1.multi_station.beam','L43784_D20120125T211154.887Z_CS007_R000_tbb.pol1.multi_station.beam']
#filenames = ['L43784_D20120125T211154.887Z_CS002_R000_tbb.pol1.multi_station.beam','L43784_D20120125T211154.887Z_CS002_R000_tbb.pol0.multi_station.beam']
#filenames=filenames[1:]

#extra = 'dynspec.pcr/data.bin' 
extra = 'nd_clean_dynspec.pcr/data.bin' 
cr.plt.ion()

zoom = 0
verbose = 1
pulse = 0

full_blocks=12224

if zoom:
    if pulse==1: 
        nblocks=560
        freq_size=1147+1#1228
        incoherent_dynspec = cr.hArray(float,[freq_size, nblocks])
    else:
        nblocks=640
        freq_size=1638
        incoherent_dynspec = cr.hArray(float,[freq_size, nblocks])
else:
    tbin=16
    nblocks=12224
    incoherent_dynspec = cr.hArray(float,[ 8193, nblocks / tbin])

for i,file in enumerate(filenames):
    beams = cr.open(file)
    print 'Working on station: '+ beams['STATION_NAME'][0]
    full_dynspec=cr.hArray(float,(beams['BEAM_SPECLEN'],full_blocks))
    filename = os.path.join(file,extra)
    full_dynspec.readfilebinary(filename,0)

    if zoom:
        zoom_dynspec=bt.cutdynspec(full_dynspec,pulse=pulse)
        zoom_dynspec.par.station = beams['STATION_NAME'][0]
        incoherent_dynspec+=zoom_dynspec/6.
#        zoom_beam=bt.cutbeam(beam)
#        zoom_dynspec=bt.rawdyncalc(zoom_beam,axis_val=True)
        if verbose:
            super_plot(zoom_dynspec,pulse=pulse,zoom=zoom)
            pdb.set_trace()     # Used for now to pause and manually save the images.
    else:
        full_dynspec = full_dynspec.Transpose()
        dynspec=cr.hArray(float,[nblocks, 8193],full_dynspec[0:nblocks].vec())
        
        if tbin>1:
            dynspec = cr.hArray_toNumpy(dynspec)
            dynspec = dynspec.reshape((dynspec.shape[0]/tbin,tbin,dynspec.shape[1]))
            dynspec = np.sum(dynspec,axis=1)
            dynspec = cr.hArray(dynspec)

        dynspec = dynspec.Transpose()

        dynspec.par.station = beams['STATION_NAME'][0]
        dynspec.par.yvalues = beams['BEAM_FREQUENCIES']/10**6
        dynspec.par.xvalues = cr.hArray(float,nblocks,list(np.arange(0,nblocks*8192*2*5e-9,8192*2*5e-9)))

        if verbose:
            super_plot(dynspec,pulse=pulse,zoom=zoom)
            pdb.set_trace()     # Used for now to pause and manually save the images.

        incoherent_dynspec+=dynspec

if zoom:
    incoherent_dynspec.par = zoom_dynspec.par
    incoherent_dynspec.par.station = 'Superterp Stations (Incoherent)'
    super_plot(incoherent_dynspec,pulse=pulse,zoom=zoom)
    pdb.set_trace()
else:
    incoherent_dynspec.par.yvalues = beams['BEAM_FREQUENCIES']/10**6
    incoherent_dynspec.par.xvalues = cr.hArray(float,nblocks,list(np.arange(0,nblocks*8192*2*5e-9,8192*2*5e-9)))
    bt.dynplot(incoherent_dynspec)    
    cr.plt.title("Dynamic Spectrum")
#    cr.plt.imshow(cr.hArray_toNumpy(incoherent_dynspec),aspect='auto',cmap=cr.plt.cm.hot,origin='lower',vmin=1e-6,vmax=0.0008) # vmin=1e-5,vmax=0.003)
    cr.plt.ylabel("Frequency [MHz]")
    cr.plt.xlabel("Time [s]")
    pdb.set_trace()