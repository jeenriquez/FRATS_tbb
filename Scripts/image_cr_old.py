#! /usr/bin/env python

import pycrtools as cr
import os
import pdb;# pdb.set_trace()
from pycrtools.tasks import imager

filenames = ['L28348_D20110612T231913.199Z_CS002_R000_tbb.h5']

f = cr.open('~/RESEARCH/VHECR/Data/'+filenames[0])

f['SELECTED_DIPOLES']='odd'
f['BLOCKSIZE']=2**10

# Find RFI and bad antennas

rfi = 0#cr.trun("FindRFI", f = f, plot_prefix ='/Users/emilio/RESEARCH/VHECR/Results/RFI',nofblocks=1000,verbose=False,startblock=5500)#,save_plots=True)

#print rfi.dirty_channels
#print rfi.good_antennas

# Flag dirty channels (from RFI excission)
##dirty_channels = list(set(station.polarization['0']["dirty_channels"] + station.polarization['1']["dirty_channels"]))
##fft_data[..., dirty_channels] = 0

if rfi:
    f['SELECTED_DIPOLES'] = rfi.good_antennas

#Imaging
output = 'out_cr.x1.no-rfi.fft.fits'

cr.trun("Imager", data = f, intgrfreq = False, nblocks=1,ntimesteps=1,startblock=5998,NAXIS1=155,NAXIS2=155,output=output,CDELT1=-0.5,CDELT2=0.5,CTYPE1 ='ALON_SIN',CTYPE2 ='ALAT_SIN') #,rfi_remove=rfi.dirty_channels,inversefft=True,CRVAL1=231.9-360,CRVAL2=63.1,CTYPE1='ALON-AZP',CTYPE2='ALAT-AZP',FREQMIN=0.3e8,FREQMAX=0.8e8,LONPOLE=180,

