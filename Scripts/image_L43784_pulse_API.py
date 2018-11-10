#! /usr/bin/env python
'''Image beams stored in
'''

import pycrtools as cr
from pytmf import *
import glob
import os
import pdb;# pdb.set_trace()
from pycrtools.tasks import imager
import Beam_Tools as bt

Obs_ID ='L43784'
Obs_ID2 ='L43784_D20120125T2111'
polarization=1

FRATS_ANALYSIS=os.environ["FRATS_ANALYSIS"].rstrip('/')+'/'

beam_suffix='*CS00[2-7]*pol%i?.new.multi_station.beam'%(polarization)
filenames=glob.glob(FRATS_ANALYSIS+'../testing/PSR_B0329.54/beam.results/beams.sub_station/'+beam_suffix)
beam_suffix2='*CS00[2-7]*pol%i?.new.beam'%(polarization)
filenames2=glob.glob(FRATS_ANALYSIS+'../testing/PSR_B0329.54/beam.results/beams.sub_station/'+beam_suffix2)

print 'Using files for calibrating: ',filenames
beams=cr.open(sorted(filenames))

print 'Using files for imaging: ',filenames2
beams2=cr.open(sorted(filenames2))

#----------------------
#Calibration

#PSR B0329+54 -- Pulsar    0.71452 sec
alpha = hms2deg(03,32,59.37); # Right assention
delta = dms2deg(54,34,44.9); # Declination
ctype1 ='RA---SIN'
ctype2 ='DEC--SIN'

beams['DM'] = 26.76
#ST_DELAYS = bt.ccBeams(beams,freq_range=[151,168],time_range=[.021,.0255])
ST_DELAYS = bt.ccBeams(beams,freq_range=[149,169],time_range=[.021299,.024494])
print 'ST_DELAYS', ST_DELAYS

#----------------------

beams2['DM'] = 26.76
beams2['CAL_DELAY']=ST_DELAYS
output = FRATS_ANALYSIS+Obs_ID2+'/products/'+'L43784.new.RA.n.DEC.fits'

#----------------------
#Imaging
#cr.trun("Imager", data = beams, intgrfreq = True, nblocks=25,ntimesteps=8,startblock=200,NAXIS1=65,NAXIS2=65,CDELT1=-.05,CDELT2=.05,FREQMIN=1.51e8,FREQMAX=1.68e8,CTYPE1=ctype1,CTYPE2=ctype2,output=output,CRVAL3=1.51e8,CRVAL1=alpha,CRVAL2=delta)
cr.trun("Imager", data = beams2, intgrfreq = True, nblocks=25,ntimesteps=1,startblock=275,NAXIS1=65,NAXIS2=65,CDELT1=-.05,CDELT2=.05,FREQMIN=1.5e8,FREQMAX=1.69e8,CTYPE1=ctype1,CTYPE2=ctype2,output=output,CRVAL3=1.5e8,CRVAL1=alpha,CRVAL2=delta)
