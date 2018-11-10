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
old_beams = 1

FRATS_ANALYSIS=os.environ["FRATS_ANALYSIS"].rstrip('/')+'/'

if old_beams:
    beam_suffix='*CS00[2-7]*pol%i?.new.multi_station.beam'%(polarization)
    filenames=glob.glob(FRATS_ANALYSIS+'../testing/PSR_B0329.54/beam.results/beams.sub_station/'+beam_suffix)
    beam_suffix2='*CS00[2-7]*pol%i?.new.beam'%(polarization)
    filenames2=glob.glob(FRATS_ANALYSIS+'../testing/PSR_B0329.54/beam.results/beams.sub_station/'+beam_suffix2)
else:
    beam_suffix='*pol%i.*HBA?.LOFAR_centered.beam'%(polarization)
    filenames=glob.glob(FRATS_ANALYSIS+Obs_ID2+'/beam.results/'+beam_suffix)
    beam_suffix2='*pol%i.*HBA?.beam'%(polarization)
    filenames2=glob.glob(FRATS_ANALYSIS+Obs_ID2+'/beam.results/'+beam_suffix2)

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

#Calculated long ago with old beams, but these values made a good image.
#ST_DELAYS = cr.hArray(float, [12], fill=[-8.0078125e-09,-8.408203125e-09,7.61718749999e-10,7.51953124997e-10,-5.0390625e-09,-5.205078125e-09,0.0,2.92968750002e-10,-4.98046875e-09,-5.419921875e-09,-4.27734375e-09,-5.400390625e-09])

#Calculated with the very old beams.
#ST_DELAYS = cr.hArray(float, [12], fill=[-7.92969e-09,-7.74414e-09,8.20313e-10,2.92969e-10,-5.66406e-09,-5.14648e-09,0,-3.32031e-10,-4.33594e-09,-4.9707e-09,-4.63867e-09,-5.10742e-09])

#Calculated from the newest beams from (wrong)
#ST_DELAYS = cr.hArray(float, [12], fill=[3.25391e-08,3.2334e-08,-4.88281e-08,-4.96387e-08,-1.55303e-07,-1.54365e-07,0,3.71094e-10,-9.41113e-08,-9.47266e-08,-1.54023e-07,-1.55303e-07])
#Calculated from the newest beams from
#hArray(float, [12], fill=[-4.67773e-09,-4.72656e-09,-7.13867e-09,-3.21289e-09,6.64063e-10,-7.05078e-09,0,6.83594e-10,5.66406e-10,1.46484e-10,-7.42188e-09,-4.08203e-09])

#Calculated from the newest beams, and after fixing the subsample residual (from clock offset calibration).
#ST_DELAYS = cr.hArray(float, [12], fill=[-2.53906e-10,-6.05469e-10,-5.56641e-10,-1.05469e-09,-1.2793e-09,-1.70898e-09,0,4.88281e-11,1.15234e-09,3.80859e-10,-6.44531e-10,-1.5332e-09]) # len=12 slice=[0:12])

#New beams, with no residuals... (so should use in this way..).
#ST_DELAYS = hArray(float, [12], fill=[-8.69141e-09,-2.68555e-09,-3.83789e-09,-4.62891e-09,-2.83203e-10,-4.30664e-09,0,6.60156e-09,8.69141e-10,2.63672e-10,-4.04297e-09,-5.30273e-09])

#----------------------

beams2['DM'] = 26.76
beams2['CAL_DELAY']=ST_DELAYS
output = FRATS_ANALYSIS+Obs_ID2+'/products/'+'L43784.new.RA.n.DEC.fits'


#----------------------
#Imaging
#cr.trun("Imager", data = beams, intgrfreq = True, nblocks=25,ntimesteps=8,startblock=200,NAXIS1=65,NAXIS2=65,CDELT1=-.05,CDELT2=.05,FREQMIN=1.51e8,FREQMAX=1.68e8,CTYPE1=ctype1,CTYPE2=ctype2,output=output,CRVAL3=1.51e8,CRVAL1=alpha,CRVAL2=delta)
cr.trun("Imager", data = beams2, intgrfreq = True, nblocks=25,ntimesteps=1,startblock=275,NAXIS1=65,NAXIS2=65,CDELT1=-.05,CDELT2=.05,FREQMIN=1.5e8,FREQMAX=1.69e8,CTYPE1=ctype1,CTYPE2=ctype2,output=output,CRVAL3=1.5e8,CRVAL1=alpha,CRVAL2=delta)
