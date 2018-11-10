#! /usr/bin/env python
'''Images from tbbs.
'''

import pycrtools as cr
from pytmf import *
import os
import pdb;# pdb.set_trace()

LOFARDATA=os.environ["LOFARDATA"].rstrip('/')+'/'
BEAMOUT=os.environ["BEAMOUT"].rstrip('/')+'/'
event_dir ='L43784_ev1/'

filename = 'L43784_D20120125T211154.887Z_CS002_R000_tbb.h5'

f = cr.open(LOFARDATA+event_dir+filename)

RA = 2
#----------------------------------------------------
#RFI excission.
# Find RFI and bad antennas

rfi = cr.trun("FindRFI", f=f,nofblocks=100,verbose=False,startblock=0,save_plots=False)

#print rfi.dirty_channels
#print rfi.good_antennas

if rfi:
    f['SELECTED_DIPOLES'] = rfi.good_antennas

#PSR B0329+54 -- Pulsar    0.71452 sec
#----------------------------------------------------
nblocks=25
ntimesteps=8
startblock=200
cdelt = -0.05
naxis = 65
#,CRVAL3=1.5e8

f['BLOCKSIZE'] = 2**16

if RA==1:
    alpha = hms2deg(03,32,59.37); # Right assention
    delta = dms2deg(54,34,44.9); # Declination
    ctype1 ='RA---SIN'
    ctype2 ='DEC--SIN'
elif RA==2:
    alpha = 0.
    delta = 90.
    ctype1 ='ALON_AIT'
    ctype2 ='ALAT_AIT'
    cdelt = -1
    naxis = 65
    nblocks=20
    ntimesteps=5
    startblock=200
    dt = f['BLOCKSIZE']*f['SAMPLE_INTERVAL'][0]
    mask = rfi.dirty_channels
else:
    alpha = 289.06051721551938
    delta = 68.773413970138805
    ctype1 ='ALON_SIN'
    ctype2 ='ALAT_SIN'


output = BEAMOUT+'out_pul.single_st.DMd.RA.n.DEC.fits'

#cr.trun("Imager", data=f, intgrfreq=True, nblocks=nblocks, ntimesteps=ntimesteps, startblock=startblock, NAXIS1=naxis, NAXIS2=naxis, CDELT1=cdelt, CDELT2=cdelt, FREQMIN=1.5e8, FREQMAX=1.69e8, CTYPE1=ctype1, CTYPE2=ctype2, output=output, CRVAL1=alpha, CRVAL2=delta)
cr.trun("ImagerDM", data=f, DM=26.76, intgrfreq=True, nblocks=nblocks, ntimesteps=ntimesteps, startblock=startblock, NAXIS1=naxis, NAXIS2=naxis, CDELT1=cdelt, CDELT2=cdelt, FREQMIN=1.5e8, FREQMAX=1.69e8, CTYPE1=ctype1, CTYPE2=ctype2, output=output, CRVAL1=alpha, CRVAL2=delta, dt=dt, mask=mask)
