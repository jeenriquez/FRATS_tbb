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

Obs_ID ='L74100'
polarization=0

FRATS_ANALYSIS=os.environ["FRATS_ANALYSIS"].rstrip('/')+'/'

##beam_suffix='*pol%i_*HBA?.beam'%(polarization)
#beam_suffix='*pol0_sample*'
beam_suffix='*CS00[2-7]*pol0_sample*'
filenames=glob.glob(FRATS_ANALYSIS+Obs_ID+'_new'+'/beam.results/'+beam_suffix)
#filenames = ['L74100_D20121107T191144.924Z_CS005_R000_tbb.pol0_sample_calibrated_HBA0.LOFAR_centered.beam', 'L74100_D20121107T191144.993Z_CS030_R000_tbb.pol0_sample_calibrated_HBA1.LOFAR_centered.beam', 'L74100_D20121107T191144.917Z_CS007_R000_tbb.pol0_sample_calibrated_HBA1.LOFAR_centered.beam', 'L74100_D20121107T191144.993Z_CS030_R000_tbb.pol0_sample_calibrated_HBA0.LOFAR_centered.beam', 'L74100_D20121107T191144.929Z_CS011_R000_tbb.pol0_sample_calibrated_HBA0.LOFAR_centered.beam', 'L74100_D20121107T191144.991Z_CS021_R000_tbb.pol0_sample_calibrated_HBA1.LOFAR_centered.beam', 'L74100_D20121107T191144.916Z_CS004_R000_tbb.pol0_sample_calibrated_HBA0.LOFAR_centered.beam', 'L74100_D20121107T191144.941Z_CS017_R000_tbb.pol0_sample_calibrated_HBA0.LOFAR_centered.beam', 'L74100_D20121107T191144.929Z_CS011_R000_tbb.pol0_sample_calibrated_HBA1.LOFAR_centered.beam', 'L74100_D20121107T191144.964Z_CS002_R000_tbb.pol0_sample_calibrated_HBA1.LOFAR_centered.beam', 'L74100_D20121107T191144.924Z_CS005_R000_tbb.pol0_sample_calibrated_HBA1.LOFAR_centered.beam', 'L74100_D20121107T191144.941Z_CS017_R000_tbb.pol0_sample_calibrated_HBA1.LOFAR_centered.beam', 'L74100_D20121107T191144.964Z_CS002_R000_tbb.pol0_sample_calibrated_HBA0.LOFAR_centered.beam', 'L74100_D20121107T191144.917Z_CS007_R000_tbb.pol0_sample_calibrated_HBA0.LOFAR_centered.beam', 'L74100_D20121107T191144.991Z_CS021_R000_tbb.pol0_sample_calibrated_HBA0.LOFAR_centered.beam', 'L74100_D20121107T191144.916Z_CS004_R000_tbb.pol0_sample_calibrated_HBA1.LOFAR_centered.beam']
#filenames = ['L74100_D20121107T191144.917Z_CS007_R000_tbb.pol0_sample_calibrated_HBA1.LOFAR_centered.beam', 'L74100_D20121107T191144.993Z_CS030_R000_tbb.pol0_sample_calibrated_HBA0.LOFAR_centered.beam','L74100_D20121107T191144.991Z_CS021_R000_tbb.pol0_sample_calibrated_HBA1.LOFAR_centered.beam',   'L74100_D20121107T191144.941Z_CS017_R000_tbb.pol0_sample_calibrated_HBA0.LOFAR_centered.beam','L74100_D20121107T191144.941Z_CS017_R000_tbb.pol0_sample_calibrated_HBA1.LOFAR_centered.beam','L74100_D20121107T191144.964Z_CS002_R000_tbb.pol0_sample_calibrated_HBA0.LOFAR_centered.beam','L74100_D20121107T191144.917Z_CS007_R000_tbb.pol0_sample_calibrated_HBA0.LOFAR_centered.beam',          'L74100_D20121107T191144.991Z_CS021_R000_tbb.pol0_sample_calibrated_HBA0.LOFAR_centered.beam']

print 'Using files for imaging: ',filenames
#f = cr.open([FRATS_ANALYSIS+Obs_ID+'_new'+'/beam.results/'+filename for filename in filenames])
f = cr.open(filenames)

#PSR B0329+54 -- Pulsar    0.71452 sec
alpha = hms2deg(03,32,59.37); # Right assention
delta = dms2deg(54,34,44.9); # Declination
ctype1 ='RA---SIN'
ctype2 ='DEC--SIN'

f['DM'] = 26.76
#f['CAL_DELAY']=cr.hArray(float, [8], fill=[1.00586e-09,6.32813e-09,-7.97852e-09,9.86328e-09,5.83984e-09,0,-4.26758e-09,3.55469e-09]) # len=24 slice=[0:24])

ST_DELAYS = bt.ccBeams(f,freq_range=[175,195],time_range=[0.346,0.3495],verbose=1)
f['CAL_DELAY']=ST_DELAYS

output = FRATS_ANALYSIS+Obs_ID+'_new'+'/L74100.cal_on.fits'

cr.trun("Imager", data = f, intgrfreq = True, nblocks=50,ntimesteps=7,startblock=252*16,NAXIS1=65,NAXIS2=65,CDELT1=-.05,CDELT2=.05,FREQMIN=1.75e8,FREQMAX=1.95e8,CTYPE1=ctype1,CTYPE2=ctype2,output=output,CRVAL3=1.5e8,CRVAL1=alpha,CRVAL2=delta)
