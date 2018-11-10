#! /usr/bin/env python

import pycrtools as cr
import os
import pdb;# pdb.set_trace()

LOFARDATA=os.environ["LOFARDATA"].rstrip('/')+'/'

filenames = ['L43784_D20120125T211154.887Z_CS002_R000_tbb.pol1.multi_station.beam','L43784_D20120125T211154.887Z_CS002_R000_tbb.pol1.multi_station.beam']

f = cr.open([LOFARDATA+fname for fname in filenames])

f['CAL_DELAY']=[0.0,1.00e-12]
f['DM'] = 26.76
output = 'out_pul.sub_st.high_res.fits'

pdb.set_trace()
cr.trun("Imager", data = f, intgrfreq = False, nblocks=30,ntimesteps=10,startblock=100,NAXIS1=33,NAXIS2=33,FREQMIN=1.5e8,FREQMAX=1.69e8,CDELT1=-.125,CDELT2=.125,CTYPE1 ='ALON_SIN',CTYPE2 ='ALAT_SIN',output=output)
