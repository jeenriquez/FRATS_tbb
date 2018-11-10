#! /usr/bin/env python
'''Image beams stored in scratch/lofar/frats/PSR_B0329.54/beam.results/beams.*/
'''

import pycrtools as cr
from pycrtools.tasks import imager
from pytmf import *
import os
import pdb;# pdb.set_trace()

BEAMDATA=os.environ["BEAMDATA"].rstrip('/')+'/'
BEAMOUT=os.environ["BEAMOUT"].rstrip('/')+'/'

beams_dir = 'beams.sub_station/'

extra = '.new.'

#filenames = ['L43784_D20120125T211154.866Z_CS005_R000_tbb.pol0.beam', 'L43784_D20120125T211154.871Z_CS006_R000_tbb.pol0.beam', 'L43784_D20120125T211154.887Z_CS004_R000_tbb.pol0.beam', 'L43784_D20120125T211154.867Z_CS003_R000_tbb.pol0.beam', 'L43784_D20120125T211154.887Z_CS002_R000_tbb.pol0.beam', 'L43784_D20120125T211154.887Z_CS007_R000_tbb.pol0.beam']
#filenames = ['L43784_D20120125T211154.866Z_CS005_R000_tbb.pol1.beam', 'L43784_D20120125T211154.871Z_CS006_R000_tbb.pol1.beam', 'L43784_D20120125T211154.887Z_CS004_R000_tbb.pol1.beam', 'L43784_D20120125T211154.867Z_CS003_R000_tbb.pol1.beam', 'L43784_D20120125T211154.887Z_CS002_R000_tbb.pol1.beam', 'L43784_D20120125T211154.887Z_CS007_R000_tbb.pol1.beam']
#filenames = ['L43784_D20120125T211154.866Z_CS005_R000_tbb.pol1.multi_station.beam', 'L43784_D20120125T211154.871Z_CS006_R000_tbb.pol1.multi_station.beam', 'L43784_D20120125T211154.887Z_CS004_R000_tbb.pol1.multi_station.beam', 'L43784_D20120125T211154.867Z_CS003_R000_tbb.pol1.multi_station.beam', 'L43784_D20120125T211154.887Z_CS002_R000_tbb.pol1.multi_station.beam', 'L43784_D20120125T211154.887Z_CS007_R000_tbb.pol1.multi_station.beam']
filenames = ['L43784_D20120125T211154.887Z_CS002_R000_tbb.pol1a'+extra+'beam','L43784_D20120125T211154.867Z_CS003_R000_tbb.pol1a'+extra+'beam','L43784_D20120125T211154.887Z_CS004_R000_tbb.pol1a'+extra+'beam','L43784_D20120125T211154.866Z_CS005_R000_tbb.pol1a'+extra+'beam','L43784_D20120125T211154.871Z_CS006_R000_tbb.pol1a'+extra+'beam','L43784_D20120125T211154.887Z_CS007_R000_tbb.pol1a'+extra+'beam','L43784_D20120125T211154.887Z_CS002_R000_tbb.pol1b'+extra+'beam','L43784_D20120125T211154.867Z_CS003_R000_tbb.pol1b'+extra+'beam','L43784_D20120125T211154.887Z_CS004_R000_tbb.pol1b'+extra+'beam','L43784_D20120125T211154.866Z_CS005_R000_tbb.pol1b'+extra+'beam','L43784_D20120125T211154.871Z_CS006_R000_tbb.pol1b'+extra+'beam','L43784_D20120125T211154.887Z_CS007_R000_tbb.pol1b'+extra+'beam']
#filenames = ['L43784_D20120125T211154.887Z_CS002_R000_tbb.pol1a'+extra+'beam','L43784_D20120125T211154.867Z_CS003_R000_tbb.pol1a'+extra+'beam','L43784_D20120125T211154.887Z_CS004_R000_tbb.pol1a'+extra+'beam','L43784_D20120125T211154.866Z_CS005_R000_tbb.pol1a'+extra+'beam','L43784_D20120125T211154.871Z_CS006_R000_tbb.pol1a'+extra+'beam','L43784_D20120125T211154.887Z_CS007_R000_tbb.pol1a'+extra+'beam']
#filenames = ['L43784_D20120125T211154.887Z_CS002_R000_tbb.pol1b'+extra+'beam','L43784_D20120125T211154.867Z_CS003_R000_tbb.pol1b'+extra+'beam','L43784_D20120125T211154.887Z_CS004_R000_tbb.pol1b'+extra+'beam','L43784_D20120125T211154.866Z_CS005_R000_tbb.pol1b'+extra+'beam','L43784_D20120125T211154.871Z_CS006_R000_tbb.pol1b'+extra+'beam','L43784_D20120125T211154.887Z_CS007_R000_tbb.pol1b'+extra+'beam']

f = cr.open([BEAMDATA+beams_dir+fname for fname in filenames])

f['CAL_DELAY']=[0.0,7.61718749999e-10,-4.98046875e-09,-8.0078125e-09,-5.0390625e-09,-4.27734375e-09,2.92968750002e-10,7.51953124997e-10,-5.419921875e-09,-8.408203125e-09,-5.205078125e-09,-5.400390625e-09] #wrt CS002a, , CC hFFTWResample , ears, new beams with modified rel.ant.pos.
#f['CAL_DELAY']=[0.0,8.00781249999e-10,-4.98046875e-09,-8.046875e-09,-5.05859375e-09,-4.21875e-09,3.12499999998e-10,7.81250000002e-10,-5.4296875e-09,-8.41796875e-09 ,-5.1953125e-09 ,-5.3515625e-09]
#f['CAL_DELAY']=[0.0,8.00781249999e-10,-4.98046875e-09,-8.046875e-09,-5.05859375e-09,-4.21875e-09]
#f['CAL_DELAY']=[3.12499999998e-10,7.81250000002e-10,-5.4296875e-09,-8.41796875e-09 ,-5.1953125e-09 ,-5.3515625e-09]
#f['CAL_DELAY']=[0.0,6.4453125e-10,-5.37109375e-09,-8.4375e-09,-5.2734375e-09,-4.98046875e-09]

RA = 1

#PSR B0329+54 -- Pulsar    0.71452 sec
if RA:
    alpha = hms2deg(03,32,59.37); # Right assention
    delta = dms2deg(54,34,44.9); # Declination
    ctype1 ='RA---SIN'
    ctype2 ='DEC--SIN'
else:
    alpha = 289.06051721551938
    delta = 68.773413970138805
    ctype1 ='ALON_SIN'
    ctype2 ='ALAT_SIN'

f['DM'] = 26.76

output = BEAMOUT+'out_pul.sub_st.cal.new.bugtest.RA.n.DEC.fits'

#cr.trun("Imager", data = f, intgrfreq = True, nblocks=25,ntimesteps=8,startblock=200,NAXIS1=65,NAXIS2=65,CDELT1=-.05,CDELT2=.05,FREQMIN=1.5e8,FREQMAX=1.69e8,CTYPE1=ctype1,CTYPE2=ctype2,output=output,CRVAL3=1.5e8,CRVAL1=alpha,CRVAL2=delta)
cr.trun("Imager", data = f, intgrfreq = True, nblocks=25,ntimesteps=1,startblock=275,NAXIS1=65,NAXIS2=65,CDELT1=-.05,CDELT2=.05,FREQMIN=1.5e8,FREQMAX=1.69e8,CTYPE1=ctype1,CTYPE2=ctype2,output=output,CRVAL3=1.5e8,CRVAL1=alpha,CRVAL2=delta)
