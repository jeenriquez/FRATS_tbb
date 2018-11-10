'''Script that adds Beams with given delays and checks for SNR.
'''

import pycrtools as cr
import numpy as np
import sys; import os; import glob
import pdb;# pdb.set_trace()
import Beam_Tools as bt


def dostats(data):
    stats = cr.hArray(float,4,fill=0)
    data = cr.hArray(float,len(data),data)

    lenght = len(data)
    mask1=lenght/5
    mask2=lenght/20

    data.sort()

    stats[0] = data[:-1*mask1].vec().mean()     # Mean without the peak.
    data-=stats[0]
    stats[1] = data[-1*mask2:].vec().mean()     # Max value
    stats[2] = data[:-1*mask1].vec().stddev()   # Stddev
    stats[3] = stats[1]/stats[2]                # SNR

    return stats

#------------------------------------------------------------------------

#Find files.
BEAMDATA=os.environ["BEAMDATA"].rstrip('/')+'/'
beams_dir = 'beams.sub_station/'

#------------------
#Options
verbose = 1
set = 2    # Set of files, cal delays and so on.
#------------------

#filenames = ['L43784_D20120125T211154.887Z_CS002_R000_tbb.pol0.multi_station.beam','L43784_D20120125T211154.867Z_CS003_R000_tbb.pol0.multi_station.beam','L43784_D20120125T211154.887Z_CS004_R000_tbb.pol0.multi_station.beam','L43784_D20120125T211154.866Z_CS005_R000_tbb.pol0.multi_station.beam','L43784_D20120125T211154.871Z_CS006_R000_tbb.pol0.multi_station.beam','L43784_D20120125T211154.887Z_CS007_R000_tbb.pol0.multi_station.beam']
#filenames = ['L43784_D20120125T211154.887Z_CS002_R000_tbb.pol1.multi_station.beam','L43784_D20120125T211154.867Z_CS003_R000_tbb.pol1.multi_station.beam','L43784_D20120125T211154.887Z_CS004_R000_tbb.pol1.multi_station.beam','L43784_D20120125T211154.866Z_CS005_R000_tbb.pol1.multi_station.beam','L43784_D20120125T211154.871Z_CS006_R000_tbb.pol1.multi_station.beam','L43784_D20120125T211154.887Z_CS007_R000_tbb.pol1.multi_station.beam']
#filenames = ['L43784_D20120125T211154.887Z_CS002_R000_tbb.pol1.multi_station.beam','L43784_D20120125T211154.867Z_CS003_R000_tbb.pol1.multi_station.beam','L43784_D20120125T211154.887Z_CS004_R000_tbb.pol1.multi_station.beam','L43784_D20120125T211154.871Z_CS006_R000_tbb.pol1.multi_station.beam','L43784_D20120125T211154.887Z_CS007_R000_tbb.pol1.multi_station.beam']
#filenames = ['L43784_D20120125T211154.887Z_CS002_R000_tbb.pol1.multi_station.beam','L43784_D20120125T211154.887Z_CS004_R000_tbb.pol1.multi_station.beam','L43784_D20120125T211154.887Z_CS007_R000_tbb.pol1.multi_station.beam']
#filenames = ['L43784_D20120125T211154.887Z_CS002_R000_tbb.pol1.multi_station.beam']#,'L43784_D20120125T211154.887Z_CS002_R000_tbb.pol1.multi_station.beam']
#filenames = ['L43784_D20120125T211154.887Z_CS002_R000_tbb.pol1a.multi_station.beam','L43784_D20120125T211154.887Z_CS002_R000_tbb.pol1b.multi_station.beam']
#filenames = filenames[0]

extra = {
        0: '.new.multi_station.',
        1: '.multi_station.',
        2: '.new.multi_station.',
        }

filenames = {
        0:['L43784_D20120125T211154.887Z_CS002_R000_tbb.pol1a'+extra[set]+'beam','L43784_D20120125T211154.867Z_CS003_R000_tbb.pol1a'+extra[set]+'beam','L43784_D20120125T211154.887Z_CS004_R000_tbb.pol1a'+extra[set]+'beam','L43784_D20120125T211154.866Z_CS005_R000_tbb.pol1a'+extra[set]+'beam','L43784_D20120125T211154.871Z_CS006_R000_tbb.pol1a'+extra[set]+'beam','L43784_D20120125T211154.887Z_CS007_R000_tbb.pol1a'+extra[set]+'beam','L43784_D20120125T211154.887Z_CS002_R000_tbb.pol1b'+extra[set]+'beam','L43784_D20120125T211154.867Z_CS003_R000_tbb.pol1b'+extra[set]+'beam','L43784_D20120125T211154.887Z_CS004_R000_tbb.pol1b'+extra[set]+'beam','L43784_D20120125T211154.866Z_CS005_R000_tbb.pol1b'+extra[set]+'beam','L43784_D20120125T211154.871Z_CS006_R000_tbb.pol1b'+extra[set]+'beam','L43784_D20120125T211154.887Z_CS007_R000_tbb.pol1b'+extra[set]+'beam'],
        1:['L43784_D20120125T211154.887Z_CS002_R000_tbb.pol1a'+extra[set]+'beam','L43784_D20120125T211154.867Z_CS003_R000_tbb.pol1a'+extra[set]+'beam','L43784_D20120125T211154.887Z_CS004_R000_tbb.pol1a'+extra[set]+'beam','L43784_D20120125T211154.866Z_CS005_R000_tbb.pol1a'+extra[set]+'beam','L43784_D20120125T211154.871Z_CS006_R000_tbb.pol1a'+extra[set]+'beam','L43784_D20120125T211154.887Z_CS007_R000_tbb.pol1a'+extra[set]+'beam','L43784_D20120125T211154.887Z_CS002_R000_tbb.pol1b'+extra[set]+'beam','L43784_D20120125T211154.867Z_CS003_R000_tbb.pol1b'+extra[set]+'beam','L43784_D20120125T211154.887Z_CS004_R000_tbb.pol1b'+extra[set]+'beam','L43784_D20120125T211154.866Z_CS005_R000_tbb.pol1b'+extra[set]+'beam','L43784_D20120125T211154.871Z_CS006_R000_tbb.pol1b'+extra[set]+'beam','L43784_D20120125T211154.887Z_CS007_R000_tbb.pol1b'+extra[set]+'beam'],
        2:['L43784_D20120125T211154.887Z_CS002_R000_tbb.pol0a'+extra[set]+'beam','L43784_D20120125T211154.867Z_CS003_R000_tbb.pol0a'+extra[set]+'beam','L43784_D20120125T211154.887Z_CS004_R000_tbb.pol0a'+extra[set]+'beam','L43784_D20120125T211154.866Z_CS005_R000_tbb.pol0a'+extra[set]+'beam','L43784_D20120125T211154.871Z_CS006_R000_tbb.pol0a'+extra[set]+'beam','L43784_D20120125T211154.887Z_CS007_R000_tbb.pol0a'+extra[set]+'beam','L43784_D20120125T211154.887Z_CS002_R000_tbb.pol0b'+extra[set]+'beam','L43784_D20120125T211154.867Z_CS003_R000_tbb.pol0b'+extra[set]+'beam','L43784_D20120125T211154.887Z_CS004_R000_tbb.pol0b'+extra[set]+'beam','L43784_D20120125T211154.866Z_CS005_R000_tbb.pol0b'+extra[set]+'beam','L43784_D20120125T211154.871Z_CS006_R000_tbb.pol0b'+extra[set]+'beam','L43784_D20120125T211154.887Z_CS007_R000_tbb.pol0b'+extra[set]+'beam'],
        }

#beams_suffix = '*pol0?'+extra[set]+'beam'
#filenames = glob.glob(BEAMDATA+beams_dir+beams_suffix)

beams = cr.open([BEAMDATA+beams_dir+file for file in filenames[set]])
beams['NCHUNKS'] = 191

#beams['CAL_DELAY']=[-5.4375e-09,-5.96875e-09,0.0,-2.9375e-09,9.375e-11,2.8125e-10]   # wrt CS007
#beams['CAL_DELAY']=[5.4375e-09,5.96875e-09,0.0,2.9375e-09,-9.375e-11,-2.8125e-10]   # wrt CS007 , negative
#beams['CAL_DELAY']=[0.0,6.875e-10,-5.4375e-09,-8.3125e-09,1.0625e-09,-4.9375e-09]   # wrt CS002 ,
#beams['CAL_DELAY']=[0.0,-6.875e-10,5.4375e-09,8.3125e-09,-1.0625e-09,4.9375e-09]   # wrt CS002 , negative
#beams['CAL_DELAY']=[0.0,-3.4375e-09,6.25000000003e-10,9.37500000001e-10,9.37500000001e-10,-1.84375e-08]   # wrt CS002 , CC values
#beams['CAL_DELAY']=[0.0,-2.5e-09,-3.34375e-09,-3.5625e-09,-3.6875e-09,-1.721875e-08]   # wrt CS002 , abs(CC) hFFTWResample
#beams['CAL_DELAY']=[0.0,2.5e-09,3.34375e-09,3.5625e-09,3.6875e-09,1.721875e-08]   # wrt CS002 , abs(CC) hFFTWResample , negative
#beams['CAL_DELAY']=[0.0,1.6875e-09,2.375e-09,-5.84375e-09,2.59375e-09,-1.3375e-08]   # wrt CS002 , CC hFFTWResample
#beams['CAL_DELAY']=[0.0,1.6859375e-09,2.3625e-09,-5.8421875e-09,2.6046875e-09,-1.3371875e-08]   # wrt CS002 , CC hFFTWResample, more acurate
#beams['CAL_DELAY']=[0.0,-1.6875e-09,-2.375e-09,5.84375e-09,-2.59375e-09,1.3375e-08]   # wrt CS002 , CC hFFTWResample , negative
#beams['CAL_DELAY']/=2.68   # Magic value :)
#beams['CAL_DELAY']=[0.0,6.4453125e-10,-5.37109375e-09,-8.4375e-09,-5.2734375e-09,-4.98046875e-09 ]   #wrt CS002 , CC hFFTWResample , upsampling
#beams['CAL_DELAY']*=-1
#beams['CAL_DELAY']=[0.0,3.12499999998e-10]  #wrt CS002a, , CC hFFTWResample , ears.
#beams['CAL_DELAY']=[0.0,7.61718749999e-10,-4.98046875e-09,-8.0078125e-09,-5.0390625e-09,-4.27734375e-09,2.92968750002e-10,7.51953124997e-10,-5.419921875e-09,-8.408203125e-09,-5.205078125e-09,-5.400390625e-09] #wrt CS002a, , CC hFFTWResample , ears, new beams with modified rel.ant.pos.


cal_delay = {
        0:[0.0,7.61718749999e-10,-4.98046875e-09,-8.0078125e-09,-5.0390625e-09,-4.27734375e-09,2.92968750002e-10,7.51953124997e-10,-5.419921875e-09,-8.408203125e-09,-5.205078125e-09,-5.400390625e-09], #wrt CS002a, , CC hFFTWResample , ears, new beams with modified rel.ant.pos.
        1:[0.0,8.00781249999e-10,-4.98046875e-09,-8.046875e-09,-5.05859375e-09,-4.21875e-09,3.12499999998e-10,7.81250000002e-10,-5.4296875e-09,-8.41796875e-09 ,-5.1953125e-09 ,-5.3515625e-09],  #wrt CS002a, , CC hFFTWResample , ears.
        2:[0.0,8.30078124997e-10,-4.326171875e-09,-7.890625e-09,-5.625e-09,-4.619140625e-09,-3.22265625e-10,3.02734375003e-10,-4.951171875e-09,-7.734375e-09,-5.107421875e-09,-5.087890625e-09] #wrt CS002a, , CC hFFTWResample , ears, new beams with modified rel.ant.pos.
        }

beams['CAL_DELAY']=cal_delay[set]

#------------------------------------------------------------------------

extra1 = 'wrt CS002a'

print 'CAL_DELAY', beams['CAL_DELAY']

TAB = bt.addBeams(beams,dm=26.76)
dynspec,cleandynspec = bt.rawdyncalc(TAB,clean=True)

freq_range=[149,169]
time_range = [0.0,0.0501]

zoom_dynspec = bt.cutdynspec(cleandynspec,freq_range=freq_range,time_range=time_range)

flat_cleandyn=bt.flatdyn(zoom_dynspec,axis='x')
stats = dostats(flat_cleandyn)
if verbose:
    print '----------------------------------'
    print 'Calculations for ' + extra1
    print '----------------------------------'
    print 'Original Mean ', stats[0]
    print 'Min', flat_cleandyn.vec().min()
    print 'Max', stats[1]
    print 'Stddev ', stats[2]
    print '....................'
    print 'For '+str(len(filenames[set]))+' stations'
    print 'SNR', stats[3]
    print '----------------------------------'

cr.plt.ion()
cr.plt.plot(flat_cleandyn)

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
