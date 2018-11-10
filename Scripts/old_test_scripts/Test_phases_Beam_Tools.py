'''Script to find station delays by adding phased Beams.
'''

interactive =0

if not interactive:
    import matplotlib
    matplotlib.use('Agg')

import pycrtools as cr
import numpy as np
import sys; import os
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

#------------------

if interactive:
    cr.plt.ion()
cr.plt.rcParams['font.size']=10
verbose = 1

#filenames = ['L43784_D20120125T211154.887Z_CS002_R000_tbb.pol0.multi_station.beam','L43784_D20120125T211154.867Z_CS003_R000_tbb.pol0.multi_station.beam','L43784_D20120125T211154.887Z_CS004_R000_tbb.pol0.multi_station.beam','L43784_D20120125T211154.866Z_CS005_R000_tbb.pol0.multi_station.beam','L43784_D20120125T211154.871Z_CS006_R000_tbb.pol0.multi_station.beam','L43784_D20120125T211154.887Z_CS007_R000_tbb.pol0.multi_station.beam']
filenames = ['L43784_D20120125T211154.887Z_CS002_R000_tbb.pol1.multi_station.beam','L43784_D20120125T211154.867Z_CS003_R000_tbb.pol1.multi_station.beam','L43784_D20120125T211154.887Z_CS004_R000_tbb.pol1.multi_station.beam','L43784_D20120125T211154.866Z_CS005_R000_tbb.pol1.multi_station.beam','L43784_D20120125T211154.871Z_CS006_R000_tbb.pol1.multi_station.beam','L43784_D20120125T211154.887Z_CS007_R000_tbb.pol1.multi_station.beam']
filenames = filenames[:2] 

beams = cr.open(filenames)
beams['NCHUNKS'] = 191

factor = .1
delay_step = beams['TBB_SAMPLE_INTERVAL'][0]*factor
delay_range = np.arange(delay_step*110,delay_step*131,delay_step)
for i,dt in enumerate(delay_range):
    if abs(dt) < 1e-15:   #Removing unfeasible time resolution in delays
        delay_range[i] = 0            
SNR = cr.hArray(float,len(delay_range),fill=0)

for i,dt in enumerate(delay_range):
    beams['CAL_DELAY'] = [0.0,dt]
    TAB = bt.addBeams(beams,dm=26.776,verbose=interactive)
    dynspec,cleandynspec = bt.rawdyncalc(TAB,clean=True,verbose=interactive)

    zoom_dynspec = bt.cutdynspec(cleandynspec,pulse=0)
    flat_cleandyn = bt.flatdyn(zoom_dynspec)

    stats = dostats(flat_cleandyn)
    SNR[i] = stats[3]

    if verbose:
        print '----------------------------------'
        print 'Original Mean ', stats[0]
        print 'Min', flat_cleandyn.vec().min()
        print 'Max', stats[1]
        print 'Stddev ', stats[2]
        print '....................'
        print 'For '+str(len(filenames))+' stations, and delay ='+str(dt)
        print 'SNR', stats[3]
        print '----------------------------------'

#        if len(delay_range)!=1:
#            cr.plt.subplot(4,np.ceil(len(delay_range)/4.),i)
        cr.plt.plot(flat_cleandyn.par.xvalues,flat_cleandyn,label='SNR = '+str(SNR[i])+', for Delay = '+str(dt))
        cr.plt.xlabel("+/- Time [s]")

cr.plt.legend(loc=1)
#        pdb.set_trace()

if not interactive:
    cr.plt.savefig('Test.phases.vs.SNR.['+str(delay_range[0])+','+str(delay_range[-1])+'].'+beams['STATION_NAME'][0]+'_'+beams['STATION_NAME'][-1]+'peaks.ps',orientation='landscape',papertype='a4')
print 'SNRs : ', SNR
print 'Delay range: ', delay_range    


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
