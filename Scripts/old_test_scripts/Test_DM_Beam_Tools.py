'''Script to test DM vs SNR.
'''

interactive =1

if not interactive:
    import matplotlib
    matplotlib.use('Agg')

import pycrtools as cr
import numpy as np
import sys; import os
import pdb;# pdb.set_trace()
import Beam_Tools as bt
#------------------

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

cr.plt.rcParams['font.size']=20

filenames = ['L43784_D20120125T211154.887Z_CS002_R000_tbb.pol1.multi_station.beam']

beams = cr.open(filenames)
beams['NCHUNKS'] = 191

verbose = 1

DM_range = np.arange(26.81,26.9,.005)
#DM_range = [26.776]
SNR = cr.hArray(float,len(DM_range),fill=0)
max = cr.hArray(float,len(DM_range),fill=0)
stddev = cr.hArray(float,len(DM_range),fill=0)

for i,dm in enumerate(DM_range):
    dynspec,cleandynspec = bt.dyncalc(beams=beams,clean=True,dm=float(dm),verbose=interactive)
    zoom_dynspec = bt.cutdynspec(cleandynspec,pulse=0)
    flat_cleandyn = bt.flatdyn(zoom_dynspec)
    flat_cleandyn -= flat_cleandyn[0:200].vec().mean()
    
    stats = dostats(flat_cleandyn)
    SNR[i] = stats[3]
    max[i] = stats[1]
    stddev[i] = stats[2]

    if verbose:
        print '----------------------------------'
        print 'Mean ', stats[0]
        print 'Min', flat_cleandyn.vec().min()
        print 'Max', stats[1]
        print 'Stddev ', stats[2]
        print '....................'
        print 'For '+str(len(filenames))+' stations, and DM ='+str(dm)
        print 'SNR', stats[3]
        print '----------------------------------'

#        if len(DM_range)!=1:
#            cr.plt.subplot(4,np.ceil(len(DM_range)/4.),i)
        cr.plt.plot(flat_cleandyn.par.xvalues,flat_cleandyn,label='SNR = '+str(SNR[i])+', for DM = '+str(dm))
cr.plt.autoscale(axis='x',tight=True)
        cr.plt.xlabel("+/- Time [s]")

cr.plt.legend()
#        pdb.set_trace()

if not interactive:
    cr.plt.savefig('Test.DM.vs.SNR.['+str(DM_range[0])+','+str(DM_range[-1])+'].ps',orientation='landscape',papertype='a4')
print 'SNRs : ', SNR
print 'DM ragne: ', DM_range    


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
