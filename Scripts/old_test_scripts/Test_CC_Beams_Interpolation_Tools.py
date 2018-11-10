'''Script to test the cross correlation of pulses.
'''

import pycrtools as cr
import numpy as np
import sys; import os
import pdb;# pdb.set_trace()
import Beam_Tools as bt
#-----------------
filenames = ['L43784_D20120125T211154.887Z_CS002_R000_tbb.pol1.multi_station.beam','L43784_D20120125T211154.867Z_CS003_R000_tbb.pol1.multi_station.beam','L43784_D20120125T211154.887Z_CS004_R000_tbb.pol1.multi_station.beam','L43784_D20120125T211154.866Z_CS005_R000_tbb.pol1.multi_station.beam','L43784_D20120125T211154.871Z_CS006_R000_tbb.pol1.multi_station.beam','L43784_D20120125T211154.887Z_CS007_R000_tbb.pol1.multi_station.beam']

filenames = [filenames[0],filenames[0]]
#filenames = filenames[0:7:5]
beams = cr.open(filenames)

beams['DM'] = 26.76

#block_range=np.arange(240,320) # Where pulse 3 is located.
block_range=np.arange(260,300) # Where pulse 3 is located.
nblocks = len(block_range)

blocklen = 16384
speclen = 8193

cross_correlation = cr.hArray(float,[1,16384]) #beams.empty('TIMESERIES_DATA')[0]
fft_matrix = cr.hArray(complex,[nblocks,speclen])   # ,len(filenames)

#-----------------
#Frequency masking
Freq_Range=[149,169]
frequencies = cr.hArray(float,speclen, list(np.arange(1e+08,2e+08+12207.03125,12207.03125)))
frequency_slice = cr.hArray(int, 2)
frequencies/=10**6
cr.hFindSequenceBetweenOrEqual(frequency_slice,frequencies,Freq_Range[0], Freq_Range[1], 0, 0)

#pdb.set_trace()
#-----------------
for i,block in enumerate(block_range):
    print 'CC calculation at {0:.2%}  \r'.format(float(block-block_range[0])/nblocks),
    sys.stdout.flush()
    fft = beams.empty('FFT_DATA')
    beams.getFFTData(fft,int(block))
    # to correlate only the frequencies where the pulse is, and not all the strong RFI.
    fft[...,:frequency_slice[0]] = 0+0j
    fft[...,frequency_slice[1]:] = 0+0j    

    #Cross correlation: vec1*conj(vec2)
    fft[1:,...].crosscorrelatecomplex(fft[0], True)
#    fft[0:1,...].crosscorrelatecomplex(fft[0], True)
    
    fft_matrix[i] = fft[1]    
#------------------------
#Mean cross correlation of blocks within the pulse time-width
fft_matrix=fft_matrix.Transpose()
mean_fft = cr.hArray(complex,[1,8193], fft_matrix[...].sum())
mean_fft /= fft_matrix.shape()[0]
print 'Initial cross corellation finished.'    
#------------------------
# Interpolating cc values.
factor = 0.1
delay_step = beams['TBB_SAMPLE_INTERVAL'][0]*factor
delay_range = np.arange(0,delay_step/factor,delay_step)
for i,dt in enumerate(delay_range):
    if abs(dt) < 1e-15:   #Removing unfeasible time resolution in delays
        delay_range[i] = 0            
    
#pdb.set_trace()
big_ts = cr.hArray(float,[len(delay_range),blocklen])

for i,dt in enumerate(delay_range):
    print 'Interpolation (phasing) calculation at {0:.2%}  \r'.format(float(i)/len(delay_range)),
    sys.stdout.flush()
    shift_fft_matrix = cr.hArray(complex,mean_fft.shape(),mean_fft)
    ts = cr.hArray(float,blocklen)

    delays = cr.hArray(float,2,fill=[dt])
    
    weights = cr.hArray(complex,speclen)
    phases=cr.hArray(float,weights,fill=0)
    phases.delaytophase(beams['BEAM_FREQUENCIES'],delays)
    weights.phasetocomplex(phases)
    shift_fft_matrix[...].mul(weights[...])
    ts.invfftw(shift_fft_matrix)

    big_ts[i] = ts

print 'Interpolation (phasing) is done.'

final_cc = cr.hArray(float,[len(delay_range)*blocklen])
interpolated_dt_range = range(len(delay_range)*mean_fft.shape()[0])

#Arranging interpolated timeseries data.
for f in range(blocklen):
    print 'Interpolation (arranging) calculation at {0:.2%}  \r'.format(float(f)/float(blocklen)),
    sys.stdout.flush()
    for dt in interpolated_dt_range:    
        final_cc[f*len(interpolated_dt_range)+dt] = big_ts[dt,f]

print 'Interpolation (arranging) is done.'
print 'Iterpolation finished.'    
#pdb.set_trace()
#------------------------
final_cc/=blocklen#*nblocks
final_cc.abs()

final_cc=cr.hArray(float,len(final_cc),fill=final_cc)

cr.plt.ion()
time = np.arange(-1*blocklen*5e-9/2,blocklen*5e-9/2,5e-9*factor)
cr.plt.plot(time,cr.hArray_toNumpy(final_cc))
cr.plt.xlabel('Time [s] ')

final_cc=cr.hArray(float,[1,len(final_cc)],fill=final_cc)
#pdb.set_trace()
#------------------------
#Fitting  cc.
range_init = speclen*len(delay_range)
range_cc = int(range_init/50)
cut_cc = cr.hArray(float,[1,len(final_cc[0,range_init-range_cc:range_init+range_cc].vec())],final_cc[0,range_init-range_cc:range_init+range_cc].vec())

cut_cc.runningaverage(25,cr.hWEIGHTS.GAUSSIAN)

mx=cr.trun("FitMaxima",cut_cc,doplot=True,peak_width=50,splineorder=2)

print '--------------------'
print 'CC of ', beams['STATION_NAME']
print 'Calibration Delay = '+str((range_init-range_cc+mx.lags[0])*beams['SAMPLE_INTERVAL'][0]*factor - blocklen*len(delay_range)*beams['SAMPLE_INTERVAL'][0]*factor/2 )+' [sec]'


'''#fft_matrix=fft_matrix.Transpose()
#mean_fft = cr.hArray(complex,[1,8193], fft_matrix[...].sum())
#cross_correlation.invfftw(mean_fft)

#cross_correlation/=nblocks
#cross_correlation.abs()'''

