'''Script to test the cross correlation of pulses.
This script will be added as a function in Beam_Tools since it seems to work well.
'''

import pycrtools as cr
import numpy as np
import sys; import os; import glob
import pdb;# pdb.set_trace()
import Beam_Tools as bt
from pycrtools.tasks import pulsecal
import scipy.signal


#-----------------
print '--------------------'
#beams_dir = '/vol/astro1/lofar/frats/tbb/testing/PSR_B0329.54/beam.results/beams.sub_station/'
#beams_suffix = '*pol0?.new.multi_station.beam'
beams_dir = '/vol/astro1/lofar/frats/tbb/analysis/L43784_D20120125T2111/beam.results/'
beams_suffix = '*pol0.*LOFAR_centered.beam'
filenames = glob.glob(beams_dir+beams_suffix)

#filenames = [filenames[6],filenames[7]]
filenames = [filenames[1],filenames[0]]

beams = cr.open(filenames)
print 'Files used : ', filenames

beams['DM'] = 26.76
#beams['CAL_DELAY'] = [0,1.38e-9]
#print 'cal-delay', beams['CAL_DELAY']

#block_range=np.arange(3662,3906) # Where pulse 2 is located.
#block_range=np.arange(7250,7280) # Where pulse 1 is located.
block_range=np.arange(260,300) # Where pulse 3 is located.
print 'Block range : ' , block_range
nblocks = len(block_range)

blocklen = 16384
speclen = 8193

mean_at_fft = 0
envelope = 0


cross_correlation = cr.hArray(float,[1,blocklen]) #beams.empty('TIMESERIES_DATA')[0]
fft_matrix = cr.hArray(complex,[nblocks,speclen])   # ,len(filenames)
ts_matrix = cr.hArray(complex,[nblocks,blocklen*2])

invfftplan = cr.FFTWPlanManyDftC2r(blocklen, 1, 1, 1, 1, 1, cr.fftw_flags.ESTIMATE)

#-----------------
#Frequency masking
#Freq_Range=[140,155]   #Pulse2
Freq_Range=[149,169]  #Pulse3
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
    cc_one = cr.hArray(float,blocklen*2,fill=0)
    beams.getFFTData(fft,int(block))
    # to correlate only the frequencies where the pulse is, and not all the strong RFI.
    fft[...,:frequency_slice[0]] = 0+0j
    fft[...,frequency_slice[1]:] = 0+0j
    #Cross correlation: vec1*conj(vec2)
    fft4cc=cr.hArray(complex,[1,8193],fft[1])
    fft[0].crosscorrelatecomplex(fft4cc, True)

#    fft[0:1,...].crosscorrelatecomplex(fft[1], True)
    #------------------------
    #Magic
    fft*=-1
    fft4cc*=-1

#    crosscorr=trerun('CrossCorrelateAntennas',"crosscorr",timeseries_data_cut_to_pulse,pardict=par,oversamplefactor=10)
    if mean_at_fft:
#        fft_matrix[i] = fft4cc
        fft_matrix[i] = fft[0]
    else:
    #------------------------
    #Upsampling
        crosscorr = cr.hArray(complex,blocklen+1,fill=0.0)
#        crosscorr[speclen:] =  fft4cc
        crosscorr[speclen:] =  fft[0]
        cc_one.invfftw(crosscorr)
        ts_matrix[i] = cc_one

#------------------------
#Mean cross correlation of blocks within the pulse time-width
if mean_at_fft:
    fft_matrix=fft_matrix.Transpose()
    mean_fft = cr.hArray(complex,[1,8193], fft_matrix[...].sum())
    mean_fft /= fft_matrix.shape()[0]

#------------------------
#Back to timeseries
    #mean_fft[...].nyquistswap(2)
    #cr.hFFTWExecutePlan(cross_correlation[...],mean_fft[...],invfftplan)

    cross_correlation.invfftw(mean_fft[0])
    #cross_correlation.abs()
else:
    ts_matrix=ts_matrix.Transpose()
    cross_correlation = cr.hArray(float,[1,blocklen*2], ts_matrix[...].sum())

cross_correlation/=blocklen

#cross_correlation*=-1

cross_correlation.reshape([cross_correlation.shape()[-1]])

print 'Initial cross corellation finished.'

#------------------------
# Resample
resample = 32
delay_step = beams['TBB_SAMPLE_INTERVAL'][0] / resample

cc_smooth = cr.hArray(float,[1,blocklen*resample])
cr.hFFTWResample(cc_smooth,cross_correlation)

# Ploting resampled
cr.plt.ion()
cc_smooth=cr.hArray(float,len(cc_smooth),fill=cc_smooth)

time = np.arange(-1*blocklen*5e-9/2,blocklen*5e-9/2,5e-9/resample)
cr.plt.plot(time,cr.hArray_toNumpy(cc_smooth))
cr.plt.xlabel('Time [s] ')

#cc_smooth.abs()
cr.plt.plot(time,cr.hArray_toNumpy(cc_smooth))

#------------------------
#Using a Hilbert envelope

if envelope:
    cr.plt.ion()
    cr.plt.figure()
    cr.plt.plot(cc_smooth)

    cc_nparray = cc_smooth.toNumpy()
    envelope_cc = abs(scipy.signal.hilbert(cc_nparray))
    cc_smooth = cr.hArray(envelope_cc)

    cr.plt.plot(cc_smooth,'r')

#------------------------
#Fitting  cc.
#range_init = speclen*len(delay_range)
#range_cc = int(range_init/50)
#cut_cc = cr.hArray(float,[1,len(final_cc[0,range_init-range_cc:range_init+range_cc].vec())],final_cc[0,range_init-range_cc:range_init+range_cc].vec())

#cut_cc.runningaverage(25,cr.hWEIGHTS.GAUSSIAN)
cc_smooth=cr.hArray(float,[1,len(cc_smooth)],fill=cc_smooth)
cut_cc = cc_smooth

mx=cr.trun("FitMaxima",cut_cc,doplot=True,peak_width=2000,splineorder=2)

print '--------------------'
print 'CC of ', beams['STATION_NAME']
print 'Calibration Delay = ' + str(mx.lags[0]*delay_step - delay_step*blocklen/2*resample)+' [sec]'


#------------------------#------------------------#------------------------#------------------------#------------------------

#260*(2**14*5e-9) = 0.0212992
#299*(2**14*5e-9) = 0.02449408

ST_DELAYS = bt.ccBeams(beams,freq_range=[149,169],time_range=[.021299,.024494],verbose=1)

print '--------------------'
print 'Now, doing the same calculation with ccBeams in Beam_Tools.py:'
print 'Calibration Delay = ' + str(ST_DELAYS[1])+' [sec]'

