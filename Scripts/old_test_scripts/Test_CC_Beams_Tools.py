'''Script to test the cross correlation of pulses.
'''

import pycrtools as cr
import numpy as np
import sys; import os
import pdb;# pdb.set_trace()
import Beam_Tools as bt

filenames = ['L43784_D20120125T211154.887Z_CS002_R000_tbb.pol1.multi_station.beam','L43784_D20120125T211154.867Z_CS003_R000_tbb.pol1.multi_station.beam','L43784_D20120125T211154.887Z_CS004_R000_tbb.pol1.multi_station.beam','L43784_D20120125T211154.866Z_CS005_R000_tbb.pol1.multi_station.beam','L43784_D20120125T211154.871Z_CS006_R000_tbb.pol1.multi_station.beam','L43784_D20120125T211154.887Z_CS007_R000_tbb.pol1.multi_station.beam']

filenames = [filenames[0],filenames[0]]
#filenames = filenames[0:7:5]
beams = cr.open(filenames)

beams['DM'] = 26.76

block_range=np.arange(240,320)
nblocks = len(block_range)

blocklen = 16384
speclen = 8193

# Do mean of complex or timeseries
mean_fft =1


#pdb.set_trace()
if mean_fft:
    cross_correlation = cr.hArray(float,[1,16384]) #beams.empty('TIMESERIES_DATA')[0]
    fft_matrix = cr.hArray(complex,[nblocks,speclen])   # ,len(filenames)
else:
    cross_correlation = cr.hArray(float,[nblocks,blocklen])

#Frequency masking
Freq_Range=[149,169]
frequencies = cr.hArray(float,speclen, list(np.arange(1e+08,2e+08+12207.03125,12207.03125)))
frequency_slice = cr.hArray(int, 2)
frequencies/=10**6
cr.hFindSequenceBetweenOrEqual(frequency_slice,frequencies,Freq_Range[0], Freq_Range[1], 0, 0)

for i,block in enumerate(block_range):
    print 'Calculation at {0:.2%}  \r'.format(float(block-block_range[0])/nblocks),
    sys.stdout.flush()
    fft = beams.empty('FFT_DATA')
    beams.getFFTData(fft,int(block))
    # to correlate only the frequencies where the pulse is, and not all the strong RFI.
    fft[...,:frequency_slice[0]] = 0+0j
    fft[...,frequency_slice[1]:] = 0+0j    

    #Cross correlation: vec1*conj(vec2)
    fft[1:,...].crosscorrelatecomplex(fft[0], True)
#    fft[0:1,...].crosscorrelatecomplex(fft[0], True)
    
    if mean_fft:
        fft_matrix[i] = fft[1]    
    else:
        cc = beams.empty('TIMESERIES_DATA')
        cc[1].invfftw(fft[1])

        cross_correlation[i] = cc[1]
    
cross_correlation = cross_correlation.Transpose()
final_cc = cross_correlation[...].sum()

print 'Need to check that the code is working, now with either doing the station mean on the complex or on the timeseries.'
stop

#------------------------
final_cc/=blocklen*nblocks
#final_cc.abs()

final_cc=cr.hArray(float,len(final_cc),fill=final_cc)

cr.plt.ion()
time = np.arange(-1*blocklen*5e-9/2,blocklen*5e-9/2,5e-9)
cr.plt.plot(time,cr.hArray_toNumpy(final_cc))
cr.plt.xlabel('Time [s] ')

final_cc=cr.hArray(float,[1,len(final_cc)],fill=final_cc)




#------------------------
#Fitting  cc.
range_cc = 20
cut_cc = cr.hArray(float,[1,len(final_cc[0,speclen-range_cc:speclen+range_cc].vec())],final_cc[0,speclen-range_cc:speclen+range_cc].vec())

mx=cr.trun("FitMaxima",cut_cc,doplot=True,peak_width=40)#,splineorder=2)

print '--------------------'
print 'CC of ', beams['STATION_NAME']
print 'Calibration Delay = '+str(mx.lags[0]*beams['SAMPLE_INTERVAL'][0] - blocklen*beams['SAMPLE_INTERVAL'][0]/2 )+' [sec]'


'''#fft_matrix=fft_matrix.Transpose()
#mean_fft = cr.hArray(complex,[1,8193], fft_matrix[...].sum())
#cross_correlation.invfftw(mean_fft)

#cross_correlation/=nblocks
#cross_correlation.abs()'''

