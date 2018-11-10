'''Script to test delay addition, and  cross correlation in synthetic data (time series with a single pulse).
'''
import numpy as np
from pycrtools import *
import pdb;# pdb.set_trace()

casa=0
nyquist=1 # 0 if ny1, 1 if ny2.

beams=open('L43784_D20120125T211154.887Z_CS002_R000_tbb.pol1.multi_station.beam')
blocklen=2**14#2**10
speclen = blocklen/2+1
ts_test = hArray(float,blocklen)#,[0,1])
a= range(5)*3
ts_test[10:15]=a
#ts_test[40:45]=a
a.reverse()
ts_test[17:22]=a
#ts_test[47:52]=a
ts_test[15]=90
ts_test[16]=559
#ts_test[46]=90
#ts_test[45]=559


fft_test = hArray(complex,blocklen/2+1)
if casa:
    hFFTCasa(fft_test,ts_test,2)
else:
    fftplan = FFTWPlanManyDftR2c(blocklen, 1, 1, 1, 1, 1, fftw_flags.ESTIMATE)
    hFFTWExecutePlan(fft_test,ts_test, fftplan)
    if nyquist:
        fft_test.nyquistswap(2)
#        hReverse(fft_test)
auto_fft = hArray(complex,[3,blocklen/2+1])
auto_fft[0] = fft_test
auto_fft[2] = fft_test

#------------------------
freqs=hArray(float,beams['BEAM_FREQUENCIES'],beams['BEAM_FREQUENCIES']-1e8+nyquist*1e8)
#hReverse(freqs)

#------------------------
t_step = 3.21e-9
#delay = np.arange(0,10e-9+t_step,t_step)
delay = t_step
delay = hArray(delay)

for i,dt in enumerate(list(delay)):
    print 'delay', dt
    if not dt:
        continue
    #------------------------
    weights = beams.empty('FFT_DATA')
    phases=hArray(float,weights,fill=0)
    phases.delaytophase(freqs,hArray(dt*(-1.)))
    weights.phasetocomplex(phases)
    auto_fft[1] = fft_test
    auto_fft[1].mul(weights)
    #------------------------
#    plt.subplot(int(np.ceil(len(delay)/4.)),4,i)
    
    ts_plot=0
    if ts_plot:
        factor=2
#        fft_single = hArray(complex,auto_fft[1].shape()[-1]*factor,auto_fft[1])
        fft_single = hArray(complex,blocklen+1,fill=0.0)
        
        if factor==2:
#            fft_single[0:speclen-1] = 0.0
            fft_single[speclen:] =  auto_fft[1]
#            fft_single = fft_single[1:]
#        fft_single.nyquistswap(2)
        ts2=hArray(float,blocklen*factor)
        ts2.invfftw(fft_single)
        ts2/=blocklen*(-1)
        resample=1
        plt.plot(ts_test,label='original')
#        ts2.abs()
        time = np.arange(0,blocklen,1./factor)
        plt.plot(time,ts2,label=str(dt)+'[s] shifted')
        plt.xlim(10,25)
        resample = 16
        cc_smooth = hArray(float,[blocklen*resample])
        hFFTWResample(cc_smooth,ts2)
        time2 = np.arange(0,blocklen,1./resample)
        plt.plot(time2,cc_smooth,label='resampled')
        expected_delay = 16+(dt)/5e-9
        plt.plot([expected_delay,expected_delay],[min(ts_test),max(ts_test)],label='expected_delay')
        plt.xlabel('Sample #')
        plt.legend(loc=3)

if ts_plot:
    assert False

#------------------------
#Cross-correlation
if casa:
    #auto_fft[1] = auto_fft[0]  * complex_conjugate(auto_fft[1])
    raise NotImplementedError("Need to do this by hand.")
else:
    auto_fft[1:,...].crosscorrelatecomplex(auto_fft[0], True)

#-------
#Magic...

auto_fft*=-1

#-------
# Auto-correlation
#crosscorr = hArray(complex,auto_fft[2].shape()[-1],auto_fft[2])
crosscorr = hArray(complex,blocklen+1,fill=0.0)
crosscorr[speclen:] =  auto_fft[2]

what_now=hArray(float,blocklen*2)
if casa:
    what_now.invfftcasa(crosscorr,2)
else:
    what_now.invfftw(crosscorr)
resample=2
time = np.arange(-1*blocklen*5e-9/2,blocklen*5e-9/2,5e-9/resample)
plt.ion()
plt.plot(time,what_now)

#-------
#Delay correlation
#crosscorr = hArray(complex,auto_fft[1].shape()[-1],auto_fft[1])
crosscorr = hArray(complex,blocklen+1,fill=0.0)
crosscorr[speclen:] =  auto_fft[1]

if casa:
    what_now.invfftcasa(crosscorr,2)
else:
    what_now.invfftw(crosscorr)
plt.plot(time,what_now)

#------------------------
# Resample
resample = 16
delay_step = beams['TBB_SAMPLE_INTERVAL'][0] / resample

cc_smooth = hArray(float,[1,blocklen*resample])
hFFTWResample(cc_smooth,what_now)

# Ploting resampled
cc_smooth=hArray(float,len(cc_smooth),fill=cc_smooth)

time = np.arange(-1*blocklen*5e-9/2,blocklen*5e-9/2,5e-9/resample)
plt.plot(time,hArray_toNumpy(cc_smooth))
plt.xlabel('Time [s] ')

cc_smooth=hArray(float,[1,len(cc_smooth)],fill=cc_smooth)
#cc_smooth.abs()

#------------------------
#Fitting  cc.
mx=trun("FitMaxima",cc_smooth,doplot=True,peak_width=200,splineorder=2)

print '--------------------'
print 'CC of ', beams['STATION_NAME']
print 'Calibration Delay = ' + str(mx.lags[0]*delay_step - delay_step*blocklen/2*resample)+' [sec]'
