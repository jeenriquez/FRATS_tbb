from pycrtools import *
beam_suffix='*pol0*CasA*'
import Beam_Tools as bt
filenames=glob.glob(beam_suffix)
beams=open(sorted(filenames))

freq_range=[41,80]
time_range=[0.6,1.1]
verbose=1
ref_station='CS002'
inv_base_fit=1


if 'LBA' in beams['ANTENNA_SET'][0]:
#        raise NotImplementedError('This code is optimized for HBA (Nyquist 2) FFTs. Since invfftw is bugged for Nyquist Zone 2. Also, reference antena is for HBA0.')
    antenna_set = beams['ANTENNA_SET'][0]
    factor = [0,1]
else:
    antenna_set = 'HBA0'
    factor = [1,2]

#Find the reference station index
try:
    i = np.arange(beams['NOF_BEAM_DATASETS'])[(np.array(beams['STATION_NAME'])==ref_station) & (np.array(beams['ANTENNA_SET'])==antenna_set)]
    ref_station_num = int(i[0])
except:
    print 'Refence station not found. Using station: ' + beams['STATION_NAME'][0]
    ref_station_num = 0

#Frequency selection indexing.
if type(freq_range) != type([]) or len(freq_range)!=2 or freq_range[0] > freq_range[1] or freq_range[0]<0 or freq_range[1]<0:
    raise ValueError('Need a list of lenght 2, with first element "<" or "=" second element. Both elements positive.')
frequencies = beams.getFrequencies()
frequency_slice = cr.hArray(int, 2)
frequencies.setUnit('M','Hz')
cr.hFindSequenceBetweenOrEqual(frequency_slice,frequencies,freq_range[0], freq_range[1], 0, 0)

#Time selection indexing.
if type(time_range) != type([]) or len(time_range)!=2 or time_range[0] > time_range[1] or time_range[0]<0 or time_range[1]<0:
    raise ValueError('Need a list of lenght 2, with first element "<" or "=" second element. Both elements positive.')
times = beams.getTimes()
times.setUnit('','s')
block_range = cr.hArray(int, 2)
cr.hFindSequenceBetweenOrEqual(block_range,times,time_range[0], time_range[1], 0, 0)
nblocks = block_range[1] - block_range[0]

#Large time series matrix
ts_matrix = cr.hArray(float,(nblocks,beams['NOF_BEAM_DATASETS'],beams['BEAM_BLOCKLEN']*factor[1]))

t0=time.clock()
#-----------------
#CC loop
for i in range(nblocks):
    if verbose:
        print 'CC calculation at {0:.2%}  \r'.format(float(i)/nblocks),
        sys.stdout.flush()
    #Correlating frequecies only where pulse is present.

    #Gaussian smothing of the frequency band.
    gaussian_weights = cr.hArray(cr.hGaussianWeights(int(beams["BLOCKSIZE"]/100), 8.0))
    bandpass_filter = cr.hArray(float,len(frequencies), fill=1)
    bandpass_filter[:frequency_slice[0]] = 0
    bandpass_filter[frequency_slice[1]:] = 0
    cr.hRunningAverage(bandpass_filter, gaussian_weights)

#        fft_matrix[...,:frequency_slice[0]] = 0+0j
#        fft_matrix[...,frequency_slice[1]:] = 0+0j

    #Gettind the data.
    fft_matrix = beams.empty('FFT_DATA')
    beams.getFFTData(fft_matrix,int(block_range[0]+i))

    #For frequency gain calibration. (Temporarily here.)
#        fft_matrix *= inv_base_fit

    # Apply bandpass
    fft_matrix[...].mul(bandpass_filter)

    #Cross correlation: vec1*conj(vec2)
    fft4cc=cr.hArray(complex,[1,beams['BEAM_SPECLEN']],fft_matrix[ref_station_num])
    for beam in np.arange(beams['NOF_BEAM_DATASETS'])[(np.array(range(beams['NOF_BEAM_DATASETS']))!=ref_station_num)]:
        fft_matrix[int(beam)].crosscorrelatecomplex(fft4cc, True)
    fft_matrix[ref_station_num].crosscorrelatecomplex(fft4cc, True)   #Autocorrelation.

    #------------------------
    #There is a difference in defintion between ffts (here fftw and fftcasa).
    #If the fftw definition is used, then one of the FFT data vectors still needs to be multiplied by -1,1,-1,...
    #to get the peak of the crosscorrelation for lag = 0 in the middle of floatvec.
    #This makes sense since the lag can be positive or negative.
    fft_matrix*=-1

    #------------------------
    #Upsampling
    ts_cc = cr.hArray(float,[beams['NOF_BEAM_DATASETS'],beams['BEAM_BLOCKLEN']*factor[1]],)
    long_fft = cr.hArray(complex,[beams['NOF_BEAM_DATASETS'],beams['BEAM_SPECLEN']*factor[1]-factor[0]])
    for beam in range(beams['NOF_BEAM_DATASETS']):
        long_fft[beam,beams['BEAM_SPECLEN']:] = fft_matrix[beam]
        ts_cc[beam].invfftw(long_fft[beam])
    ts_matrix[i] = ts_cc

if verbose:
    print "Finished - ccLoop in %f sec"%(time.clock()-t0)

CAL_DELAY = cr.hArray(float,beams['NOF_BEAM_DATASETS'])

#-----------------
#Finding time delay - loop.
for beam in range(beams['NOF_BEAM_DATASETS']):
    #------------------------
    #Single block CC averaging.
    cross_correlation = cr.hArray(float,[beams['NOF_BEAM_DATASETS'],beams['BEAM_BLOCKLEN']*factor[1]])
    cc = cr.hArray(float,[1,beams['BEAM_BLOCKLEN']*factor[1]])

    for block in range(nblocks):
        cc += ts_matrix[block,beam]
    cross_correlation[beam] = cc / nblocks

    #------------------------
    # Resampling
    resample = 32
    delay_step = beams['TBB_SAMPLE_INTERVAL'][0] / resample

    cc_smooth = cr.hArray(float,[1,beams['BEAM_BLOCKLEN']*resample])
    cr.hFFTWResample(cc_smooth,cross_correlation[beam])

    #------------------------
    #Fitting  cc.
    mx=cr.trun("FitMaxima",cc_smooth,peak_width=200,splineorder=2,doplot=verbose,newfigure=verbose)
    delay = mx.lags[0]*delay_step - delay_step*beams['BEAM_BLOCKLEN']/2*resample
    CAL_DELAY[beam] = delay
    if verbose:
        print '--------------------'
        print 'CC of %s with %s'%(beams['FILENAMES'][beam],beams['FILENAMES'][ref_station_num])
        print 'Calibration Delay = ' + str(mx.lags[0]*delay_step - delay_step*beams['BEAM_BLOCKLEN']/2*resample)+' [sec]'

