"""
This module has a number of tools useful in managing TBB beamformed data (beams). These tools are used by the FRATS pipelines, but can also be used independently.
"""

import math; import time
import sys; import os
import pycrtools as cr
import numpy as np
import pdb;# pdb.set_trace()
import copy
from pycrtools.tasks import pulsecal


#NOTE:
#I'm organizing these functions in a better way.
#Options:
#   Make a taks (class) wich contains them.
#   Make a taks for each function
#   Make one module with multilpe tasks

def addBeams(beams, dyncalc=False,save_file=False,fraction=False,clean=False,tbin=1,dm=0,verbose=1,incoherent=False):
    '''Add complex beams toguether.
    If requested:
        - Calculates a dynamic spectrum.
        - Dedisperses the beam.
        - Bins the dynspec in integer time blocks (not for the TAB, which stays in complex).

    =============== ===== ===================================================================
    *beams*               Input beams, opened by the beam.py interphase.
    *dyncalc*       True  If True it calculates the dynamic spectrum
    *dm*            0     Dedispersion Measure
    *save_file*     False Saves a file in hArray format with additional information besides the array values.
    *fraction*      None  If not None, then a list of the form [x,y] such that extracting the fraction x/y of the data. with x>-1, and y>=x. This also makes tbin=1
    *tbin*          1     If >1 integrates over this number of blocks. Or time binning.
    *clean*         False If True it calculates the cleaned spectrum.
    *verbose*       1     If want to print status and performance.
    *incoherent*    False Adding beams coherently (default) or not.
    =============== ===== ===================================================================

    Example::

        import Beam_Tools as bt
        TAB,dynspec,cleandynspec = bt.addBeams(beams,dyncalc=True,tbin=16,dm=26.76,clean=True)

    '''
    if incoherent:
        print 'Incoherent addition of the following beams: '
    else:
        print 'Coherent addition of the following beams: '

    for file in beams['FILENAMES']:
        print file
    print '------_-------'

    if save_file:
        raise KeyError("Keyword is invalid for now: "+key)

    if fraction:
        if tbin > 1:
            print 'WARNING: tbin has been reset to 1, since using fraction.'
        tbin = 1

    t0=time.clock()

    speclen = beams['BLOCKSIZE']/2+1
    block_duration = beams['BLOCKSIZE']*beams['SAMPLE_INTERVAL'][0]
    nblocks = beams['NCHUNKS']*beams['BEAM_NBLOCKS']
    beam = beams.empty('FFT_DATA')
    total_beam = cr.hArray(complex,(nblocks,beam.shape()[1]))
    fft_matrix = cr.hArray(complex,beam.shape()[1])
    if not beams['DM'] and dm:
        beams['DM'] = dm

    if not fraction:
        fraction=[1,1]
    else:
        if type(fraction) != type([]) or len(fraction)!=2 or fraction[0] > fraction[1] or fraction[0]<0 or fraction[1]<0:
            raise ValueError('Need a list of lenght 2, with first element "<" or "=" second element. Both elements positive.')
        if fraction[0]==0: fraction[0]=1
        nblocks = int(nblocks/fraction[1])

    if dyncalc:
        dynspec = cr.hArray(float,(nblocks,beam.shape()[1]))

    if incoherent:
        for block in range(nblocks):
            bm = cr.hArray(complex,beam.shape()[1])
            if verbose:
                print 'Addition beams at {0:.2%}  \r'.format(float(block)/nblocks),
            sys.stdout.flush()
            beams.getFFTData(beam,block)
            dynspec[block].spectralpower2(beam[...])
        total_beam = dynspec
    else:
        for block in range(nblocks):
            bm = cr.hArray(complex,beam.shape()[1])
            if verbose:
                print 'Addition beams at {0:.2%}  \r'.format(float(block)/nblocks),
            sys.stdout.flush()
            beams.getFFTData(beam,block)
            for b in range(beam.shape()[0]):
                fft_matrix[:] = beam[b].vec()/float(b+1.)
                if b>0:
                    bm*=(b)/(b+1.)
                bm+=fft_matrix   #Around here I could possibly change the code to be saving each station at the time... just like beamformer does.. to help not eating too much memory.
            total_beam[block] = bm
            if dyncalc:
                dynspec[block].spectralpower2(bm)

    if verbose:
        print "Finished - addBeams in %f sec"%(time.clock()-t0)

    if dyncalc:
        start_time=(fraction[0]-1)*(block_duration*nblocks)
        end_time=fraction[0]*(block_duration*nblocks)

        #Create a frequency vector
        frequencies = beams['BEAM_FREQUENCIES']

        #Create a time vector
        times = cr.hArray(float,int(round((end_time-start_time)/(block_duration*tbin))),name="Time",units=("","s"))
        times.fillrange(start_time,block_duration*tbin)

        #Time integration.
        if tbin>1:
            dynspec = cr.hArray_toNumpy(dynspec)
            dynspec = dynspec.reshape((dynspec.shape[0]/tbin,tbin,dynspec.shape[1]))
            dynspec = np.sum(dynspec,axis=1)
            dynspec = cr.hArray(dynspec)

        #Cleaning dynamic spectrum.
        if clean:
            avspec=cr.hArray(float,speclen,fill=0.0)
            dynspec[...].addto(avspec)
            cleandynspec = dynspec/avspec

        #Transposing arrays.
        dynspec = dynspec.Transpose()
        if clean:
            cleandynspec = cleandynspec.Transpose()

        #Adding parameters
        dynspec.par.yvalues= frequencies
        dynspec.par.xvalues= times
        dynspec.par.tbin = tbin
        if clean:
            cleandynspec.par.yvalues= frequencies
            cleandynspec.par.xvalues= times
            cleandynspec.par.tbin = tbin

        if clean:
            return total_beam,dynspec, cleandynspec
        else:
            return total_beam,dynspec

    else:
        return total_beam

def cutDynspec(dynspec,freq_range=None,time_range=None,hdr=None):
    '''Cuts a dynamic spectrum around the dedispersed peak.

    ============== ===== ===================================================================
    *dynspec*            Dynamic spectrum produced either by addBeams or rawdyncalc.
    *freq_range*   None  Frequency range in MHz
    *time_range*   None  Time range in seconds
    ============== ===== ===================================================================

    Example::

        import Beam_Tools as bt
        zoom_dynspec = bt.cutDynspec(cleandynspec,freq_range=[175,195],time_range=[0.1,0.2])

    '''
    #-----------
    #Log of things to check:
    #1)tbin divisiion
    #-----------

    info = copy.deepcopy(dynspec.par)
    dynspec = dynspec.Transpose()
    nblocks = dynspec.shape()[0]

    if hdr:
        speclen = hdr['FFTSIZE']
        blocklen = hdr['BLOCKSIZE']*info.tbin
        dt_sample = cr.asval(hdr['SAMPLE_INTERVAL'])
        f0 = hdr['BEAM_FREQUENCIES'][0]
        f1 = hdr['BEAM_FREQUENCIES'][-1]
        df = cr.asval(hdr['FREQUENCY_INTERVAL'])
    else:
        speclen = 8193
        blocklen = (8192)*2
        dt_sample = 5e-9
        f0 = 1e+08; f1 = 2e+08; df = 12207.03125

    dt = dt_sample*blocklen*info.tbin

    frequencies = info.yvalues
    if not freq_range:
        frequency_slice  = [0,dynspec.shape()[1]]
        nchannels = dynspec.shape()[1]
    else:
        freq_range = cr.hArray(float,[2],fill=freq_range,units=('M','Hz')).setUnit("","Hz")
        frequency_slice = cr.hArray(int, 2)
        cr.hFindSequenceBetweenOrEqual(frequency_slice,frequencies,freq_range[0], freq_range[1], 0, 0)
        nchannels = frequency_slice[1] - frequency_slice[0]

    times = info.xvalues
    if not time_range:
        time_slice = [0,dynspec.shape()[0]]
        nblocks = dynspec.shape()[0]
    else:
        time_slice = cr.hArray(int, 2)
        cr.hFindSequenceBetweenOrEqual(time_slice,times,time_range[0], time_range[1], 0, 0)
        nblocks = time_slice[1] - time_slice[0]

    #First cutting the time axis.
    sliced_dynspec = cr.hArray(float,[nblocks,speclen],dynspec[time_slice[0]:time_slice[1]])

    sliced_dynspec = sliced_dynspec.Transpose()

    #Then cutting the frequency axis.
    zoom_dynspec = cr.hArray(float,[nchannels,nblocks],sliced_dynspec[frequency_slice[0]:frequency_slice[1]])

    times = cr.hArray(float,nblocks,list(np.arange(time_slice[0]*dt,time_slice[1]*dt,dt)))
    frequencies = cr.hArray(float,nchannels,frequencies[frequency_slice[0]:frequency_slice[1]].vec())

    zoom_dynspec.par.xvalues = times
    zoom_dynspec.par.yvalues = frequencies
    zoom_dynspec.par.tbin = info.tbin

    return zoom_dynspec

def ccBeams(beams,ref_station='CS002',freq_range=[130,160],time_range=[0,0.1],verbose=False,antenna_set='',inv_base_fit=1):
    ''' Cross correlate the dedispersed pulses.

    ============== ========= ===================================================================
    *ref_station*  CS002     Reference station to which all the other will be cross correlated.
    *freq_range*   [130,160] List with the begin and end of frequency selection. [f1,f2] in MHz
    *time_range*   [0,0.1]   List with the begin and end of time selection. [t1,t2] in s
    *verbose*      False     Prints and plots releavant information.
    *antenna_set*  None      Antenna set of the station to which all the other will be CCed, specially used for HBA1  substations.
    ============== ========= ===================================================================

    Example::
        import Beam_Tools as bt
        ST_DELAYS = bt.ccBeams(beams,freq_range=[175,195],time_range=[.345,.355],verbose=1)

    '''

    if not antenna_set:
        if 'LBA' in beams['ANTENNA_SET'][0]:
            raise NotImplementedError('This code is optimized for HBA (Nyquist 2) FFTs. Since invfftw is bugged for Nyquist Zone 2. Also, reference antena is for HBA0.')
            antenna_set = beams['ANTENNA_SET'][0]
            factor = [0,1,0]
        else:
            antenna_set = 'HBA0'
            factor = [1,2,beams['BEAM_SPECLEN']]

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
    freq_range = cr.hArray(float,[2],fill=freq_range,units=('M','Hz')).setUnit("","Hz")
    frequencies = beams.getFrequencies()
    frequency_slice = cr.hArray(int, 2)
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
            long_fft[beam,factor[2]:] = fft_matrix[beam]
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

    return CAL_DELAY

def doStats(data,verbose=True,nof_St=0):
    ''' Do simple SNR calculation.

    ============== ===== ===================================================================
    *data*               One dimentional vector with the frequency integrated spectrum.
    *verbose*      True  Print the SNR results.
    *nof_St*          0  Use str(beams['NOF_STATION_BEAM']), to show the number of stations.
    ============== ===== ===================================================================

    Example::

        import Beam_Tools as bt
        stats = bt.doStats(flat_cleandyn,verbose=True,nof_St=str(beams['NOF_BEAM_DATASETS']))

    '''
    stats = cr.hArray(float,6,fill=0)
    data = cr.hArray(float,len(data),data)

    lenght = len(data)
    mask1=int(lenght/5)
    mask2=[int(lenght/100) if lenght/100 >3 else 3 ][0]

    data.sort()

    stats[0] = data[:-1*mask1].vec().mean()     # Mean without the peak.
    data-=stats[0]
    stats[1] = data[-1*mask2:].vec().mean()     # Max value
    stats[2] = data[:-1*mask1].vec().stddev()   # Stddev
    stats[3] = stats[1]/stats[2]                # SNR
    stats[4] = data.vec().min()
    stats[5] = data.vec().max()

    if verbose:
        print '----------------------------------'
        print 'SNR calculations.'
        print '----------------------------------'
        print 'Original Mean ', stats[0]
        print 'Min - OM', stats[4]
        print 'Max - OM', stats[1]
        print 'Stddev ', stats[2]
        print '....................'
        print 'For %i stations'%(nof_St)
        print 'SNR (max_range/stddev)', stats[3]
        print 'SNR (max/stddev)', stats[5]/stats[2]
        print '----------------------------------'

    return stats

#=================BAR OF GREATNESS ========= anything below needs to be double checked / cleanned / optimized.

#cr.plt.rcParams['font.size']=20

def flatDynspec(dynspec,axis='x',verbose=False, down_sample=0):
    '''Integrates the dynspec over time or freq.

    ================== ===== ===================================================================
    *dynspec*                Dynamic spectrum produced either by addBeams, cutDynspec or rawdyncalc.
    *axis*             'x'   Axis to integrate over, 'x' for freq integration, 'y' for time integration.
    *verbose*          False To plot or print extra info.
    *down_sample*      0     Downsample factor, if 0, then not downsampling.
    ================== ===== ===================================================================

    Example::

        import Beam_Tools as bt
        flat_cleandyn=bt.flatDynspec(zoom_dynspec,axis='x', verbose=True, down_sample=20)

    '''
    #-----------
    #Things to check log:
    #1)
    #-----------

    ax = cr.hArray()

    if axis =='y':
        ax.par.xvalues =dynspec.par.yvalues
        x_label = 'Frequency [Hz]'
    elif axis =='x':
        ax.par.xvalues =dynspec.par.xvalues
        dynspec = dynspec.Transpose()
        x_label = 'Time [s]'

    ax.par.flattype = axis

    flat_dynspec = cr.hArray(dynspec[...].sum())
    flat_dynspec.par = ax.par

    if down_sample:
        #Downsampling
        flat_sampled = cr.hDownsample(flat_dynspec.vec(),down_sample)
        flat_sampled = cr.hArray(float,len(flat_dynspec.par.xvalues),flat_sampled)
        flat_sampled.par = ax.par
        flat_dynspec = flat_sampled

    if verbose:
        cr.plt.ion()
        cr.plt.clf()
        cr.plt.plot(flat_dynspec.par.xvalues,flat_dynspec.vec())
        cr.plt.xlabel(x_label)

    return flat_dynspec

def rawdyncalc(TAB,fraction=None,tbin=1,clean=False,axis_val=False,verbose=1):
    '''
    Calculates the dynamic spectrum of a raw beam (no TBB info).

    =============== ===== ===================================================================
    *TAB*                 Input complex array.
    *fraction*      None  If not None, then a list of the form [x,y] such that extracting the fraction x/y of the data. with x>-1, and y>=x.
    *tbin*          1     If >1 integrates over this number of blocks. Or time binning.
    *clean*         False If True it calculates the cleaned spectrum.
    *verbose*       1     If want to print status and performance.
    =============== ===== ===================================================================

    Example::

        import Beam_Tools as bt
        dynspec = bt.rawdyncalc(filename)

    or::

        import Beam_Tools as bt
        dynspec,cleandynspec = bt.rawdyncalc(TAB,tbin=16,clean=True)

    '''

    t0=time.clock()

    speclen = TAB.shape()[1]
    block_duration = (TAB.shape()[1]-1)*2*5e-09
    nblocks = TAB.shape()[0]

    #Create a frequency vector
    if axis_val:
        frequencies = TAB.par.yvalues
    else:
        frequencies = cr.hArray(float,speclen, list(np.arange(1e+08,2e+08+12207.03125,12207.03125)))

    if not fraction:
        fraction=[1,1]
    else:
        if type(fraction) != type([]) or len(fraction)!=2 or fraction[0] > fraction[1] or fraction[0]<0 or fraction[1]<0:
            raise ValueError('Need a list of lenght 2, with first element "<" or "=" second element. Both elements positive.')
        if fraction[0]==0: fraction[0]=1
        nblocks = int(nblocks/fraction[1])

    start_time=(fraction[0]-1)*(block_duration*nblocks)/fraction[1]
    end_time=fraction[0]*(block_duration*nblocks)/fraction[1]

    dynspec = cr.hArray(float,[nblocks,speclen])
    block_range = range((fraction[0]-1)*nblocks,(fraction[0])*nblocks)
    beam = cr.hArray(complex,speclen)

    #Reading TAB.
    for block in block_range:
        if verbose:
            print 'Dynspec calculation at {0:.2%}  \r'.format(float(block)/nblocks),
        sys.stdout.flush()
        beam = TAB[block]
        dynspec[block].spectralpower2(beam)

    #Time integration.
    if tbin>1:
        dynspec = cr.hArray_toNumpy(dynspec)
        dynspec = dynspec.reshape((dynspec.shape[0]/tbin,tbin,dynspec.shape[1]))
        dynspec = np.sum(dynspec,axis=1)
        dynspec = cr.hArray(dynspec)

    #Cleaning dynamic spectrum.
    if clean:
        avspec=cr.hArray(float,speclen,fill=0.0)
        dynspec[...].addto(avspec)
        cleandynspec = dynspec/avspec

    #Transposing arrays.
    dynspec = dynspec.Transpose()
    if clean:
        cleandynspec = cleandynspec.Transpose()

    #Create a time vector
    if axis_val:
        times = TAB.par.xvalues
    else:
        times = cr.hArray(float,int(round((end_time-start_time)/(block_duration*tbin))),name="Time",units=("","s"))
        times.fillrange(start_time,block_duration*tbin)

    #Adding parameters
    dynspec.par.yvalues= frequencies
    dynspec.par.xvalues= times
    dynspec.par.tbin = tbin
    if clean:
        cleandynspec.par.yvalues= frequencies
        cleandynspec.par.xvalues= times
        cleandynspec.par.tbin = tbin

    if verbose:
        print "Finished - dyncalc time used:",time.clock()-t0,"s."

    if clean:
        return dynspec, cleandynspec
    else:
        return dynspec

def calcIndividualStationDynspecs(beams,dm=0,save_file=False,fraction=None,tbin=1,clean=False,verbose=1,time_range=None):
    '''
    Calculates the dynamic spectrum of individual beams.

    ============ ===== ===================================================================
    *beams*            Input array.
    *dm*         0     Dedispersion Measure
    *save_file*  False Saves a file in hArray format with additional information besides the array values.
    *fraction*   None  If not None, then a list of the form [x,y] such that extracting the fraction x/y of the data. with x>-1, and y>=x.
    *tbin*       1     If >1 integrates over this number of blocks. Or time binning.
    *clean*      False If True it calculates the cleaned spectrum.
    *verbose*    1     If want to print status and performance.
    *time_range* None  List with the begin and end of time selection. [t1,t2] in s
    ============ ===== ===================================================================

    Example::

        import Beam_Tools as bt
        dynspecs = bt.calcIndividualStationDynspecs(beams)

    or::

        import Beam_Tools as bt
        dynspecs,cleandynspecs = bt.calcIndividualStationDynspecs(beams,tbin=16,clean=True,save_file=True,time_range=[0,.1])

    '''
    #-----------
    #Log of things to check:
    #1)tbin divisiion
    #-----------

    t0=time.clock()

#    filename_bin = os.path.join(filename,"data.bin")   #Maybe needed in the future if using "from_file" option.

    speclen = beams['BEAM_SPECLEN']
    block_duration = beams['BLOCKSIZE']*beams['SAMPLE_INTERVAL'][0]
    nblocks = beams['NCHUNKS']*beams['BEAM_NBLOCKS']

    if not beams['DM'] and dm:
        beams['DM'] = dm

    if time_range:
        #Time selection indexing.
        if type(time_range) != type([]) or len(time_range)!=2 or time_range[0] > time_range[1] or time_range[0]<0 or time_range[1]<0:
            raise ValueError('Need a list of lenght 2, with first element "<" or "=" second element. Both elements positive.')
        times = beams.getTimes()
        block_range = cr.hArray(int, 2)
        cr.hFindSequenceBetweenOrEqual(block_range,times,time_range[0], time_range[1], 0, 0)
        nblocks = block_range[1] - block_range[0]
        fraction = False
        block_range=range(block_range[0],block_range[1])

    if not fraction:
        fraction=[1,1]
        if not time_range:
            block_range = range(0,nblocks)
    else:
        if type(fraction) != type([]) or len(fraction)!=2 or fraction[0] > fraction[1] or fraction[0]<0 or fraction[1]<0:
            raise ValueError('Need a list of lenght 2, with first element "<" or "=" second element. Both elements positive.')
        if fraction[0]==0: fraction[0]=1
        nblocks = int(nblocks/fraction[1])
        block_range=range((fraction[0]-1)*nblocks,(fraction[0])*nblocks)

    beam = beams.empty('FFT_DATA')
    dynspec = cr.hArray(float,[beam.shape()[0],nblocks,beam.shape()[1]])

    #Dynamic spectrum calculation.
    for i,block in enumerate(block_range):
        if verbose:
            print 'Dynamic spectrum calculation at {0:.2%}  \r'.format(float(i)/nblocks),
        sys.stdout.flush()
        beams.getFFTData(beam,block)
        for b in range(beam.shape()[0]):
            dynspec[b,i].spectralpower2(beam[b])

    #Time integration.
    if tbin>1:
        dynspec = cr.hArray_toNumpy(dynspec)
        dynspec = dynspec.reshape((dynspec.shape[0],dynspec.shape[1]/tbin,tbin,dynspec.shape[2]))
        dynspec = np.sum(dynspec,axis=2)
        dynspec = cr.hArray(dynspec)

    #Cleaning dynamic spectrum.
    if clean:
        cleandynspec=dynspec*0
        for b in range(beam.shape()[0]):
            avspec = cr.hArray(float,speclen,fill=0.0)
            dynspec[b,...].addto(avspec)
            cleandynspec[b] = (dynspec[b]*dynspec.shape()[1]/avspec)[0]

    #Transposing arrays.
    dynspec_T = dynspec.Transpose()*0
    for b in range(beam.shape()[0]):
        dynspec_T[b] = dynspec[b].Transpose()
    dynspec = dynspec_T

    if clean:
        cleandynspec_T = cleandynspec.Transpose()*0
        for b in range(beam.shape()[0]):
            cleandynspec_T[b] = cleandynspec[b].Transpose()
        cleandynspec = cleandynspec_T

    #Create a frequency vector
    frequencies = beams['BEAM_FREQUENCIES']

    #Create a time vector
    start_time=block_range[0]*block_duration
    end_time=block_range[-1]*block_duration
    times = cr.hArray(float,int(round((end_time-start_time)/(block_duration*tbin))),name="Time",units=("","s"))
    times.fillrange(start_time,block_duration*tbin)

    #Adding parameters
    dynspec.par.yvalues= frequencies
    dynspec.par.xvalues= times
    dynspec.par.tbin = tbin
    if clean:
        cleandynspec.par.yvalues= frequencies
        cleandynspec.par.xvalues= times
        cleandynspec.par.tbin = tbin

    #Saving file(s).
    if save_file:
        dynspec.write(os.path.join(filename,"dynspec"),nblocks=1,block=0,clearfile=True)
        print 'Saving binary in %s' %os.path.join(filename,"dynspec.pcr")
        if clean:
            cleandynspec.write(os.path.join(filename,"clean_dynspec"),nblocks=1,block=0,clearfile=True)
            print 'Saving binary in %s' %os.path.join(filename,"clean_dynspec.pcr")

    if verbose:
        print "Finished - dyncalc time used:",time.clock()-t0,"s."

    if clean:
        return dynspec, cleandynspec
    else:
        return dynspec

def showPulseProfiles(beams,freq_range=[130,160],time_range=[0,0.1],doplot=False,tbin=16):
    ''' Show pulse profile for each beam.

    ============== ===== ===================================================================
    *beams*              Input beams, opened by the beam.py interphase.
    *freq_range*   None  Frequency range in MHz
    *time_range*   None  Time range in seconds
    *tbin*         16    If >1 integrates over this number of blocks. Or time binning.
    *doplot*       False If True then it plots the flat spectra of each beam.
    ============== ===== ===================================================================

    Example::

        import Beam_Tools as bt
        dynspecs,flat_dynspecs = bt.showPulseProfiles(beams,freq_range=[175,195],time_range=[0.1,0.2])

    '''

    Obs_ID = beams['TBB_FILENAME'].split('/')[-1].split('_')[0]

    #Dynspec calc.
    print 'Calculating single beam dynamic spectra.'
    dynspec, cleandynspec = calcIndividualStationDynspecs(beams,tbin=tbin,clean=True)
    dynspec = cleandynspec

    SNR = cr.hArray(float,dynspec.shape()[0])

    if doplot:
        cr.plt.ion()

    #Cut and Integrate frequencies.
    all_flat_dynspec = cr.hArray(float,[dynspec.shape()[0],dynspec.shape()[2]])
    for dyn in range(dynspec.shape()[0]):
        zoom_dynspec = cutDynspec(dynspec[dyn],freq_range=freq_range,time_range=time_range)
        flat_dynspec = flatDynspec(zoom_dynspec,axis='x')

        stats = doStats(flat_dynspec,verbose=False)
        SNR[dyn] = stats[5]/stats[2]
        print '  '*2
        print 'The SNR in '+beams['STATION_NAME'][dyn]+'_'+beams['ANTENNA_SET'][dyn]+' is:', SNR[dyn]

        #ploting
        if doplot:
            times = zoom_dynspec.par.xvalues
            cr.plt.plot(times,flat_dynspec/flat_dynspec.vec().mean() + dyn*0.1,label= beams['STATION_NAME'][dyn]+'_'+beams['ANTENNA_SET'][dyn]+', SNR = '+str(SNR[dyn]))

        all_flat_dynspec[dyn,...] = flat_dynspec.vec()

    if doplot:
        cr.plt.rcParams['font.size']=20
        cr.plt.rcParams['axes.titlesize']=30
        cr.plt.xlabel('Relative Time [s]')
        cr.plt.ylabel('Power [Relative units]')
        cr.plt.title(Obs_ID+': Pulse profile for individual beams.')
        cr.plt.legend()

    return dynspec,all_flat_dynspec

def findRFI(dynspec,return_more=False, single_plot=False):
    '''
    Find RFI in dynspec by sum(square(v)) / square(sum(v))
    with the sum over the blocks for each frequency.

    ============== ===== ===================================================================
    *dynspec*
    *return_more*
    *plot*                 List with beam index [n],
    ============== ===== ===================================================================

    Example::

        import Beam_Tools as bt
        RFI_channels = bt.findRFI(dynspec)

    '''

    if len(dynspec.shape()) == 2:
        dynspec.reshape([1,dynspec.shape()[0],dynspec.shape()[1]])

    dynspec_sqr = cr.hArray(float,dynspec.shape(),dynspec)
    dynspec_sqr.square()

    sum_sqr = cr.hArray(float,[dynspec.shape()[0],dynspec.shape()[1]])
    sqr_sum = cr.hArray(float,[dynspec.shape()[0],dynspec.shape()[1]])

    #adding blocks
    for beam in range(dynspec.shape()[0]):
        sqr_sum[beam,...] = dynspec[beam,...].sum()
        sum_sqr[beam,...] = dynspec_sqr[beam,...].sum()

    #Squaring the sums
    sqr_sum.square()

    #Define ratio
    dynspec_ratio = sum_sqr / sqr_sum

    if dynspec.shape()[1] > 10:
        channel_cut = dynspec.shape()[1] / 10
    else:
        channel_cut = 1

    RFI_channels = []

    #Define thresholds
    for beam in range(dynspec.shape()[0]):
        scrash_vec = cr.hArray(float,dynspec_ratio.shape()[1],sorted(dynspec_ratio[beam].vec()))
        sigma = scrash_vec[channel_cut:-1*channel_cut].stddev()

        RFI_threshold1 = scrash_vec.median()+4*sigma
        RFI_threshold2 = scrash_vec.median()-4*sigma

        #Find RFI
        RFI_1 = dynspec_ratio[beam].vec().Find('>',RFI_threshold1[0])
        RFI_2 = dynspec_ratio[beam].vec().Find('<',RFI_threshold2[0])
        RFI_channels.append(sorted(list(set(RFI_1) | set(RFI_2))))

        if beam == single_plot[0]:
            cr.plt.ion()
            cr.plt.plot(dynspec_ratio[beam].vec())
            cr.plt.plot([0,dynspec.shape()[1]],[RFI_threshold1[0],RFI_threshold1[0]])
            cr.plt.plot([0,dynspec.shape()[1]],[RFI_threshold2[0],RFI_threshold2[0]])

    if return_more:
        return RFI_channels, dynspec_ratio
    else:
        return RFI_channels

def calibSpectrumGain():
    '''Performs the frequency gain calibration.
    '''
    raise NotImplementedError


def fitBaseline(flat_dynspec,order=10):
    '''This program fits the baseline.
    '''

    flat_numpy_dynspec = cr.hArray(float,flat_dynspec,flat_dynspec)
    flat_numpy_dynspec = flat_numpy_dynspec.toNumpy()
    x_numpy_values = flat_dynspec.par.xvalues.toNumpy()
    numpy_coeffs = np.polyfit(x_numpy_values,flat_numpy_dynspec,order)
    numpy_fit = np.polyval(numpy_coeffs,x_numpy_values)

    flat_fit = cr.hArray(numpy_fit)
    flat_fit.sqrt()
    flat_fit = 1/flat_fit

    return flat_fit



def dynDM(beams,dynspec,DM=0,Ref_Freq=None,from_file=False,verbose=False,save_file=False):
    """
    Do dedispersion by integer shifting.
    Calculate dedispersed time series.

    ============== ===== ===================================================================
    *beams*              Input beam
    *cleandynspec* None  Array with cleaned dynamic spectrum.
    *DM*           0     Dispersion Measure.
    *Ref_Freq*     None  Reference frequencies in Hz
    *from_file*    False Read cleandynspec from file.
    *verbose*      False If true then prints extra information, plot the dedispersed dynspec, and calculates/plots the time series.
    *save_file*    False Saves a file in hArray format with additional information besides the array values.
    ============== ===== ===================================================================

    Example::

        dedispersed_dynspec = Task.dynDM(cleandynspec=cleandynspec,DM=26.83,Ref_Freq=[151e6,170e6])

    or::

        Task.dynDM(dynspec,DM=26.83,Ref_Freq=170e6,from_file=True,verbose=True)

    Preferably use a cleaned dynamic spectrum as input.
    """

    raise NotImplementedError

    Ref_Freq = cr.asval(Ref_Freq)

    block_duration = beams['BLOCKSIZE']*beams['SAMPLE_INTERVAL'][0]
    tbin = cleandynspec.par.tbin

    #Create a frequency vector
    frequencies = cleandynspec.par.yvalues

    #Create a time vector
    times = cleandynspec.par.xvalues

    #Dedispersion parameters
    dedispersed_dynspec = cr.hArray(float,cleandynspec,fill=0)
    dt = block_duration*tbin

    #Calculate the relative shifts in samples per frequency channels
    #Constant value comes from "Handbook of Pulsar Astronomy - by Duncan Ross Lorimer , Section 4.1.1, pagina 86 (in google books)"
    shifts = ( 4.148808e-3*DM/dt) * 1e9**2 * (Ref_Freq**-2 - frequencies**-2)

    #Integer offsets to reference frequency (shift to center)
    offsets = cr.Vector(int,frequencies,fill=shifts)
#    offsets += times.shape()[0]/2

    #Now do the actual dedispersion by integer shifting ... that's all
    dedispersed_dynspec[...].shift(cleandynspec[...],offsets.vec())

    #Adding parameters
    dedispersed_dynspec.par.yvalues= frequencies
    dedispersed_dynspec.par.xvalues= times

    if save_file:
        dedispersed_dynspec.write(os.path.join(filename,"dedispersed_dynspec"),nblocks=1,block=0,clearfile=True)
        print 'Saving binary in %s' %os.path.join(filename,"dedispersed_dynspec.pcr")

    return dedispersed_dynspec

#=================REAL BAR OF DEPRECATION ========= anything below needs to go, or mayorly debugged.

def dyncalc_multibeam(filename=None,beams=None,nbeam=0,save_file=False,from_file=False,fraction=None,tbin=1,clean=False):
    '''
    Calculates the dynamic spectrum.

    =============== ===== ===================================================================
    *beams*         None  Input array.
    *nbeam*         0     Beam to work with, if ()beams has stored multiple ones.
    *save_file*     False Saves a file in hArray format with additional information besides the array values.
    *from_file*    False Read cleandynspec from file.
    *fraction*      None  If not None, then a list of the form [x,y] such that extracting the fraction x/y of the data. with x>-1, and y>=x.
    *tbin*          1     If >1 integrates over this number of blocks. Or time binning.
    *clean*         False If True it calculates the cleaned spectrum.
    =============== ===== ===================================================================

    Example::

        import Beam_Tools as bt
        dynspec = bt.dyncalc(filename)

    or::

        import Beam_Tools as bt
        dynspec,cleandynspec = bt.dyncalc(beams,tbin=16,clean=True,save_file=True)

    The regular and clean dynamic spectra (all blocks) are returned and stored in ``Task.dynspec`` and ``Task.cleandynspec`` respectively.
    '''

    raise NotImplementedError


    if beams==None and filename==None:
        raise ValueError('Need to provide either the filename or the opened .beam file')

    if nbeam!=0 or from_file or save_file:
        raise KeyError("Keyword is invalid for now: "+key)

    t0=time.clock()

    if beams==None:
        beams=cr.open(filename)
    elif filename==None:
        filename=beams['FILENAMES'][0]

#    filename_bin = os.path.join(filename,"data.bin")   #Maybe needed in the future if using "from_file" option.

    speclen = beams['BLOCKSIZE']/2+1
    block_duration = beams['BLOCKSIZE']*beams['SAMPLE_INTERVAL'][0]
    nblocks = beams['NCHUNKS']*beams['BEAM_NBLOCKS']

    if not fraction:
        fraction=[1,1]
    else:
        if type(fraction) != type([]) or len(fraction)!=2 or fraction[0] > fraction[1] or fraction[0]<0 or fraction[1]<0:
            raise ValueError('Need a list of lenght 2, with first element "<" or "=" second element. Both elements positive.')
        if fraction[0]==0: fraction[0]=1
        nblocks = int(nblocks/fraction[1])

    start_time=(fraction[0]-1)*(block_duration*nblocks)/fraction[1]
    end_time=fraction[0]*(block_duration*nblocks)/fraction[1]

    beam = beams.empty('FFT_DATA')
    tm = cr.hArray(float,beam)
    dynspec = cr.hArray(float,[nblocks,beam.shape()[0],speclen])
    block_range=range((fraction[0]-1)*nblocks,(fraction[0])*nblocks)
    pdb.set_trace()
    #Reading beams.
    for block in block_range:
        print ' Calculation at {0:.2%}  \r'.format(float(block)/len(block_range)),
        sys.stdout.flush()
        beams.getFFTData(beam,block)
        tm[...].spectralpower2(beam[...])
        dynspec[block] = tm

    #Time integration.
    if tbin>1:
        dynspec = cr.hArray_toNumpy(dynspec)
        dynspec = dynspec.reshape((dynspec.shape[0]/tbin,tbin,dynspec.shape[1],dynspec.shape[2]))
        dynspec = np.sum(dynspec,axis=1)
        dynspec = cr.hArray(dynspec)

    #Rearranging dimensions.
        dynspec = cr.hArray_toNumpy(dynspec)
#        dynspec = np.rollaxis(dynspec,0,3)
        dynspec = np.swapaxes(dynspec,0,1)
        dynspec = np.swapaxes(dynspec,1,2)
        dynspec = cr.hArray(dynspec.copy())  #The .copy() is quick&dirty way to avoid a bug.

    #Cleaning dynamic spectrum.
    if clean:
        avspec=cr.hArray(float,speclen,fill=0.0)
        cleandynspec=cr.hArray(float,dynspec,fill=0.0)
        for dim in range(dynspec.shape()[0]):
            single_dynspec = cr.hArray(float,dynspec.shape()[1:],dynspec[dim])
            single_dynspec[...].addto(avspec)
            cleandynspec[dim] = single_dynspec/avspec

    #Create a frequency vector
    frequencies = beams['BEAM_FREQUENCIES']

    #Create a time vector
    times = cr.hArray(float,int(round((end_time-start_time)/(block_duration*tbin))),name="Time",units=("","s"))
    times.fillrange(start_time,block_duration*tbin)

    #Adding parameters
    dynspec.par.yvalues= frequencies
    dynspec.par.xvalues= times
    dynspec.par.tbin = tbin
    if clean:
        cleandynspec.par.yvalues= frequencies
        cleandynspec.par.xvalues= times
        cleandynspec.par.tbin = tbin

    #Saving file(s).
    if save_file:
        dynspec.write(os.path.join(filename,"dynspec"),nblocks=1,block=0,clearfile=True)
        print 'Saving binary in %s' %os.path.join(filename,"dynspec.pcr")
        if clean:
            cleandynspec.write(os.path.join(filename,"clean_dynspec"),nblocks=1,block=0,clearfile=True)
            print 'Saving binary in %s' %os.path.join(filename,"clean_dynspec.pcr")

    print "Finished - total time used:",time.clock()-t0,"s."

    if clean:
        return dynspec, cleandynspec
    else:
        return dynspec

def cutbeam(beam,axis='xy',freq_range=[0,0],startblock=0,nblocks=1):
    '''Cuts a beam around the dedispersed peak.
    '''

    raise NotImplementedError


    if not np.any(freq_range):
        freq_range=[100,200]

    if not nblocks:
        if axis=='xy':
            nblocks=640
        else:
            nblocks=beam.shape()[0]

    speclen = 8193

    freq_range = cr.hArray(float,[2],fill=freq_range,units=('M','Hz')).setUnit("","Hz")
    frequencies = cr.hArray(float,speclen, list(np.arange(1e+08,2e+08+12207.03125,12207.03125)))
    frequency_slice = cr.hArray(int, 2)
    cr.hFindSequenceBetweenOrEqual(frequency_slice,frequencies,freq_range[0], freq_range[1], 0, 0)

    if axis=='xy':
#        sliced_beam = cr.hArray(complex,[nblocks,beam.shape()[1]],beam)
        sliced_beam = cr.hArray(float,[nblocks,beam.shape()[1]],beam[startblock:nblocks+startblock])
    else:
        sliced_beam = beam

    T_beam = sliced_beam.Transpose()

    zoom_beam = cr.hArray(complex,[frequency_slice[1]-frequency_slice[0],nblocks],T_beam[frequency_slice[0]:frequency_slice[1]])

    zoom_beam=zoom_beam.Transpose()

    times = cr.hArray(float,nblocks,list(np.arange(startblock*8192*2*5e-9,(startblock+nblocks)*8192*2*5e-9,8192*2*5e-9)))
    frequencies = cr.hArray(float,frequency_slice[1]-frequency_slice[0],frequencies[frequency_slice[0]:frequency_slice[1]].vec())

    zoom_beam.par.xvalues = times
    zoom_beam.par.yvalues = frequencies

    return zoom_beam

def getCAL_DELAY(stations,antenna_set):
    '''Will return a delay(s) for the corresponding station(s).
    '''

    raise NotImplementedError


    CAL_DELAY=[0.0,1.6859375e-09,2.3625e-09,-5.8421875e-09,2.6046875e-09,-1.3371875e-08]
    CAL_DELAY=[0.0,8.00781249999e-10,-4.98046875e-09,-8.046875e-09,-5.05859375e-09,-4.21875e-09,3.12499999998e-10,7.81250000002e-10,-5.4296875e-09,-8.41796875e-09 ,-5.1953125e-09 ,-5.3515625e-09]


def dynplot(dynspec,log=True):
    '''Plot dynamic spectrum.
    '''

    raise NotImplementedError



    #Create a frequency vector
    frequencies = dynspec.par.yvalues

    #Create a time vector
    times = dynspec.par.xvalues

    #Plotings....

    cr.plt.ion()
    cr.plt.clf()
    if not log:
        cr.plt.imshow(cr.hArray_toNumpy(dynspec),aspect='auto',origin='lower',cmap=cr.plt.cm.hot,vmin=1e-3,vmax=0.03,extent=(times[0],times[-1],frequencies[0],frequencies[-1]))
    else:
        cr.plt.imshow(np.log10(cr.hArray_toNumpy(dynspec)),aspect='auto',origin='lower',vmin=-3,vmax=-1.5,cmap=cr.plt.cm.hot,extent=(times[0],times[-1],frequencies[0],frequencies[-1]))
    cr.plt.xlabel("+/- Time [s]")
    cr.plt.ylabel("Frequency [MHz]")

def addBeams_old_idea():
    """Returns the sum of beams in comple (All Chunks).

    Required Arguments:

    ============= =================================================
    Parameter     Description
    ============= =================================================
    *data*        data array to write FFT data to.
    *chunk*       index of chunk to return data from.
    ============= =================================================

    Output:


    """
    raise NotImplementedError

    #Reading beams.
    for file in filename_bin:
        print file
    for chunk in chunk_range:
        print ' Calculation at {0:.2%}  \r'.format(float(chunk)/len(chunk_range)),
        sys.stdout.flush()
        for nb in range(nbeams):
#                    print filename_bin[nb]
            bm.readfilebinary(filename_bin[nb],chunk*speclen*nbeams*nblocks)
            bm/=float(nb+1.)
            if nb>1:
                beams*=(nb)/(nb+1.)
            beams+=bm
            bm*=0
