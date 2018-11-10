#------------------------------------------------------------------------------------------------
# Testing the wavelenght gain calibration
#------------------------------------------------------------------------------------------------

from pycrtools import *
import Beam_Tools as bt
from pycrtools.tasks import fitbaseline, findrfi
import matplotlib.pyplot as plt
import socket

cr= 0
plt.ion()

#----------------------------------
if cr:
    filename_lofar = LOFARSOFT + "/data/lofar/VHECR_example.h5"

    datafile = open(filename_lofar)

    fftdata = datafile.empty("FFT_DATA")
    nblocks = datafile["MAXIMUM_READ_LENGTH"] / datafile["BLOCKSIZE"]
    avspectrum = hArray(float, dimensions=fftdata, name="Average spectrum")

    # Calculate spectral power
    for block in range(nblocks):
        datafile["BLOCK"] = block
        fftdata.read(datafile, "FFT_DATA")
        hSpectralPower(avspectrum[...], fftdata[...])

    frequencies = datafile["FREQUENCY_DATA"].setUnit("M", "")

    # Plotting
    plt.figure()
    avspectrum.par.xvalues = frequencies
    avspectrum.par.title = "Average spectrum"
    avspectrum[0].plot(logplot="y", clf=False)

#----------------------------------
#Baseline fiting
    flat_dynspec = hArray(float,[513],avspectrum[0].vec())
    flat_dynspec.par = avspectrum.par

    fitbaseline=trerun("FitBaseline","",spectrum=flat_dynspec,fittype='POLY',numin=10.,numax=90.,doplot=3,ncoeffs=12)
    #calcbaseline=trerun("CalcBaseline",fitbaseline,spectrum=fitbaseline.spectrum,normalize=True,doplot=3)
    calcbaseline=trerun("CalcBaseline","",fitbaseline.spectrum,doplot=3,normalize=False)
    applybaseline=trun("ApplyBaseline",calcbaseline.spectrum,doplot=3)

#----------------------------------
#Applyting bandpass in complex array, test.

    avspectrum_complex = hArray(float, dimensions=fftdata, name="Average spectrum")

    # Calculate spectral power
    for block in range(nblocks):
        datafile["BLOCK"] = block
        fftdata.read(datafile, "FFT_DATA")
        fftdata*= calcbaseline.baseline
        hSpectralPower(avspectrum_complex[...], fftdata[...])

    # Plotting
    plt.figure()
    avspectrum.par.xvalues = frequencies
    avspectrum.par.title = "Average spectrum, comparison"
    plt.plot(calcbaseline.spectrum.vec(),label='Before')
    plt.plot(avspectrum_complex[0].vec(),label='After')
    plt.legend()

else:
#----------------------------------
#Testing with FRATS data (HBA)
    node_local = socket.gethostname()
    if node_local != 'locus013':
        filename = '/vol/astro1/lofar/frats/tbb/data/L43784_D20120125T2111/L43784_D20120125T211154.887Z_CS002_R000_tbb.h5'
        nobl=200
    else:
#        filename = '/data/TBB/FRATS/tbb/data/L211994_D20140316T1454/L211994_D20140316T145444.000Z_CS002_R000_tbb.h5'
#        nobl=10
        filename = '/data/TBB/FRATS/tbb/data/L199858_D20140119T1641/L199858_D20140119T164157.000Z_CS002_R002_tbb.h5'
        nobl=200

    file=open(filename)
    file['BLOCKSIZE'] = 2**14
    file['SELECTED_DIPOLES'] = 'even'

    sbl= 0
    freq_range = (100, 200)
    clean = 0

    rfi = trun("FindRFI", f=file,nofblocks=nobl,startblock=sbl,freq_range=freq_range,save_plots=False,verbose=False)

    rfi.median_cleaned_spectrum.par.xvalues=rfi.f['FREQUENCY_DATA'].setUnit('M','Hz')
    rfi.median_average_spectrum.par.xvalues=rfi.f['FREQUENCY_DATA'].setUnit('M','Hz')

    plt.figure()
    plt.plot(np.log10(rfi.median_average_spectrum.toNumpy()))

#----------------------------------
#Baseline fiting
    if clean:
        fitbaseline_frat_med=trerun("FitBaseline","",spectrum=rfi.median_cleaned_spectrum,fittype='POLY',numin=110.,numax=190.,doplot=3)#,ncoeffs=18)
        calcbaseline_frat_med=trerun("CalcBaseline","",fitbaseline_frat_med.spectrum,doplot=3,normalize=False)
        applybaseline_frat_med=trun("ApplyBaseline",calcbaseline_frat_med.spectrum,doplot=3)

    else:
        fitbaseline_frat_avg=trerun("FitBaseline","",spectrum=rfi.median_average_spectrum,fittype='POLY',numin=110.,numax=190.,doplot=3)#,ncoeffs=18)
        calcbaseline_frat_avg=trerun("CalcBaseline","",fitbaseline_frat_avg.spectrum,doplot=3,normalize=False)
        plt.figure()
        plt.plot(calcbaseline_frat_avg.spectrum*calcbaseline_frat_avg.baseline)
        plt.title('calcbaseline_frat_avg.spectrum*calcbaseline_frat_avg.baseline')
        plt.figure()
        applybaseline_frat_avg=trun("ApplyBaseline",calcbaseline_frat_avg.spectrum,doplot=3)



















'''
From cr_event.py
FitBaseline = dict(ncoeffs=80,numin=30,numax=85,fittype="POLY",splineorder=3)
ApplyBaseline=dict(rmsfactor=7)


fitbaseline=trerun("FitBaseline","",averagespectrum_good_antennas,extendfit=0.5,pardict=par,doplot=3 if Pause.doplot else 0)
CalcBaseline = dict(baseline=False), # Make sure baseline is recreated when the task is run a second time
print "---> Calculate a smooth version of the spectrum which is later used to set amplitudes."
calcbaseline1=trerun("CalcBaseline",1,averagespectrum_good_antennas,pardict=par,invert=False,HanningUp=False,normalize=False,doplot=0)
amplitudes=hArray(copy=calcbaseline1.baseline)
print "---> Calculate baseline again, but now for multiplication with data to flatten the spectrum."
calcbaseline_flat=trerun("CalcBaseline","flat",averagespectrum_good_antennas,pardict=par,invert=True,normalize=False,doplot=Pause.doplot)
print "---> Calculate a baseline with Galactic powerlaw"
calcbaseline_galactic=trerun("CalcBaseline","galactic",averagespectrum_good_antennas,pardict=par,invert=True,normalize=False,powerlaw=0.5,doplot=Pause.doplot)

applybaseline=trerun("ApplyBaseline","",station_spectrum,baseline=station_gaincurve,pardict=par,doplot=Pause.doplot)
fft_data[...].randomizephase(applybaseline.dirty_channels[...,[0]:applybaseline.ndirty_channels.vec()],amplitudes[...])
'''

'''
TAB,dynspec,cleandynspec = bt.addBeams(beams,dyncalc=True,clean=True)
zoom_dynspec = bt.cutDynspec(dynspec,time_range=[0.0,0.2])
flat_dynspec = bt.flatDynspec(zoom_dynspec,axis='y',verbose=True)
#---
dynspecs = bt.calcSingleStationDynspec(beams,time_range=[0,.2])
dynspec0=hArray(float,[8193,2442],dynspecs[0])
dynspec1=hArray(float,[8193,2442],dynspecs[1])
dynspec0.par  = dynspecs.par
dynspec1.par  = dynspecs.par
flat_dynspec0 = bt.flatDynspec(dynspec0,axis='y',verbose=True)
plt.figure()
flat_dynspec1 = bt.flatDynspec(dynspec1,axis='y',verbose=True)
flat_dynspecs = hArray(float,[2,8193])
flat_dynspecs[0] = flat_dynspec0.vec()
flat_dynspecs[1] = flat_dynspec1.vec()
flat_dynspecs.par = flat_dynspec0.par
#---

fitbaseline=trerun("FitBaseline","",spectrum=flat_dynspec,doplot=3,fittype='POLY',extendfit=0.5)
fitbaseline=trerun("FitBaseline","",spectrum=flat_dynspec,extendfit=0.5,doplot=1)
fitbaseline=trerun("FitBaseline","",spectrum=flat_dynspec,extendfit=0.5,doplot=0)

fitbaseline=trerun("FitBaseline","",spectrum=flat_dynspec,fittype='POLY',extendfit=0.5,doplot=3)
fitbaseline3=trerun("FitBaseline","",spectrum=flat_dynspec,fittype='POLY',extendfit=0.5,doplot=3,ncoeffs=2)

calcbaseline=trerun("CalcBaseline",'flat',flat_dynspec,invert=False,HanningUp=False,normalize=False,doplot=3,addHanning=False)
calcbaseline_flat=trerun("CalcBaseline","flat",flat_dynspec,invert=True,normalize=False)
calcbaseline1=trerun("CalcBaseline",1,flat_dynspec,invert=False,HanningUp=False,normalize=False,doplot=1)


'''
