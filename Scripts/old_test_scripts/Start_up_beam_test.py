#! /usr/bin/env python
'''Start up a set of beams.
'''


from pycrtools import *
from pytmf import *
import Beam_Tools as bt
import glob
import os
import pdb;# pdb.set_trace()

Obs_ID ='L43784'
polarization=0

FRATS_ANALYSIS=os.environ["FRATS_ANALYSIS"].rstrip('/')+'/'

beam_suffix='*pol%i_*HBA?.LOFAR_centered.beam'%(polarization)
beam_suffix='*CS00*pol%i_*HBA?.LOFAR_centered.beam'%(polarization)
filenames=glob.glob(FRATS_ANALYSIS+Obs_ID+'/beam.results/'+beam_suffix)

print 'Using files : ',filenames


beams = open(filenames)



stop
################################################################################################################


#beam_suffix='*CS00[2-7]*pol0.sample*'
#beam_suffix='*CS002*pol0_sample*'


from pycrtools import *
import Beam_Tools as bt
beam_suffix='*CS??[2-7]*pol0*HBA?.LOFAR_centered.beam'
filenames=glob.glob(beam_suffix)
beams=open(sorted(filenames))

#------
beams['DM']=26.76
ST_DELAYS = bt.ccBeams(beams,freq_range=[151,168],time_range=[.021,.0255],verbose=1)
beams['CAL_DELAY']=ST_DELAYS
#------
TAB,dynspec,cleandynspec = bt.addBeams(beams,dyncalc=True,tbin=16,dm=26.76,clean=True)
#------
TAB=hArray(complex,[191*64,8193])
TAB.readfilebinary('L43784_D20120125T211154.Superterp.pol0.TAB.bin',0)
dynspec,cleandynspec = bt.rawdyncalc(TAB,tbin=16,clean=True)
#------
time_range=[0, 0.05]
freq_range=[152, 167]
zoom_dynspec = bt.cutDynspec(cleandynspec,freq_range=freq_range,time_range=time_range)
flat_cleandyn=bt.flatDynspec(zoom_dynspec,axis='x')
stats = bt.doStats(flat_cleandyn,verbose=True,nof_St=beams['NOF_BEAM_DATASETS'])

#----------
#For L74100
time_range=[0, 0.05]
freq_range=[175,195]
zoom_dynspec_CS002 = bt.cutDynspec(cleandynspec_CS002,freq_range=freq_range,time_range=time_range)
flat_cleandyn_CS002=bt.flatDynspec(zoom_dynspec_CS002,axis='x')
stats_CS002 = bt.doStats(flat_cleandyn_CS002,verbose=True,nof_St=str(beams['NOF_BEAM_DATASETS']))
#------------------------------------------------------------------------------------------------
times = dynspec.par.xvalues
freqs = dynspec.par.yvalues
extent=(times[0],times[-1],freqs[0],freqs[-1])
plt.imshow(np.log10(hArray_toNumpy(cleandynspec)),aspect='auto',origin='lower',cmap=plt.cm.hot,extent=extent,vmax=-2.5,vmin=-3.5)
plt.imshow(np.log10(hArray_toNumpy(dynspec)),aspect='auto',origin='lower',cmap=plt.cm.hot,extent=extent,vmax=7,vmin=3)
#------------------------------------------------------------------------------------------------


beam_suffix='*CS???*pol0*HBA?.LOFAR_centered.beam'
filenames = glob.glob(beam_suffix)
beams=open(sorted(filenames))

cleandynspec_dm=hArray(float,[380,8193],fill=0)
cleandynspec_dm.readfilebinary('L198910_D20140115T1226.cleandynspec_dm.core.bin')
cleandynspec_dm.par.yvalues = beams.getFrequencies()
cleandynspec_dm.par.xvalues=beams.getTimes()
cleandynspec_dm.par.tbin = 16

plt.imshow(np.log10(hArray_toNumpy(cleandynspec_dm)),aspect='auto',origin='lower',cmap=plt.cm.hot,extent=extent,vmax=-2.,vmin=-3.5)

dynspec_dm_140_150 = bt.cutDynspec(cleandynspec_dm,freq_range=[140,150])
flat_dynspec_140_150 = bt.flatDynspec(dynspec_dm_140_150,axis='x',verbose=True)

#------------------------------------------------------------------------------------------------
#Jupiter
from pycrtools import *
beam_suffix='*pol0*TauA*'
beam_suffix='*pol0*CasA*'
filenames=glob.glob(beam_suffix)
beams=open(sorted(filenames))
beamsTAB,dynspec,cleandynspec = bt.addBeams(beams,dyncalc=True,tbin=16,clean=True,incoherent=True)
beams
import Beam_Tools as bt
TAB,dynspec,cleandynspec = bt.addBeams(beams,dyncalc=True,tbin=16,clean=True,incoherent=True)
ST_DELAYS = bt.ccBeams(beams,freq_range=[41,80],time_range=[0.6,1.1],verbose=1)



#------------------------------------------------------------------------------------------------
# Wavelenght gain calibration

from pycrtools import *
import Beam_Tools as bt
from pycrtools.tasks import fitbaseline
beam_suffix='*CS002*pol0_sample*'
filenames=glob.glob(beam_suffix)
beams=open(filenames)
#---
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

#------------------------------------------------------------------------------------------------
#---
TAB_pol0_HBA0,dynspec_pol0_HBA0,cleandynspec_pol0_HBA0 = bt.addBeams(beams_pol0_HBA0,dyncalc=True,clean=True)
zoom_dynspec = bt.cutDynspec(dynspec,time_range=[0.0,0.2])
flat_dynspec = bt.flatDynspec(zoom_dynspec,axis='y',verbose=True)
#---
TAB_pol0_HBA1,dynspec_pol0_HBA1,cleandynspec_pol0_HBA1 = bt.addBeams(beams_pol0_HBA1,dyncalc=True,clean=True)
zoom_dynspec = bt.cutDynspec(dynspec,time_range=[0.0,0.2])
flat_dynspec = bt.flatDynspec(zoom_dynspec,axis='y',verbose=True)
#---
TAB_pol1_HBA0,dynspec_pol1_HBA0,cleandynspec_pol1_HBA0 = bt.addBeams(beams_pol1_HBA0,dyncalc=True,clean=True)
zoom_dynspec = bt.cutDynspec(dynspec,time_range=[0.0,0.2])
flat_dynspec = bt.flatDynspec(zoom_dynspec,axis='y',verbose=True)
#---
TAB_pol1_HBA1,dynspec_pol1_HBA1,cleandynspec_pol1_HBA1 = bt.addBeams(beams_pol1_HBA1,dyncalc=True,clean=True)
zoom_dynspec = bt.cutDynspec(dynspec,time_range=[0.0,0.2])
flat_dynspec = bt.flatDynspec(zoom_dynspec,axis='y',verbose=True)




#------------------------------------------------------------------------------------------------
# findRFI
from pycrtools.tasks import findrfi

rfi = cr.trun("FindRFI", f=file,nofblocks=nobl,verbose=False,startblock=sbl,freq_range=freq_range,save_plots=True,plot_prefix=outdir+'/RFI/FindRFI.'+str(Obs_ID2)+'.st_'+file['STATION_NAME'][0]+'.sbl_'+str(sbl)+'.nobl_'+str(nobl)+'.pol_'+str(polarization)+'.')
rfi=trun("FindRFI",f=tbb,nofblocks=2,verbose=True,startblock=0,sigma=1e20,save_plots=False)

#------------------------------------------------------------------------------------------------
#correlator (Jana)

python $PROS/CorrelatorTBB2.py $FRATS_DATA/L103773_ev1/L103773_D20130317T104515.341Z_CS003_R000_tbb.h5 -o Test_correlation_as_given -f100 -e 200 -n 1000 -b 16

ipython
import os
PROS=os.environ["PROS"].rstrip('/')+'/'
run PROS/CalibrationTBB_public Test_correlation_on_pulsar -s 141MHz -l l

PROS/CalibrationTBB_public Solar_flare_orrelation -s 141MHz -l l

#------------------------------------------------------------------------------------------------
#beam pipeline
python $PROS/frats_event.py $FRATS_DATA/L43784_ev1/L43784_D20120125T211154.867Z_CS003_R000_tbb.h5


#------------------------------------------------------------------------------------------------
#One antenna beam -  test
filename=['L43784_D20120125T211154.867Z_CS003_R000_tbb.h5']
pointings = [{'az': 0, 'el': 90}]
from pycrtools.tasks import beamformer
from pycrtools import *
trun('BeamFormer',filenames=filename,blocklen=2**14,pointings=pointings,detail_name='Beam_test_1_antenna',NyquistZone=2,antenna_list=[0])

#------------------------------------------------------------------------------------------------
#pipeline

python frats_event.py /data/TBB/FRATS/tbb/data/L183000_D20131027T1607/L183000_D20131027T160747.488Z_CS007_R005_tbb.h5 -r 1  -m lotaas -p 1

#------------------------------------------------------------------------------------------------
#Imaging
from pycrtools import *
from pycrtools.tasks import imager,beamformer
tbb=open('L103773_D20130317T104515.341Z_CS003_R000_tbb.h5')
tbb['BLOCKSIZE']=2**14
tbb['SELECTED_DIPOLES']='even'
tbb['SELECTED_DIPOLES']=range(0,47,2)
trun("Imager",data=tbb,intgrfreq=True,nblocks=200,ntimesteps=1,phase_calibrate=True,startblock=6000,NAXIS1=91,NAXIS2=91,CDELT1=-2.,CDELT2=2.,FREQMIN=1.34e8,FREQMAX=1.36e8,output='L103773_peculiar_all_sky')
trun("Imager",data=tbb,intgrfreq=True,nblocks=200,ntimesteps=1,phase_calibrate=True,startblock=4000,NAXIS1=181,NAXIS2=181,CDELT1=-1.,CDELT2=1.,FREQMIN=1.34e8,FREQMAX=1.36e8,output='L103773_peculiar_all_sky')
trun("Imager",data=tbb,intgrfreq=True,nblocks=200,ntimesteps=1,phase_calibrate=True,startblock=6000,NAXIS1=161,NAXIS2=161,CDELT1=-1.,CDELT2=1.,FREQMIN=1.35e8,FREQMAX=1.36e8,output='L103773_peculiar_all_sky')
trun("Imager",data=tbb,intgrfreq=True,nblocks=300,ntimesteps=1,phase_calibrate=True,startblock=10000,NAXIS1=161,NAXIS2=161,CDELT1=-1.,CDELT2=1.,FREQMIN=1.35e8,FREQMAX=1.36e8,output='L103773_peculiar_all_sky')

trun("Imager",data=tbb,intgrfreq=True,nblocks=200,ntimesteps=1,phase_calibrate=True,startblock=1000,NAXIS1=161,NAXIS2=161,CDELT1=-1.,CDELT2=1.,FREQMIN=1.12e8,FREQMAX=1.13e8,output='L103773_peculiar_all_sky_ear1_pol0_112-113')
