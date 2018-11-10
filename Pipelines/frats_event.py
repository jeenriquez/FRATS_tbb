#! /usr/bin/env python
"""
**Usage:**
execfile(PYP+"pipelines/frats_event.py")

This is the basic LOFAR FRATS event processing script, including RFI excision and beamforming.
To be added: gain calibration, pulse finding?

**Parameters:**
======================= =================================================================
*fname*                 filename of raw data to be processed.

*polarization*          either 0 or 1 for selecting even or odd antennas

*inner_tiles*           If True then it uses the inner 4 tiles of a HBA station, then substation is True.

*substation*            If True, then using HBA substations (HBA0 and HBA1)

*station_centered*      True if need phase center at (sub)station center, or False at LOFAR center.

*blocklen*              Block lenght, number of samples in powers of two.

*flag_antenna*          List of antenna numbers to flag.

*rfi_find*              If 1, then performing findrfi task only; 0 for beamform only; 2 for both.

*default_input*         Extra name for input file with different DM, pointing, etc..

*ref_station*           File name of reference station, the one which will be used for time sample-calibration. Usually the one with the latest time stamp.

*ref_time*              Reference time, at which the reference station started recording, this is used in time sample-calibration.

======================= =================================================================

For more help on more parameters run ``frats_event.py --help``.

Uses beamformer in simple-calibrated TBB data.

Example::

python frats_event.py filename_tbb.h5

Note:
It needs environment variables: $BF_PARSETS, $FRATS_ANALYSIS, and $FRATS_TRIGGER

    In Coma:
        $BF_PARSETS   = /vol/astro1/lofar/frats/bf/parsets
        $FRATS_ANALYSIS = /vol/astro1/lofar/frats/tbb/analysis
        $FRATS_TRIGGER = /vol/astro1/lofar/frats/tbb/analysis/trigger_info/

    In locus013:
        $BF_PARSETS   = /globalhome/lofarsystem/log/
        $FRATS_ANALYSIS = /data/TBB/FRATS/tbb/analysis/
        $FRATS_TRIGGER = /staging1/frats/tbb/analysis/trigger_info/

.. moduleauthor:: J. Emilio Enriquez <e.enriquez 'at' astro.ru.nl>

Revision History:
V1.0 created by E. Enriquez, Feb 2013
V1.1 modified by E. Enriquez, Aug 2013
V1.2 modified by E. Enriquez, Oct 2014
"""

import matplotlib
matplotlib.use("Agg")

from optparse import OptionParser
import pycrtools as cr
from pytmf import *
import numpy as np
from pycrtools import tools
from pycrtools import bfdata as bf
from pycrtools import metadata as md
from pycrtools.tasks import beamformer, findrfi, fitbaseline
import os; import sys; import glob
import matplotlib.colors as colors
import matplotlib.cm as cmx
import pdb;# pdb.set_trace()

#------------------------------------------------------

def getRADEC(allpar,beam):
    keyRA="Observation.Beam["+str(beam)+"].angle1"
    keyDEC="Observation.Beam["+str(beam)+"].angle2"
    return float(allpar[keyRA]),float(allpar[keyDEC])

def make_list(option, opt_str, value, parser):
    setattr(parser.values, option.dest, value.replace('[','').replace(']','').split(','))

def make_list_int(option, opt_str, value, parser):
    setattr(parser.values, option.dest, map(float,value.replace('[','').replace(']','').split(',')))

def get_sample_offset(fullparsetname,file,antenna_set,max_sample_number,sample_calibrated=False,ref_station='',ref_time=None, ):
    #----------------------------------------------------
    #Calibrate station times.

    f_clock_offset = float(md.getClockCorrectionParset(fullparsetname, file['STATION_NAME'][0], antennaset=antenna_set))
    blocklen = file['BLOCKSIZE']

    if not sample_calibrated:
        #First step: Rounding the sample_number used to the nearest whole block since second start, including CLOCK_OFFSET.
        sample_offset = cr.hArray(int,1,max_sample_number)
        cr.hModulus(sample_offset, blocklen)

        sample_offset = cr.asval(blocklen - sample_offset + int(f_clock_offset/file['SAMPLE_INTERVAL'][0]))
    else:
        if not ref_time:
            if not ref_station:
                f0=file
            else:
                f0=cr.open(ref_station)

            print 'WARNING, using a ref_station assumes the dumps were within the same second, if not then give a ref_time.'
            f0_clock_offset = float(md.getClockCorrectionParset(fullparsetname, f0['STATION_NAME'][0], antennaset='HBA0'))
            t0 = max(f0['SAMPLE_NUMBER'])*f0['SAMPLE_INTERVAL'][0]+f0_clock_offset
        else:
            t0 = ref_time[1]

        t = max_sample_number+f_clock_offset
        sample_offset = [int((t0-t)/file['SAMPLE_INTERVAL'][0]), (t0-t)/file['SAMPLE_INTERVAL'][0] - int((t0-t)/file['SAMPLE_INTERVAL'][0])]

    return sample_offset

def getRADEC_2_AZEL(alpha, delta, utc, angl_offset=0.0):
    ''' Quick funciton to change coordinates from RADEC to AZEL.
    '''

    phi = deg2rad(52.915122495) #(LOFAR Superterp)
    L = deg2rad(6.869837540)

    angl_offset = deg2rad(angl_offset)
    alpha += angl_offset
    #delta = delta+angl_offset

    #----------------------------------------------------
    # Make input and output arrays for conversion
    equatorial = cr.hArray([alpha, delta])
    horizontal = equatorial.new()

    #----------------------------------------------------
    # Convert all coordinates in the input array
    # assumes difference between UTC and UT is 0 (hence ut1_utc=0.)
    cr.hEquatorial2Horizontal(horizontal, equatorial, utc, 0., L, phi)

    if horizontal[0] < 0:
        horizontal[0] += 2. * np.pi # Need to check the definitions used for positive azimut angles.

    azel = [horizontal[0],horizontal[1]]
    pointings = [{'az': azel[0], 'el': azel[1]}]

    return pointings

#------------------------------------------------------
#Command line options
#------------------------------------------------------
parser = OptionParser()

parser.add_option("-u","--substation",action="store_false",default=True,help="If True, then using HBA substations (HBA0 and HBA1).")
parser.add_option("-b","--blocklen",type="int",default=2**14,help="Block lenght")
parser.add_option("-i","--inner_tiles",action="store_true",default=False,help="If True then it uses the inner 4 tiles of a HBA station, then substation is True.")
parser.add_option("-s","--station_centered",action="store_true",default=False,help="True if need phase center at (sub)station center, or False at LOFAR center.")
parser.add_option("-p","--polarization",type="int",default=0,help="Polarization, 0 for even, 1 for odd antennas.")
parser.add_option("-f","--flag_antenna",type="str",action='callback',default=None,callback=make_list,help="List of antenna numbers to flag.")
parser.add_option("-r","--rfi_find",type="int",default=2,help="If 1, then performing findrfi task only; 0 for beamform only; 2 for both.")
parser.add_option("--default_input",type="str",default='trigger',help='Extra name for input file with different DM, pointing, etc..')
parser.add_option("--ref_station",type="str",default='',help='File name of reference station, the one which will be used for time sample-calibration. Usually the one with the latest time stamp.')
parser.add_option("--ref_time",type="str",default=None,action='callback',callback=make_list_int,help='Reference time, at which the reference station started recording, this is used in time sample-calibration.')
parser.add_option("--parset_time",action="store_false",default=True,help="Using this parameter to either use the clock_offsets from the parset files or the static values.")
parser.add_option("-m","--obs_mode",type="str",default='',help='Different observing modes, e.g. "msec" in case of milisecond file for RFI testing, "lotaas" to have the correct frequency range during lotaas observing.')

(options, args) = parser.parse_args()

if not parser.get_prog_name()=="frats_event.py":
    #   Program was run from within python
    substation = 1
    inner_tiles = 0
    polarization = 0
    station_centered = False
    blocklen = 2**14
    rfi_find = 2
    default_input = 'trigger'
    ref_station = ''
    ref_time = None
    obs_mode = ''
    parset_time = True
else:
    substation = options.substation
    inner_tiles = options.inner_tiles
    polarization = options.polarization
    station_centered = options.station_centered
    blocklen = options.blocklen
    flag_antenna = options.flag_antenna
    rfi_find = options.rfi_find
    default_input = options.default_input
    ref_station = options.ref_station
    ref_time = options.ref_time
    obs_mode = options.obs_mode
    parset_time = options.parset_time


#----------------------------------------------------
#Open given TBB file.
fname = args[0]
file = cr.open(fname)

#----------------------------------------------------
#General parameters logistics.

FRATS_ANALYSIS=os.environ["FRATS_ANALYSIS"].rstrip('/')+'/'
Obs_ID = fname.split('/')[-1].split('_')[0]
Obs_ID2 = fname.split('/')[-1].split('.')[0][:-2] # Minute accuracy.

outdir = FRATS_ANALYSIS+Obs_ID2

# Create output directories, if not already present.
if not os.path.isdir(outdir) or not os.path.isdir(outdir+'/beam.results/'):
    os.makedirs(outdir+'/beam.results/')
if not os.path.isdir(outdir+'/RFI/'):
    os.mkdir(outdir+'/RFI/')

#Naming
if station_centered:
    detail_name=''
    detail_name2=''
else:
    detail_name='.LOFAR_centered'
    detail_name2='.LOFAR_centered'

if 'msec' in obs_mode:
    rfi_find = 1

#Maximum sample number calculation (NOTE: units are in seconds!!)
#Done here to avoid difference value from specific selection (ie. polarization)
#Also takes into account dipoles with different second.
if ref_time:
    time_tbb = cr.hArray(float,len(file['SAMPLE_NUMBER']),fill=file['TIME'])
    time_tbb -= ref_time[0]
    sample_tbb = cr.hArray(float,len(file['SAMPLE_NUMBER']),fill=file['SAMPLE_NUMBER'])
    sample_tbb *= file['SAMPLE_INTERVAL'][0]
    sample_tbb += time_tbb
    max_sample_number = max(sample_tbb)
else:
    max_sample_number = max(file['SAMPLE_NUMBER']) *file['SAMPLE_INTERVAL'][0]

if 'LBA' in file['ANTENNA_SET']:
    inner_tiles = 0
    substation = 0

if inner_tiles:
    substation = 1
    RCU_selection = [17,19,29,31]
    RCU_selection2 = list(cr.hArray(int,(4),RCU_selection)+48)
    file['SELECTED_DIPOLES']=RCU_selection+RCU_selection2
else:
    poli=['even','odd']
    file['SELECTED_DIPOLES']=poli[polarization]

if flag_antenna:
    if flag_antenna[0] == '-1': #To flag the last one (it is a relativelly usual problem that this one has not all the data recorded.)
        file['SELECTED_DIPOLES'] = file['SELECTED_DIPOLES_INDEX'][:-1]
    else:
        antenna_list = np.array(file['SELECTED_DIPOLES_INDEX'])
        for flag in flag_antenna:
            antenna_list=antenna_list[np.arange(len(antenna_list))[(antenna_list!=int(flag))]]
        file['SELECTED_DIPOLES'] = list(antenna_list)

sample_calibrated=1     #Temporary parameter, until block lvl calibration works.

if substation:
    if not sample_calibrated:
        detail_name = '.pol%i.HBA0%s'%(polarization,detail_name)  #HBA0
        detail_name2 = '.pol%i.HBA1%s'%(polarization,detail_name2)   #HBA1
    else:
        detail_name = '.pol%i.sample_calibrated.HBA0%s'%(polarization,detail_name)  #HBA0
        detail_name2 = '.pol%i.sample_calibrated.HBA1%s'%(polarization,detail_name2)   #HBA1
else:
    detail_name = '.pol%i%s'%(polarization,detail_name)
    detail_name2 = ''

file['BLOCKSIZE'] = int(blocklen)

utc = tools.strdate2jd(file['TIME_HR'][0])

#--------------------------------------------------------------------------------------------------------
#Calibration
#--------------------------------------------------------------------------------------------------------
#
if substation and not inner_tiles:
    RCU_selection = [n for n in file['SELECTED_DIPOLES_INDEX'] if n<48]
    RCU_selection2 = [n for n in file['SELECTED_DIPOLES_INDEX'] if n>47]

    if len(RCU_selection)<1 or len(RCU_selection2)<1:
        substation = False
        if len(RCU_selection)<1:
            antenna_set = 'HBA1'
        else:
            antenna_set = 'HBA0'

#---------------------------------------
#RFI excission.
# Find RFI and bad antennas

rfi_channels = [0]  #Always flagging the first channel in case rfi_find=0
freqs = cr.hArray(float,len(file['FREQUENCY_DATA']),fill=file['FREQUENCY_DATA']).setUnit("","Hz")
freqs.setUnit('M','Hz')

if 'LBA' in file['ANTENNA_SET']:
    freq_range = (10, 90)
else:
    freq_range = (110, 190)

if rfi_find:
    rfi_text = '\n**************************** \nSTATION ..........   ' + file['STATION_NAME'][0] + '\nPolarization        ' + str(polarization) + '\n'

    nobl=200  # Standard number of blocks
    sbl=0    # Standard start block

    if 'lotaas' in obs_mode:
        freq_range = (119, 151)

    if 'msec' in obs_mode or file['MAXIMUM_READ_LENGTH']/file['BLOCKSIZE'] < 2*nobl or cr.hArray(file['DATA_LENGTH']).median()[0] * file['SAMPLE_INTERVAL'][0] < 0.04:
        nobl= file['MAXIMUM_READ_LENGTH']/file['BLOCKSIZE']  # Number of blocks
        obs_mode='msec'

    if rfi_find == 2:
        splots = False
    else:
        splots = True

    rfi = cr.trun("FindRFI", f=file,nofblocks=nobl,verbose=False,startblock=sbl,freq_range=freq_range,save_plots=splots,plot_prefix=outdir+'/RFI/FindRFI.'+str(Obs_ID2)+'.st_'+file['STATION_NAME'][0]+'.sbl_'+str(sbl)+'.nobl_'+str(nobl)+'.pol_'+str(polarization)+'.')

#Ploting the average spectrum for each bad antenna.
    color = cr.plt.get_cmap('hot')
    cNorm  = colors.Normalize(vmin=-10, vmax=rfi.nantennas+15)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=color)

    zero_scale = rfi.average_spectrum.toNumpy()[np.where(rfi.average_spectrum.toNumpy()>0.)].min() # to avoid issues with log10(0.)
    average_spectrum = cr.hArray(rfi.average_spectrum+zero_scale)
    median_average_spectrum = cr.hArray(rfi.median_average_spectrum+zero_scale)
    average_spectrum.log10()
    median_average_spectrum.log10()
    cr.plt.figure()
    cr.plt.plot([freq_range[0],freq_range[0]],[min(average_spectrum.vec()),max(average_spectrum.vec())],'k--',label='LOFAR observation')
    cr.plt.plot([freq_range[1],freq_range[1]],[min(average_spectrum.vec()),max(average_spectrum.vec())],'k--')
    if rfi.bad_antennas:
        for i,n in enumerate(file['DIPOLE_NAMES']):
            if n in rfi.bad_antennas:
                cr.plt.plot(freqs,average_spectrum[i/2].vec(),label=str(i)+':'+file['SELECTED_DIPOLES'][i/2],color=scalarMap.to_rgba(i/2))
    #Plotting median spectrum for comparison.
    cr.plt.plot(freqs,median_average_spectrum.vec(),label='Median average spectrum')
    if rfi.bad_antennas:
        cr.plt.legend(loc=2,ncol=4,mode="expand",prop={'size':8,},title='Bad/dead antennas')
    else:
        cr.plt.legend(loc=2,mode="expand",prop={'size':8,},title='No bad/dead antennas.')
    cr.plt.xlabel('Frequency [MHz]')
    cr.plt.ylabel('Log- Power Spectra [ADU]')
    cr.plt.title('average_spectrum - First RFI region.')
    pp = outdir+'/RFI/FindRFI.'+str(Obs_ID2)+'.st_'+file['STATION_NAME'][0]+'.sbl_'+str(sbl)+'.nobl_'+str(nobl)+'.pol_'+str(polarization)+'.antenna_average_spectrum.png'
    cr.plt.savefig(pp)

    if 'msec' not in obs_mode:
        sbl = min(file['DATA_LENGTH'])/file['BLOCKSIZE'] -2*nobl    # Start block
        rfi2 = cr.trun("FindRFI", f=file,nofblocks=nobl,verbose=False,startblock=sbl,freq_range=freq_range,save_plots=splots,plot_prefix=outdir+'/RFI/FindRFI.'+str(Obs_ID2)+'.st_'+file['STATION_NAME'][0]+'.sbl_'+str(sbl)+'.nobl_'+str(nobl)+'.pol_'+str(polarization)+'.')

#Ploting the average spectrum for each bad antenna.
        zero_scale2 = rfi2.average_spectrum.toNumpy()[np.where(rfi2.average_spectrum.toNumpy()>0.)].min() # to avoid issues with log10(0.)
        average_spectrum2 = cr.hArray(rfi2.average_spectrum+zero_scale2)
        median_average_spectrum2 = cr.hArray(rfi2.median_average_spectrum+zero_scale2)
        average_spectrum2.log10()
        median_average_spectrum2.log10()
        cr.plt.plot([freq_range[0],freq_range[0]],[min(average_spectrum2.vec()),max(average_spectrum2.vec())],'k--',label='LOFAR observation frequency range')
        cr.plt.plot([freq_range[1],freq_range[1]],[min(average_spectrum2.vec()),max(average_spectrum2.vec())],'k--')
        cr.plt.figure()
        if rfi2.bad_antennas:
            for i,n in enumerate(file['DIPOLE_NAMES']):
                if n in rfi2.bad_antennas:
                    cr.plt.plot(freqs,average_spectrum2[i/2].vec(),label=str(i)+':'+file['SELECTED_DIPOLES'][i/2],color=scalarMap.to_rgba(i/2))
        #Plotting median spectrum for comparison.
        cr.plt.plot(freqs,median_average_spectrum2.vec(),label='Median average spectrum')
        if rfi2.bad_antennas:
            cr.plt.legend(loc=2,ncol=4,mode="expand",prop={'size':8,},title='Bad/dead antennas')
        else:
            cr.plt.legend(loc=2,mode="expand",prop={'size':8,},title='No bad/dead antennas.')
        cr.plt.xlabel('Frequency [MHz]')
        cr.plt.ylabel('Log- Power Spectra [ADU]')
        cr.plt.title('average_spectrum - Second RFI region')
        pp2 = outdir+'/RFI/FindRFI.'+str(Obs_ID2)+'.st_'+file['STATION_NAME'][0]+'.sbl_'+str(sbl)+'.nobl_'+str(nobl)+'.pol_'+str(polarization)+'.antenna_average_spectrum.png'
        cr.plt.savefig(pp2)
    else:
        bad_antennas_id2=[]
        dead_antennas_id2=[]

#Ploting the single dipole integrated spectral power, and pointing out bad antennas.
    bad_antennas_id = [i for i in range(rfi.nantennas) if rfi.f["SELECTED_DIPOLES"][i] in rfi.bad_antennas and rfi.antennas_cleaned_power[i]>0.]
    dead_antennas_id = [i for i in range(rfi.nantennas) if rfi.f["SELECTED_DIPOLES"][i] in rfi.bad_antennas and rfi.antennas_cleaned_power[i]==0.]
    if 'msec' not in obs_mode:
        bad_antennas_id2 = [i for i in range(rfi2.nantennas) if rfi2.f["SELECTED_DIPOLES"][i] in rfi2.bad_antennas and rfi2.antennas_cleaned_power[i]>0.]
        dead_antennas_id2 = [i for i in range(rfi2.nantennas) if rfi2.f["SELECTED_DIPOLES"][i] in rfi2.bad_antennas and rfi2.antennas_cleaned_power[i]==0.]
    rfi_text+= 'Bad antennas  \n'
    rfi_text+=str(sorted(list(set(bad_antennas_id) & set(bad_antennas_id2))))
    rfi_text+='\nDead antennas  \n'
    rfi_text+=str(sorted(list(set(dead_antennas_id) & set(dead_antennas_id2))))
    cr.plt.figure()
    if 'msec' not in obs_mode:
        rfi.antennas_cleaned_power.sqrt()
        rfi2.antennas_cleaned_power.sqrt()
        cr.plt.plot(rfi.antennas_cleaned_power,label='First RFI region.')
        cr.plt.plot(rfi2.antennas_cleaned_power,label='Second RFI region.')
    else:
        rfi.antennas_cleaned_power.sqrt()
        cr.plt.plot(rfi.antennas_cleaned_power,label='Full time of msec file')
    if len(bad_antennas_id)>0:
        cr.plt.plot(bad_antennas_id,rfi.antennas_cleaned_power[bad_antennas_id],'x', c='r', label = 'Bad antennas', markersize=8)
    if len(bad_antennas_id2)>0:
        if len(bad_antennas_id)>0:
            cr.plt.plot(bad_antennas_id2,rfi2.antennas_cleaned_power[bad_antennas_id2],'x', c='r', markersize=8)
        else:
            cr.plt.plot(bad_antennas_id2,rfi2.antennas_cleaned_power[bad_antennas_id2],'x', c='r', label = 'Bad antennas', markersize=8)
    if len(dead_antennas_id)>0:
        cr.plt.plot(dead_antennas_id,rfi.antennas_cleaned_power[dead_antennas_id],'D', c='k', label = 'Dead antennas', markersize=8)
    if len(dead_antennas_id2)>0:
        if len(dead_antennas_id)>0:
            cr.plt.plot(dead_antennas_id2,rfi2.antennas_cleaned_power[dead_antennas_id2],'D', c='k', markersize=8)
        else:
            cr.plt.plot(dead_antennas_id2,rfi2.antennas_cleaned_power[dead_antennas_id2],'D', c='k', label = 'Dead antennas', markersize=8)
    cr.plt.legend(loc=3,prop={'size':8,})
    cr.plt.xlabel('Antenna Index')
    cr.plt.ylabel('Integrated Spectral Power  [sqrt(ADU)]')
    cr.plt.title('antennas_cleaned_power')
    cr.plt.xticks(range(len(file['SELECTED_DIPOLES'])),[file['SELECTED_DIPOLES_INDEX'][i] for i in range(len(file['SELECTED_DIPOLES']))], rotation='vertical')
    p2 = outdir+'/RFI/FindRFI.'+str(Obs_ID2)+'.st_'+file['STATION_NAME'][0]+'.pol_'+str(polarization)+'.antennas_cleaned_power.png'
    cr.plt.savefig(p2)

#Ploting a waterfall with freq vs antenna for both RFI regions.
    if 'msec' in obs_mode:
        cr.plt.figure()
    else:
        cr.plt.figure(figsize=(8, 12))
        cr.plt.subplot(211)
    cr.plt.imshow(cr.hArray_toNumpy(average_spectrum-median_average_spectrum),aspect='auto',origin='lower',interpolation='nearest',cmap=cr.plt.cm.jet,vmin=-1,vmax=1,extent=(freqs[0],freqs[-1],0,rfi.nantennas))
    cb=cr.plt.colorbar(pad=0.0,fraction=.1)
    cb.set_label('log( Power spectrum per antenna  / median)')
    cr.plt.xlabel('Frequency [MHz]')
    cr.plt.yticks(range(len(file['SELECTED_DIPOLES'])),[file['SELECTED_DIPOLES_INDEX'][i] for i in range(len(file['SELECTED_DIPOLES']))],fontsize=8)
    cr.plt.ylabel('Antenna Index')
    cr.plt.title('First RFI region.')
    if 'msec' not in obs_mode:
        cr.plt.subplot(212)
        cr.plt.imshow(cr.hArray_toNumpy(average_spectrum2-median_average_spectrum2),aspect='auto',origin='lower',interpolation='nearest',cmap=cr.plt.cm.jet,vmin=-1,vmax=1,extent=(freqs[0],freqs[-1],0,rfi.nantennas))
        cb=cr.plt.colorbar(pad=0.0,fraction=.1)
        cb.set_label('log( Power spectrum per antenna / median )')
        cr.plt.xlabel('Frequency [MHz]')
        cr.plt.yticks(range(len(file['SELECTED_DIPOLES'])),[file['SELECTED_DIPOLES_INDEX'][i] for i in range(len(file['SELECTED_DIPOLES']))],fontsize=8)
        cr.plt.ylabel('Antenna Index')
        cr.plt.title('Second RFI region.')
    p3 = outdir+'/RFI/FindRFI.'+str(Obs_ID2)+'.st_'+file['STATION_NAME'][0]+'.pol_'+str(polarization)+'.waterfall_antennas_power.png'
    cr.plt.savefig(p3)

#Continuing with RFI annalysis.
    if 'msec' in obs_mode:
        rfi2 = rfi

    if rfi or rfi2:
        if rfi_find > 1:
            file['SELECTED_DIPOLES'] = list(set(rfi.good_antennas) & set(rfi2.good_antennas))
        rfi_channels = list(set(rfi.dirty_channels) | set(rfi2.dirty_channels))

#Saving simple log of RFI analysis.
    rfi_text+='\nNumber of dirty channels   '+str(len(rfi_channels))
    channels = (file['BLOCKSIZE']/2+1) * (freq_range[1] - freq_range[0])/100.
    dirtyness = (len(rfi_channels) - (len(freqs) - channels)) / channels *100
    rfi_text+='\nDirtyness  %i %% \n'%(dirtyness)
    if int(dirtyness) > 15:
        rfi_text+='**Warning**: dirtyness > 15%  \n'

    if rfi.bad_antennas or rfi2.bad_antennas:
        #This is a sort of place holder for future developement.
        print 'FindRFI found some bad antennas: ', sorted(list(set(rfi.bad_antennas) | set(rfi2.bad_antennas)))

    rfi_info = open(outdir+'/RFI/rfi_info_summary.txt','a')
    rfi_info.write(rfi_text)
    rfi_info.close()

if rfi_find==1:
    raise GeneratorExit('WARNING: rfi_find=1, thus only performing rfi excision, selected rfi_find=2, if also want to do beamforming.')

#To remove the strong RFI line at 170 MHz
#Anna uses a gausian convolution to avoid having "ripples". This happens to them once converting back to timeseries, and they also use all the frequency range form 110 to 190.
#For our case maybe less important. Another way to deal with it (if gaussian convolution is necessary) is to add it to the bandpass.
if 'HBA' in file['ANTENNA_SET'] and rfi_find:
    bad_channels = freqs.Find('between',169.,171.)
    rfi_channels = list(set(rfi_channels) | set(bad_channels))

#---------------------------------------
#Grain calibration
#This calibration comes from .../lofarsoft/data/lofar/StaticMetaData/CalTables , but is too old (2012) to be accurate. Thus, it will not be developed until this is updated.
#Also, this calibration is at the RCU level, and most of the corrections are arond 5%. Thus, it is not be a big issue for us.

#---------------------------------------
#Bandpass Correction
bandpass = 1

if bandpass and rfi_find:
    detail_name +='.bandpass'
    detail_name2 +='.bandpass'

    if rfi:
        median_average_spectrum = cr.hArray(rfi.median_average_spectrum)   # Could also use the mean of rfi and rfi2
        median_average_spectrum.par.xvalues=rfi.f['FREQUENCY_DATA'].setUnit('M','Hz')
    else:
        median_average_spectrum = cr.hArray(rfi2.median_average_spectrum)   #This in case there is no rfi. But, could also use the mean of rfi and rfi2
        median_average_spectrum.par.xvalues=rfi2.f['FREQUENCY_DATA'].setUnit('M','Hz')

    fitbaseline_frat_avg = cr.trerun("FitBaseline","",spectrum=median_average_spectrum,fittype='POLY',numin=freq_range[0],numax=freq_range[1],doplot=3)
    calcbaseline_frat_avg = cr.trerun("CalcBaseline","",fitbaseline_frat_avg.spectrum,doplot=3,normalize=False)

    #Looking for NaNs.
    if np.isfinite(np.sum(calcbaseline_frat_avg.baseline.toNumpy())):
        bandpass = calcbaseline_frat_avg.baseline*calcbaseline_frat_avg.spectrum.median()
    else:
        bandpass = cr.hArray(float, [file['FFTSIZE']], fill=1.0 , name="bandpass")
        print 'WARNING: There are NaNs in the bandpass calibration!!'

    bandpass.sqrt()   # Since multiplying to the complex stokes, and those we see the square when looking at the power spectrum.

no_beamforming = 0
if no_beamforming:
    raise GeneratorExit('STOPING. Since no_beamforming = 1. No bother me anymore.')

#--------------------------------------------------------------------------------------------------------
#Beamforming
#--------------------------------------------------------------------------------------------------------
#Get pointing and phase center info.

FRATS_TRIGGER=os.environ["FRATS_TRIGGER"].rstrip('/')+'/'
BF_PARSETS=os.environ["BF_PARSETS"].rstrip('/')+'/'
fullparsetname=BF_PARSETS+Obs_ID+'.parset'

Trigger_info_file = glob.glob(FRATS_TRIGGER+Obs_ID2+'*_'+default_input+'.npy')

if len(Trigger_info_file) < 1:
    raise IOError('No such file: ' + FRATS_TRIGGER+Obs_ID2+'*_'+default_input+'.npy')

trigger_info = np.load(Trigger_info_file[0])
beamnr = int(trigger_info['beam'])
dm = float(trigger_info['DM'])

try:
    alpha = trigger_info['RA'][0] # Right assention
    delta = trigger_info['DEC'][0] # Declination
except:
    par=bf.get_parameters_new(fullparsetname,useFilename=True)
    alpha = par['beam'][beamnr]['RA'] # Right assention
    delta = par['beam'][beamnr]['DEC'] # Declination

angl_offset = 0.0   # Used for testing of offset pointings (in alpha) in deg.
pointings = getRADEC_2_AZEL(alpha, delta, utc,angl_offset=angl_offset)

cal_delays = dict(zip(file["DIPOLE_NAMES"],file["DIPOLE_CALIBRATION_DELAY"]))

#for testing in coma.
if 'L74100' in fullparsetname and not 'globalhome' in fullparsetname:
    outdir = FRATS_ANALYSIS+Obs_ID+'_local'

if not parset_time:  #used for clock correction offset.
    fullparsetname=''

if angl_offset:
    detail_name += '_offset_%f'%(angl_offset)
    detail_name2 += '_offset_%f'%(angl_offset)

if 'trigger' not in default_input:
    detail_name += '.%s'%(default_input)
    detail_name2 += '.%s'%(default_input)

#----------------------------------------------------
if substation:
    antenna_set = 'HBA0'
    file['SELECTED_DIPOLES'] = RCU_selection
    sample_offset = get_sample_offset(fullparsetname,file, antenna_set,max_sample_number,sample_calibrated=sample_calibrated, ref_station = ref_station, ref_time = ref_time)
    bm = cr.trun("BeamFormer", filenames = [fname,], output_dir = outdir+'/beam.results/', pointings = pointings, cal_delays = cal_delays, blocklen = blocklen,detail_name=detail_name,sample_offset=sample_offset,station_centered=station_centered,antenna_set=antenna_set,antenna_list=file['SELECTED_DIPOLES_INDEX'],rfi_channels=rfi_channels, bandpass = bandpass)

    antenna_set = 'HBA1'
    file['SELECTED_DIPOLES'] = RCU_selection2
    sample_offset = get_sample_offset(fullparsetname,file, antenna_set,max_sample_number,sample_calibrated=sample_calibrated, ref_station = ref_station, ref_time = ref_time)
    bm2 = cr.trun("BeamFormer", filenames = [fname,], output_dir = outdir+'/beam.results/', pointings = pointings, cal_delays = cal_delays, blocklen = blocklen,detail_name=detail_name2,sample_offset=sample_offset,station_centered=station_centered,antenna_set=antenna_set,antenna_list=file['SELECTED_DIPOLES_INDEX'],rfi_channels=rfi_channels, bandpass = bandpass)
else:
    sample_offset = get_sample_offset(fullparsetname,file, file['ANTENNA_SET'],max_sample_number,sample_calibrated=sample_calibrated, ref_station = ref_station, ref_time = ref_time)
    bm = cr.trun("BeamFormer", filenames = [fname,], output_dir = outdir+'/beam.results/', pointings = pointings, cal_delays = cal_delays, blocklen = blocklen,detail_name=detail_name,sample_offset=sample_offset,station_centered=station_centered,antenna_list=file['SELECTED_DIPOLES_INDEX'],rfi_channels=rfi_channels, bandpass = bandpass)

#----------------------------------------------------

print 'If you can read this, then all went good and you have some new beams waiting for you.'
