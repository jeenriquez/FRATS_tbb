#! /usr/bin/env python
"""
**Usage:**
execfile(PYP+"pipelines/localize_frats.py")

This pipeline does runs a series of scripts and tasks to make a high angular resolution image of the pulse.

**Parameters:**
======================= =================================================================
*obs_name*              Observation name of raw data to be processed.

*substation*            If True, then using HBA substations (HBA0 and HBA1)

======================= =================================================================

For more help on more parameters run ``localize_frats.py --help``.

Examples::

python localize_frats.py L184417_D20131103T1421 -r 0 -c
python localize_frats.py L184417_D20131103T1421 -v 1 -s core -x '-m lotaas'

Note:
It needs environment variables: $BF_PARSETS, $FRATS_ANALYSIS, $FRATS_TRIGGER $FRATS_DATA

    In Coma:
        $BF_PARSETS   = /vol/astro1/lofar/frats/bf/parsets
        $FRATS_ANALYSIS = /vol/astro1/lofar/frats/tbb/analysis
        $FRATS_TRIGGER = /vol/astro1/lofar/frats/tbb/analysis/trigger_info/

    In locus013:
        $BF_PARSETS   = /globalhome/lofarsystem/log/
        $FRATS_ANALYSIS = /data/TBB/FRATS/tbb/analysis/
        $FRATS_TRIGGER = /home/enriquez/FRATS/trigger_info/
        $FRATS_DATA = /data/TBB/FRATS/tbb/data/

.. moduleauthor:: J. Emilio Enriquez <e.enriquez 'at' astro.ru.nl>

Revision History:
V1.0 created by E. Enriquez, Oct 2014
"""

import matplotlib
matplotlib.use("Agg")

from optparse import OptionParser
import pycrtools as cr
from pytmf import *
import subprocess
import numpy as np
from pycrtools import tools
from pycrtools import bfdata as bf
from pycrtools import metadata as md
import os; import sys; import glob
import pdb;# pdb.set_trace()
import Beam_Tools as bt

#------------------------------------------------------

def make_list(option, opt_str, value, parser):
    setattr(parser.values, option.dest, value.replace('[','').replace(']','').split(','))

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

parser.add_option("-c","--station_centered",action="store_true",default=False,help="True if need phase center at (sub)station center for imaging, or False at LOFAR center for dynamic spectra.")
parser.add_option("-s","--stations",type="str",action='callback',default=['superterp'],callback=make_list,help="Stations for which to make a beam.")
parser.add_option("-d","--dump_version",type="str",default='R00?',help="Use different dump versions (eg. R000,R001,...)")
parser.add_option("-r","--rfi_find",type="int",default=2,help="If 1, then performing findrfi task only; 0 for beamform only; 2 for both.")
parser.add_option("-v","--verbose",type="int",default=2,help="If 1, then only printing without doing; 0 for only doing without printing; 2 for both.")
parser.add_option("--parset_time",action="store_false",default=True,help="Using this parameter to either use the clock_offsets from the parset files or the static values.")
parser.add_option("-p","--polarization",type="int",default=0,help="Polarization, 0 for even, 1 for odd antennas.")
parser.add_option("-u","--substation",action="store_false",default=True,help="If True, then using HBA substations (HBA0 and HBA1).")
parser.add_option("-b","--blocklen",type="int",default=2**14,help="Block lenght")
parser.add_option("-t","--tbin",type="int",default=1,help="Time binning.")
parser.add_option("--default_input",type="str",default='trigger',help='Extra name for input file with different DM, pointing, etc..')
parser.add_option("-a","--analysis",type="str",default='TAB',help="Type of analysis, could be TAB- for Tied Array, IMG - for imaging, NONE - to just have the beams loaded. ")
parser.add_option("--incoherence",action="store_true",default=False,help="False if using the phase information for the final beam.")

(options, args) = parser.parse_args()

#----------------------------------------------------
if not parser.get_prog_name()=="frats_event_localize.py":
    #   Program was run from within python
    station_centered = False
    stations = ['superterp']
    dump_version = 'R000'
    rfi_find = 2
    parset_time = True
    verbose = 2
    polarization = 0
    rfi_find = 0
    substation = 1
    blocklen = 2**14
    tbin = 16
    freq_range = [130,165]
    time_range = [0, 0.05]
    default_input = 'trigger'
    analysis = 'TAB'
    incoherence = False
else:
    station_centered = options.station_centered
    stations = options.stations
    dump_version = options.dump_version
    rfi_find = options.rfi_find
    parset_time = options.parset_time
    verbose = options.verbose
    substation = options.substation
    polarization = options.polarization
    blocklen = options.blocklen
    rfi_find = options.rfi_find
    tbin = options.tbin
    default_input = options.default_input
    analysis = options.analysis
    incoherence = options.incoherence

#----------------------------------------------------
#Chooosing BEAM files.
obs_name = args[0]
FRATS_ANALYSIS=os.environ["FRATS_ANALYSIS"].rstrip('/')+'/'
FRATS_TRIGGER=os.environ["FRATS_TRIGGER"].rstrip('/')+'/'
BF_PARSETS=os.environ["BF_PARSETS"].rstrip('/')+'/'
fullparsetname=BF_PARSETS+obs_name.rsplit('_')[0]+'.parset'

# Create output directories, if not already present.
if not os.path.isdir(FRATS_ANALYSIS+obs_name):
    raise ValueError('Observation '+FRATS_ANALYSIS+obs_name+' not found.')
if not os.path.isdir(FRATS_ANALYSIS+obs_name+'/all_sky/'):
    os.makedirs(FRATS_ANALYSIS+obs_name+'/all_sky/')
if not os.path.isdir(FRATS_ANALYSIS+obs_name+'/TAB_products/'):
    os.mkdir(FRATS_ANALYSIS+obs_name+'/TAB_products/')

detail_name='.pol%i.'%(polarization)

if station_centered:
    detail_name+=''
else:
    detail_name+='*.LOFAR_centered'

if stations[0] == 'superterp':
    beam_suffix = '*CS00[2-7]_'+dump_version+'*'+detail_name
elif stations[0] == 'core':
    beam_suffix = '*CS*'+dump_version+'*'+detail_name
elif stations[0] == 'all_nl':
    beam_suffix = '*[C,R]S???*'+dump_version+'*'+detail_name
else:
    beam_suffix = '*'+stations[0]+'*'+detail_name

filenames=glob.glob(FRATS_ANALYSIS+obs_name+'*/beam.results/'+beam_suffix+'*.beam')

if 'HBA0' in filenames[0] or 'HBA1' in filenames[0]:
    factor = 2
else:
    factor = 1

if stations[0] == 'superterp' and len(filenames) > 6*factor:
        raise ValueError('There seem to be several dump versions (R00?), please remove unwanted files.')
elif stations[0] == 'core' and len(filenames) > 24*factor:
        raise ValueError('There seem to be several dump versions (R00?), please remove unwanted files.')
elif stations[0] == 'all_nl' and len(filenames) > 38*factor:
        raise ValueError('There seem to be several dump versions (R00?), please remove unwanted files.')

if verbose:
    print "Running pipeline on: \n"
    for name in filenames:
        print name
    print '------_-------'

#----------------------------------------------------
#General parameters.

try:
    Trigger_info_file = glob.glob(FRATS_TRIGGER+obs_name+'*_'+default_input+'.npy')
except:
    raise IOError('No such file: ' + FRATS_TRIGGER+obs_name+'*_'+default_input+'.npy')

if len(Trigger_info_file) < 1:
    raise IOError('No such file: ' + FRATS_TRIGGER+obs_name+'*_'+default_input+'.npy')

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

angl_offset = 0.0  #Placeholder for future developement. If changing the pointing within the incoherent beam.

#------
#Loading beams
beams = cr.open(sorted(filenames))
#------

utc = tools.strdate2jd(beams['TBB_TIME_HR'][0])
pointings = getRADEC_2_AZEL(alpha, delta, utc,angl_offset=angl_offset)

#----------------------------------------------------
#Running Pipeline
#Pulse localization
#IDEAS:
#   - could integrate over the same frequency bands that sander uses and find the peak positions.

#----------------------------------------------------
#DM selection
beams['DM'] = dm
print 'Using DM = ', dm

#------
#Inter Station time-delay self-calibration by Cross Correlation.
if not incoherence:
    ST_DELAYS = bt.ccBeams(beams,freq_range=freq_range,time_range=time_range,verbose=verbose)
    beams['CAL_DELAY'] = ST_DELAYS

#------
#BEAM addiction (coherent if having ST_DELAYS, otherwise using only the amplitudes.)

if 'TAB' in analysis and not station_centered:

    if incoherence:
        TAB,dynspec,cleandynspec = bt.addBeams(beams,dyncalc=True,tbin=tbin,dm=dm,clean=True,incoherent=True)
    else:
        TAB,dynspec,cleandynspec = bt.addBeams(beams,dyncalc=True,tbin=tbin,dm=dm,clean=True)

#print "The dynspecs are saved here:"

#------
#Beam identification
# Pulse SNR calculation.
'''
    zoom_dynspec = bt.cutDynspec(cleandynspec,freq_range=freq_range,time_range=time_range)
    flat_cleandyn = bt.flatDynspec(zoom_dynspec,axis='x')
    stats = bt.doStats(flat_cleandyn,verbose=True,nof_St=beams['NOF_BEAM_DATASETS'])


    extent=(cleandynspec.par.xvalues[0],cleandynspec.par.xvalues[-1],cleandynspec.par.yvalues[0],cleandynspec.par.yvalues[-1])
    cr.plt.imshow(np.log10(cr.hArray_toNumpy(cleandynspec)),aspect='auto',origin='lower',cmap=cr.plt.cm.hot,extent=extent)
'''
    #save fits, maybe also png?

#------
#Imaging with the beamager
if 'IMG' in analysis and station_centered:
        raise NotImplementedError

