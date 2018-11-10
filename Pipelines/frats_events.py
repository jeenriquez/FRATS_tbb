#! /usr/bin/env python
"""
**Usage:**
execfile(PYP+"pipelines/frats_events.py")

This is the top layer to run the frats_event pipeline in multiple files from the same observation.

**Parameters:**
======================= =================================================================
*obs_name*              Observation name of raw data to be processed.

*substation*            If True, then using HBA substations (HBA0 and HBA1)

======================= =================================================================

For more help on more parameters run ``frats_events.py --help``.

Examples::

python frats_events.py L184417_D20131103T1421 -r 0 -c
python frats_events.py L184417_D20131103T1421 -v 1 -s core -x '-m lotaas'

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
V1.0 created by E. Enriquez, Jan 2014
V1.1 modified by E. Enriquez, Aug 2014
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
#------------------------------------------------------

def make_list(option, opt_str, value, parser):
    setattr(parser.values, option.dest, value.replace('[','').replace(']','').split(','))

def getRef_time(filelist,fullparsetname=''):
    ''' Gets the reference time for a list of tbb files (station which start observing last).
    '''

    times = cr.hArray(float,len(filelist),0.0)
    seconds = cr.hArray(float,len(filelist),0.0)

    for i,file in enumerate(filelist):
        tbb = cr.open(file)

        if fullparsetname:
            f_clock_offset = float(md.getClockCorrectionParset(fullparsetname, tbb['STATION_NAME'][0], antennaset=tbb['ANTENNA_SET']))
        else:
            f_clock_offset = float(md.getClockCorrection(tbb['STATION_NAME'][0], antennaset=tbb['ANTENNA_SET']))

        time_tbb = cr.hArray(float,len(tbb['SAMPLE_NUMBER']),fill=tbb['TIME'])
        time_tbb -= min(tbb['TIME'])
        sample_tbb = cr.hArray(float,len(tbb['SAMPLE_NUMBER']),fill=tbb['SAMPLE_NUMBER'])
        sample_tbb *= tbb['SAMPLE_INTERVAL'][0]
        sample_tbb += time_tbb

        times[i] = max(sample_tbb)+f_clock_offset
        seconds[i] = min(tbb['TIME'])

        tbb.close()

    if verbose:
        print 'Reference time: ' , max(times) , 'from file ', filelist[cr.hFindBetweenOrEqual(times,max(times),max(times))[0]]
    return [min(seconds),max(times)]

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
parser.add_option("-x","--extra_parameters",type="str",default='',help="Extra parameters from the pipeline.")

(options, args) = parser.parse_args()

#----------------------------------------------------
if not parser.get_prog_name()=="frats_events.py":
    #   Program was run from within python
    station_centered = False
    stations = ['superterp']
    dump_version = 'R000'
    rfi_find = 2
    extra_parameters=''
    parset_time = True
    verbose = 2
else:
    station_centered = options.station_centered
    stations = options.stations
    dump_version = options.dump_version
    rfi_find = options.rfi_find
    extra_parameters = options.extra_parameters
    parset_time = options.parset_time
    verbose = options.verbose

#----------------------------------------------------
#Chooosing TBB files.

obs_name = args[0]
FRATS_DATA=os.environ["FRATS_DATA"].rstrip('/')+'/'
BF_PARSETS=os.environ["BF_PARSETS"].rstrip('/')+'/'
fullparsetname=BF_PARSETS+obs_name.rsplit('_')[0]+'.parset'

if stations[0] == 'superterp':
    beam_suffix = '*CS00[2-7]_'+dump_version+'*'
elif stations[0] == 'core':
    beam_suffix = '*CS*'+dump_version+'*'
elif stations[0] == 'all_nl':
    beam_suffix = '*[C,R]S???*'+dump_version+'*'
else:
    beam_suffix = '*'+stations[0]+'*'

filenames=glob.glob(FRATS_DATA+obs_name+'*/'+beam_suffix)

if stations[0] == 'superterp' and len(filenames) > 6:
        raise ValueError('There seem to be several dump versions (R00?), please remove unwanted files.')
elif stations[0] == 'core' and len(filenames) > 24:
        raise ValueError('There seem to be several dump versions (R00?), please remove unwanted files.')
elif stations[0] == 'all_nl' and len(filenames) > 38:
        raise ValueError('There seem to be several dump versions (R00?), please remove unwanted files.')

if verbose:
    print "Running pipeline on: \n"
    for name in filenames:
        print name
    print '------_-------'

#----------------------------------------------------
#General parameters.

if parset_time:
    ref_time = getRef_time(filenames,fullparsetname=fullparsetname)
else:
    ref_time = getRef_time(filenames,fullparsetname='')
    extra_parameters += ' -p'

if station_centered:
    station_centered='-s'
else:
    station_centered=''

#----------------------------------------------------
#Running Pipeline
#----------------------------------------------------
#Looping over files, polarization, rfi types? and phase center.
for file in filenames:
    for pol in range(2):
        if verbose:
            print "python"+" frats_event.py "+file+" -p "+str(pol)+" -r "+str(rfi_find)+" --ref_time "+str(ref_time).replace(' ','')+' '+station_centered+' '+extra_parameters
        if verbose !=1:
            proc=subprocess.Popen(["ionice","-c","3","python","frats_event.py",file,"-p",str(pol),"-r",str(rfi_find),"--ref_time",str(ref_time).replace(' ',''),station_centered,extra_parameters])

#run $PROS/frats_event.py $FRATS_DATA/L74100_ev1/L74100_D20121107T191144.964Z_CS002_R000_tbb.h5
#run $PROS/frats_event.py $FRATS_DATA/L83047_ev1/L83047_D20130116T003817.425Z_CS002_R000_tbb.h5
#run $PROS/frats_event.py $FRATS_DATA/L83047_ev1/L83047_D20130116T003817.425Z_CS002_R000_tbb.h5 -f -1 -p 1 -s
#run $PROS/frats_event.py $FRATS_DATA/L83047_ev1/L83047_D20130116T003817.425Z_CS002_R000_tbb.h5 -s
