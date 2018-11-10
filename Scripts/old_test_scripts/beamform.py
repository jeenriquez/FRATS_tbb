#! /usr/bin/env python

"""Run beamformer on Pulsar: /vol/astro/lofar/frats/tbb/data/.
"""

import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pycrtools as cr
from pytmf import *
from pycrtools import tools
from pycrtools.tasks import beamformer

#----------------------------------------------------
fname = sys.argv[1]
f = cr.open(fname)

utc = tools.strdate2jd(f['TIME_HR'][0])

#----------------------------------------------------
outdir = '../beam.results/beams.sub_station/'
substation = 1
inner = 0
detail_name = '.pol0a.new.multi_station'
detail_name2 = '.pol0b.new.multi_station'
#----------------------------------------------------
#Pointing and LOFAR posisition info.
alpha = hms2rad(3,32,59.37) # Right assention
delta = dms2rad(54,34,44.9) # Declination

phi = deg2rad(52.915122495) #(LOFAR Superterp)
L = deg2rad(6.869837540)
#angl_offset = deg2rad(0.0)
#delta = delta+angl_offset
#alpha = alpha+angl_offset

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
cal_delays = dict(zip(f["DIPOLE_NAMES"],f["DIPOLE_CALIBRATION_DELAY"]))
#----------------------------------------------------
#Calibrate station times.
f0=cr.open('/vol/astro/lofar/frats/tbb/data/L43784_D20120125T211154.887Z_CS007_R000_tbb.h5')

t0 = max(f0['SAMPLE_NUMBER'])*f0['SAMPLE_INTERVAL'][0]+f0['CLOCK_OFFSET'][0]
t = max(f['SAMPLE_NUMBER'])*f['SAMPLE_INTERVAL'][0]+f['CLOCK_OFFSET'][0]

sample_offset = int((t0-t)/f['SAMPLE_INTERVAL'][0])

if inner:
    antenna_list = cr.hArray(int,(4),[17,19,29,31])

#----------------------------------------------------

if substation:
    antenna_set = 'HBA0'
    if inner:
        bm = cr.trun("BeamFormer", filenames = [fname,], output_dir = outdir, pointings = pointings, FarField = True, NyquistZone = 2, cal_delays = cal_delays, qualitycheck = False, maxchunklen = 2**20, blocklen = 2**14, dohanning = True,detail_name=detail_name,sample_offset=sample_offset,single_station=True,antenna_set=antenna_set,antenna_list=list(antenna_list))
    else:
        bm = 0#cr.trun("BeamFormer", filenames = [fname,], output_dir = outdir, pointings = pointings, FarField = True, NyquistZone = 2, cal_delays = cal_delays, qualitycheck = False, nantennas_start = 0,  nantennas_stride = 2, maxnantennas = f["NOF_DIPOLE_DATASETS"]/2, maxchunklen = 2**20, blocklen = 2**14, dohanning = True,detail_name=detail_name,sample_offset=sample_offset,single_station=False,antenna_set=antenna_set)
    antenna_set = 'HBA1'
    if inner:
        bm2 = cr.trun("BeamFormer", filenames = [fname,], output_dir = outdir, pointings = pointings, FarField = True, NyquistZone = 2, cal_delays = cal_delays, qualitycheck = False, maxchunklen = 2**20, blocklen = 2**14, dohanning = True,detail_name=detail_name2,sample_offset=sample_offset,single_station=True,antenna_set=antenna_set,antenna_list=list(antenna_list+48))
    else:
        bm2 = cr.trun("BeamFormer", filenames = [fname,], output_dir = outdir, pointings = pointings, FarField = True, NyquistZone = 2, cal_delays = cal_delays, qualitycheck = False, nantennas_start = 48,  nantennas_stride = 2, maxnantennas = f["NOF_DIPOLE_DATASETS"], maxchunklen = 2**20, blocklen = 2**14, dohanning = True,detail_name=detail_name2,sample_offset=sample_offset,single_station=False,antenna_set=antenna_set)
else:
    bm = cr.trun("BeamFormer", filenames = [fname,], output_dir = outdir, pointings = pointings, FarField = True, NyquistZone = 2, cal_delays = cal_delays, qualitycheck = False, nantennas_start = 1,  nantennas_stride = 2, maxnantennas = f["NOF_DIPOLE_DATASETS"], maxchunklen = 2**20, blocklen = 2**14, dohanning = True,detail_name=detail_name,sample_offset=sample_offset,single_station=False)
