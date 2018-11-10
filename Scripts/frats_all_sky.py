'''This script creates an all sky image.
'''

from optparse import OptionParser
import pycrtools as cr
from pycrtools.tasks import beamformer, imager, findrfi
import os; import sys; import glob
from pytmf import *


beamforming = 0
imaging = 1

#Adding if statement, in case no need to do the "beams" and only the imaging part, ugly fix for now.

#----------------------------------------------------
#Hard coded stuff that will go
filename='/vol/astro/lofar/frats/tbb/data/L43784_ev1/L43784_D20120125T211154.887Z_CS002_R000_tbb.h5'

polarization=0
blocklen = 2**14
dm =26.76
freq_range = [166,168]

#----------------------------------------------------
#File selection
#Open given TBB file.
fname = filename# args[0]
file = cr.open(fname)

#----------------------------------------------------
#General parameters.

FRATS_ANALYSIS=os.environ["FRATS_ANALYSIS"].rstrip('/')+'/'
Obs_ID = fname.split('/')[-1].split('_')[0]

outdir = FRATS_ANALYSIS+Obs_ID+'/all_sky/'

# Create output directory, if not already present.
if not os.path.isdir(outdir) or not os.path.isdir(outdir+'/beams/'):
    os.makedirs(outdir+'/beams/')
if not os.path.isdir(outdir+'/RFI/'):
    os.mkdir(outdir+'/RFI/')

file['BLOCKSIZE'] = blocklen

#Antenna Selection
poli=['even','odd']
file['SELECTED_DIPOLES']=poli[polarization]
antenna_list = file['SELECTED_DIPOLES_INDEX']

if beamforming:
    #----------------------------------------------------
    #RFI excission.
    #----------------------------------------------------
    # Find RFI and bad antennas
    nobl=200  # Number of blocks
    sbl=0    # Start block
    rfi = cr.trun("FindRFI", f=file,nofblocks=nobl,verbose=False,startblock=sbl,freq_range=(110, 190),save_plots=True,plot_prefix=outdir+'/RFI/FindRFI_sbl_'+str(sbl)+'_nobl_'+str(nobl)+'_st_'+file['STATION_NAME'][0]+'_pol_'+str(polarization)+'_')

    sbl= min(file['DATA_LENGTH'])/file['BLOCKSIZE'] -2*nobl    # Start block
    rfi2 = cr.trun("FindRFI", f=file,nofblocks=nobl,verbose=False,startblock=sbl,freq_range=(110, 190),save_plots=True,plot_prefix=outdir+'/RFI/FindRFI_sbl_'+str(sbl)+'_nobl_'+str(nobl)+'_st_'+file['STATION_NAME'][0]+'_pol_'+str(polarization)+'_')

    if rfi:
        file['SELECTED_DIPOLES'] = list(set(rfi.good_antennas) & set(rfi.good_antennas))
        rfi_channels = list(set(rfi.dirty_channels) | set(rfi.dirty_channels))
    else:
        rfi_channels = None

    #----------------------------------------------------
    #Beamform
    #----------------------------------------------------

    pointings = [{'az': deg2rad(0), 'el': deg2rad(90)}]

    print 'Running beamformer in : ',fname

    #Creating single antenna "beams".
    for ant in antenna_list:
        detail_name='_test_single_antenna_%s'%(ant)
        file['SELECTED_DIPOLES']=[ant]
        cal_delays = dict(zip(file["DIPOLE_NAMES"],file["DIPOLE_CALIBRATION_DELAY"]))
        phase_center = list(file.getRelativeAntennaPositionsNew().vec())
        cr.trun('BeamFormer', filenames = [fname], output_dir = outdir+'beams/', blocklen = blocklen, pointings = pointings, cal_delays = cal_delays, detail_name = detail_name,antenna_list=[ant],phase_center=phase_center,rfi_channels=rfi_channels)


if imaging:
    #----------------------------------------------------
    #All-Sky imaging.
    #----------------------------------------------------

    beam_suffix = '*_test_single_antenna*'
    filenames = glob.glob(outdir+'beams/'+beam_suffix)

    beams = cr.open(filenames)

    beams['DM'] = dm

    startblock = 275
    nblocks = 25

    output = outdir+Obs_ID+'_all_sky__startblock_'+str(startblock)+'_dm_'+str(dm)+'_pol_'+str(polarization)+'_'+str(freq_range[0])+'-'+str(freq_range[1])

    cr.trun("Imager", data = beams, intgrfreq = True, nblocks=nblocks,ntimesteps=1,startblock=startblock,NAXIS1=161,NAXIS2=161,CDELT1=-1.,CDELT2=1.,FREQMIN=freq_range[0]*1e6,FREQMAX=freq_range[1]*1e6,output=output,CRVAL3=freq_range[0]*1e6)

#total runtime: 1693.02941322 s
#[Imager] completed in 40345.580 s
