--------------------------------------------------------
#Example script (CR)

import pycrtools as cr
from pycrtools.tasks import findrfi
import matplotlib.pyplot as plt

datafile = cr.open('VHECR_example.h5')


plt.ion()
plt.figure()

--------------------------------------------------------
#FindRFI

rfi = cr.trun("FindRFI", f = datafile, plot_prefix ='',nofblocks=10,verbose=False,startblock=0,save_plots=True)

if rfi:
    datafile['SELECTED_DIPOLES'] = rfi.good_antennas

--------------------------------------------------------
#Imaging
from pycrtools.tasks import imager

output = 'out_cr.test.fits'

cr.trun("Imager", data = datafile, intgrfreq = True, nblocks=1,ntimesteps=3,startblock=8,NAXIS1=155,NAXIS2=155,output=output,CDELT1=-0.5,CDELT2=0.5, rfi_remove=rfi.dirty_channels,CRVAL1=231.9,CRVAL2=63.1,CTYPE1='ALON-AZP',CTYPE2='ALAT-AZP', FREQMIN=0.3e8,FREQMAX=0.8e8)
----------
