''' Script to test the cross correlation of pulses. Old version. 
    (Incorrect since adding from the dynspec instead of the complex beams...but need the phases information too....)
'''

import pycrtools as cr
import numpy as np
import sys; import os
import pdb;# pdb.set_trace()
import Beam_Tools as bt

#filenames = ['L43784_D20120125T211154.887Z_CS002_R000_tbb.pol0.multi_station.beam','L43784_D20120125T211154.867Z_CS003_R000_tbb.pol0.multi_station.beam','L43784_D20120125T211154.887Z_CS004_R000_tbb.pol0.multi_station.beam','L43784_D20120125T211154.866Z_CS005_R000_tbb.pol0.multi_station.beam','L43784_D20120125T211154.871Z_CS006_R000_tbb.pol0.multi_station.beam','L43784_D20120125T211154.887Z_CS007_R000_tbb.pol0.multi_station.beam']
filenames = ['L43784_D20120125T211154.887Z_CS002_R000_tbb.pol1.multi_station.beam','L43784_D20120125T211154.867Z_CS003_R000_tbb.pol1.multi_station.beam','L43784_D20120125T211154.887Z_CS004_R000_tbb.pol1.multi_station.beam','L43784_D20120125T211154.866Z_CS005_R000_tbb.pol1.multi_station.beam','L43784_D20120125T211154.871Z_CS006_R000_tbb.pol1.multi_station.beam','L43784_D20120125T211154.887Z_CS007_R000_tbb.pol1.multi_station.beam']
#filenames = ['L43784_D20120125T211154.887Z_CS002_R000_tbb.pol1.multi_station.beam','L43784_D20120125T211154.867Z_CS003_R000_tbb.pol1.multi_station.beam','L43784_D20120125T211154.887Z_CS004_R000_tbb.pol1.multi_station.beam']

#extra = 'dd_dynspec.pcr/data.bin' 
extra = 'clean_dynspec.pcr/data.bin' 
cr.plt.ion()

verbose=0

nblocks=640
freq_size=1638
full_blocks=12224

N=len(filenames)
label=['']

flat_matrix=cr.hArray(float,[N,nblocks])

for i,file in enumerate(filenames):
    beams = cr.open(file)
    print 'Working on station: '+ beams['STATION_NAME'][0]
    full_dynspec=cr.hArray(float,(beams['BEAM_SPECLEN'],full_blocks))
    label.append(beams['STATION_NAME'][0])
    
    filename = os.path.join(file,extra)
    full_dynspec.readfilebinary(filename,0)
    pdb.set_trace()
    zoom_dynspec=bt.cutdynspec(full_dynspec,pulse=0)
    zoom_dynspec.par.station = beams['STATION_NAME'][0]

    flat_dyn=bt.flatdyn(zoom_dynspec)

    #Create a time vector
    times = zoom_dynspec.par.xvalues

    if verbose:
        #Create a frequency vector
        frequencies = zoom_dynspec.par.yvalues    

        cr.plt.clf()
        cr.plt.title("Pulse Profile from : " + zoom_dynspec.par.station)
        cr.plt.plot(flat_dyn.par.xvalues,flat_dyn)
        cr.plt.xlabel("+/- Time [s]")
        cr.plt.autoscale(axis='x',tight=True)
        pdb.set_trace()

    flat_dyn=cr.hArray(float,len(flat_dyn),flat_dyn)
    flat_matrix[i] = flat_dyn
    

#flat_matrix[...].runningaverage(15,cr.hWEIGHTS.GAUSSIAN)
#flat_matrix[...].plot()

pdb.set_trace()
#------------------------------------------------------
#Cross Correlation using fftw
    
fft_matrix = cr.hArray(complex,[N,nblocks/2+1])
fft_matrix[...].fftw(flat_matrix[...])
    
#Cross correlation: vec1*conj(vec2)
fft_matrix[1:,...].crosscorrelatecomplex(fft_matrix[0], True)
fft_matrix[0:1,...].crosscorrelatecomplex(fft_matrix[0], True)

croscorr=cr.hArray(float,[N,nblocks])
croscorr[...].invfftw(fft_matrix[...])
croscorr /= nblocks
croscorr.abs()
#croscorr.runningaverage(15,cr.hWEIGHTS.GAUSSIAN)

offset=cr.hArray(float,nblocks,range(nblocks))-nblocks/2.

verbose=0
if verbose:
    cr.plt.rcParams['font.size']=10
    cr.plt.ion()
    
    for i in range(N):
        cr.plt.subplot(N+1,1,i+1)
        cr.plt.title(label[i+1])
        cr.plt.plot(times,flat_matrix[i].vec())
        cr.plt.autoscale(axis='x',tight=True)    
    cr.plt.xlabel("+/- Time [s]")
    
    cr.plt.subplot(N+1,1,N+1)

for i in range(N):
    cr.plt.plot(offset,croscorr[i].vec(),label='Station 2 vs '+label[i+1])

if not verbose:
    cr.plt.legend()
cr.plt.plot([0,0],[croscorr.vec().min(),croscorr.vec().max()])
cr.plt.autoscale(axis='x',tight=True)

