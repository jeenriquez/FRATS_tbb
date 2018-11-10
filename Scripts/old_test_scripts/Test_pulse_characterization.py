'''Preliminary Pipeline for characterization of pulses.
 CAL_DELAYs are calculated for L43784_D20120125T211154 only.
'''

import pycrtools as cr
import numpy as np
import sys; import os
import glob
import pdb;# pdb.set_trace()
import Beam_Tools as bt


def getPulseInfo(pulse=0,tbin=1):
    '''Used to get the freq and time ranges to cut the pulse.
    '''
    dt_pulse = .005

    if pulse == 1:
        freq_range = [132,146]   # 130#145
        time_range = [0.57,0.6201]
        pulse_range = [.5915,.5915+dt_pulse]
        #pulse_range = [220/tbin,270/tbin]
        #frequency_slice = [2458,3687]# [2622,3769]
#        freq_size=1147+1 #1228
#        frequency_slice = [2622,3770] #[2458,3686]
#            nblocks=560
#            startblock=7000
    elif pulse == 3:       #pulse 3
        freq_range=[149,169]#[155,160]#
        time_range = [0.0,0.0501]
        pulse_range = [.0205,.0205+dt_pulse]
        #pulse_range = [250/tbin,300/tbin]
        #frequency_slice = [4015,5653]    [4506,4916]
#        nblocks=640
#        startblock=0
    elif pulse == 2:       #pulse 2
        freq_range = [140,155]
        time_range = [0.28,0.3301]
        pulse_range = [.308,.308+dt_pulse]
    elif pulse == 4:       #pulse 4
        freq_range = [126,136]
        time_range = [.86,.9101]
        pulse_range = [.881,.881+dt_pulse]
    elif pulse == 5:       #pulse 5
        freq_range = [106,117]
        time_range = [0.715,0.7651]
        pulse_range = [.739,.739+dt_pulse]
    elif pulse == 6:       #pulse 6
        freq_range = [178,190]
        time_range = [0.415,0.4651]
        pulse_range = [.435,.435+dt_pulse]
    elif pulse == 7:       #pulse 7
        freq_range = [126,136]

    info_pulse= cr.hArray()
    info_pulse.par.freq_range = freq_range
    info_pulse.par.time_range = time_range
    info_pulse.par.pulse_range = pulse_range

    return info_pulse

def super_plot(zoom_dynspec,pulse=0,zoom=1,info=None):
    '''Do some awesome ploting.
    '''
    #Create a frequency vector
    frequencies = zoom_dynspec.par.yvalues
    #Create a time vector
    times = zoom_dynspec.par.xvalues

    pdb.set_trace()
    cr.plt.clf()
    cr.plt.ion()
    zoom=1
    if zoom:
        cr.plt.rcParams['font.size']=20
        cr.plt.rcParams['axes.titlesize']=30
        cr.plt.subplot(2,2,1)
        cr.plt.title("Pulse Profile from : " + zoom_dynspec.par.station)
        cr.plt.xlabel("+/- Time [s]")
        flat_dyn=bt.flatdyn(zoom_dynspec,axis='x')
        cr.plt.plot(flat_dyn.par.xvalues,flat_dyn)
        cr.plt.autoscale(axis='x',tight=True)

        cr.plt.subplot(2,2,4)
        cr.plt.xlabel("Frequency [MHz]")
        flat_dyn,flat_dyn2=bt.flatdyn(zoom_dynspec,axis='y',pulse_range=info.par.pulse_range,background_range=[0,50/info.par.tbin])
        cr.plt.plot(flat_dyn.par.xvalues,flat_dyn)
        cr.plt.plot(flat_dyn.par.xvalues,flat_dyn2,'r')
        cr.plt.autoscale(axis='x',tight=True)

        cr.plt.subplot(2,2,3)
    #Time integration.
    tbin=1
    if tbin>1:
        zoom_dynspec = cr.hArray_toNumpy(zoom_dynspec)
        zoom_dynspec = zoom_dynspec.reshape((zoom_dynspec.shape[0]/tbin,tbin,zoom_dynspec.shape[1]))
        zoom_dynspec = np.sum(zoom_dynspec,axis=1)
        zoom_dynspec = cr.hArray(zoom_dynspec)
    zoom_dynspec = zoom_dynspec.Transpose()
    tbin=4
    if tbin>1:
        zoom_dynspec = cr.hArray_toNumpy(zoom_dynspec)
        zoom_dynspec = zoom_dynspec[:-3]
        zoom_dynspec = zoom_dynspec.reshape((zoom_dynspec.shape[0]/tbin,tbin,zoom_dynspec.shape[1]))
        zoom_dynspec = np.sum(zoom_dynspec,axis=1)
        zoom_dynspec = cr.hArray(zoom_dynspec)
    zoom_dynspec = zoom_dynspec.Transpose()
    cr.plt.imshow(np.log10(cr.hArray_toNumpy(zoom_dynspec)),aspect='auto',origin='lower',cmap=cr.plt.cm.hot,extent=(times[0],times[-1],frequencies[0],frequencies[-1]))
    cr.plt.xlabel("+/- Time [s]")
    cr.plt.ylabel("Frequency [MHz]")

def calib_delay(substation=0):
    if substation:
    #For low-res
#        cal_delay=[-4.98046875e-09,0.0,8.00781249999e-10,-8.046875e-09,-4.21875e-09,-5.05859375e-09,7.81250000002e-10,-5.4296875e-09,-8.41796875e-09, 3.12499999998e-10,-5.1953125e-09 ,-5.3515625e-09]
    #For low-res
#        cal_delay=[0.0,-4.98046875e-09,8.00781249999e-10,-8.046875e-09,-4.21875e-09,-5.05859375e-09]
    #For multi-station
#        cal_delay=[0.0,-4.98046875e-09,-8.046875e-09,-5.05859375e-09,8.00781249999e-10,-4.21875e-09,-5.1953125e-09,7.81250000002e-10 ,-5.3515625e-09,-8.41796875e-09,3.12499999998e-10,-5.4296875e-09 ]
    #For multi-station pol1.n.pol0
#        cal_delay=[0.0,-4.98046875e-09,-8.046875e-09,-5.05859375e-09,8.00781249999e-10,-4.21875e-09,-5.1953125e-09,7.81250000002e-10 ,-5.3515625e-09,-8.41796875e-09,3.12499999998e-10,-5.4296875e-09,0.0,-5.05859375e-09,8.00781249999e-10,-4.98046875e-09,-8.046875e-09,-4.21875e-09,-8.41796875e-09,7.81250000002e-10,-5.1953125e-09,3.12499999998e-10,-5.3515625e-09,-5.4296875e-09]
    #For multi-station pol0
        cal_delay=[0.0,-5.05859375e-09,8.00781249999e-10,-4.98046875e-09,-8.046875e-09,-4.21875e-09,-8.41796875e-09,7.81250000002e-10,-5.1953125e-09,3.12499999998e-10,-5.3515625e-09,-5.4296875e-09]
    #Original Order:
    ##cal_delay=[0.0,8.00781249999e-10,-4.98046875e-09,-8.046875e-09,-5.05859375e-09,-4.21875e-09,3.12499999998e-10,7.81250000002e-10,-5.4296875e-09,-8.41796875e-09 ,-5.1953125e-09 ,-5.3515625e-09]  #wrt CS002a, , CC hFFTWResample , ears.
    else:
    # low res?
    #    cal_delay=[0.0,-1.3371875e-08,1.6859375e-09,2.3625e-09,-5.8421875e-09,2.6046875e-09]
    #off1
    #    cal_delay=[0.0,1.6859375e-09,2.3625e-09,-5.8421875e-09,-1.3371875e-08,2.6046875e-09]
    #off5
    #    cal_delay=[0.0,2.3625e-09,-5.8421875e-09,1.6859375e-09,-1.3371875e-08,2.6046875e-09]
    #high_res
        cal_delay=[-5.8421875e-09,1.6859375e-09,2.6046875e-09,0.0,2.3625e-09,-1.3371875e-08]
    #Original Order:
    #    cal_delay=[0.0,1.6859375e-09,2.3625e-09,-5.8421875e-09,2.6046875e-09,-1.3371875e-08]

    return cal_delay

#------------------------------------------------------------------------------------------------------------------------
#-------------------------------PREPARE FOR SUPER SCRIPT BELOW-----------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------

#Find beam files.
BEAMDATA=os.environ["BEAMDATA"].rstrip('/')+'/'
BEAMOUT=os.environ["BEAMOUT"].rstrip('/')+'/'

beams_dir = 'beams.sub_station/'

adding_beams=0

#Target parameters:
DM = 26.76
tbin = 1

if adding_beams==1:
    beams_suffix = '*pol0*.multi_station.beam'

    filenames = glob.glob(BEAMDATA+beams_dir+beams_suffix)

    #Open files.
    beams = cr.open(filenames)

    #Setting file parameters
    nblocks = beams['NCHUNKS']*beams['BEAM_NBLOCKS']
    block_tstep = beams['BEAM_SAMPLE_INTERVAL']*beams['BEAM_BLOCKLEN']
    substation = 1

    #Set finned tuned station calibration.
    cal_delay = calib_delay(substation=substation)
    beams['CAL_DELAY']=cal_delay

    #Add beams
    TAB = bt.addBeams(beams,dm=DM)

    #Calculate dynamic spectrum.
    dynspec,cleandynspec = bt.rawdyncalc(TAB,clean=True,tbin=tbin)

#---------------------------
#If saved binary.
elif adding_beams ==2:
    filename = BEAMDATA+'pol1.cleandynspec.bin'
    cleandynspec = cr.hArray(complex, [8193L, 764L])
    cleandynspec.readfilebinary(filename,0)
elif adding_beams ==3:
    filename_pol1 = BEAMDATA +beams_dir+ 'L43784_D20120125T211154.Superterp.pol1.TAB.bin'
    TAB_pol1 = cr.hArray(complex, [12224, 8193])
    TAB_pol1.readfilebinary(filename_pol1,0)

    filename_pol0 = BEAMDATA +beams_dir+ 'L43784_D20120125T211154.Superterp.pol0.TAB.bin'
    TAB_pol0 = cr.hArray(complex, [12224, 8193])
    TAB_pol0.readfilebinary(filename_pol0,0)

    #Calculate dynamic spectrum.
#    dynspec_pol1,cleandynspec_pol1 = bt.rawdyncalc(TAB_pol1,clean=True,tbin=tbin)
    dynspec_pol0,cleandynspec_pol0 = bt.rawdyncalc(TAB_pol0,clean=True,tbin=tbin)

    cleandynspec = cleandynspec_pol0 #+ cleandynspec_pol1

    cleandynspec.par = cleandynspec_pol0.par

else:
#    filename = BEAMDATA+beams_dir+'pol1s.multi.bin'
    filename = BEAMDATA + 'Superterp.incoherent.dd_beam.bin'
    TAB = cr.hArray(complex, [12224, 8193])
    TAB.readfilebinary(filename,0)

    #Calculate dynamic spectrum.
    dynspec,cleandynspec = bt.rawdyncalc(TAB,clean=True,tbin=tbin)


pulse_cut = 1
title = 'Superterp Stations (Coherent)-pol0'

if not pulse_cut:
    #Plot.
#    frequencies = beams['BEAM_FREQUENCIES']/10**6
#    times = cr.hArray(float,nblocks,list(np.arange(0,nblocks*block_tstep,block_tstep)))
    frequencies = cleandynspec.par.yvalues
    times = cleandynspec.par.xvalues

    cr.plt.clf()
    cr.plt.ion()
    cr.plt.imshow(np.log10(cr.hArray_toNumpy(cleandynspec)),aspect='auto',origin='lower',cmap=cr.plt.cm.hot,extent=(times[0],times[-1],frequencies[0],frequencies[-1]))
    cr.plt.rcParams['font.size']=20
    cr.plt.ylabel("Frequency [MHz]")
    cr.plt.xlabel("Time [s]")
#    cr.plt.title(beams_suffix)
    cr.plt.title(title)
else:

    pulse = 1

    # Get pulse information
    if adding_beams==1:
        hdr = beams._BeamData__files[0].par.hdr
    else:
        hdr = None
    info_pulse = getPulseInfo(pulse=pulse)
    info_pulse.par.tbin = tbin

    #Cut pulse
    zoom_dynspec=bt.cutdynspec(cleandynspec,freq_range=info_pulse.par.freq_range,time_range=info_pulse.par.time_range,hdr=hdr)
    zoom_dynspec.par.station = title#beams['STATION_NAME'][0]

    #plot
    super_plot(zoom_dynspec,info=info_pulse)


