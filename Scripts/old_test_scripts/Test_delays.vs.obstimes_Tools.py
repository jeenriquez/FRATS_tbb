''' This script reads the clockcorrection from parsests and plots its evolution with time. (HBA and LBA)
'''

import os
import pycrtools as cr
from pycrtools import metadata as md
from pycrtools import bfdata as bf
import pdb;# pdb.set_trace()
from pycrtools import tools

#dirc = '/vol/astro/lofar/frats/bf/parsets/'
dirc = '/scratch/lofar/frats/parsets/'
filenames=os.listdir(dirc)[1:]
filenames.insert(0,'L43784.parset')
stations  = ['CS002','CS003','CS004','CS005','CS006','CS007']
antenna_set = 'HBA'


cr.plt.ion()
cr.plt.clf()
cr.plt.rcParams['font.size']=20

for ii,station in enumerate(stations):
    times =[]
    obs_time = []
    times_minus = []
    print 'station : ' , station
    count = 0
    for i,file in enumerate(filenames):
        try:
            a = bf.get_parameters_new(dirc+file,useFilename=True)
            obs_time.append(tools.strdate2jd(a['starttime']))
            times.append(md.getClockCorrectionParset(dirc+file,station,antenna_set))       
    #        times.par.xvalues = obs_time    
            times_minus.append(md.getClockCorrectionParset(dirc+file,'CS002',antenna_set))       
            times_minus[count] = float(times_minus[count]) - float(times[count])
            if not count:
                tm = times_minus[0]
                t0=obs_time[0]
            times_minus[count] = float(times_minus[count]) - float(tm)    
            count+=1
        except :
            while len(times) > len(times_minus):
                times = times[:-1]
            while len(obs_time) > len(times):
               obs_time = obs_time[:-1]            
            pass

    cr.plt.plot(obs_time,times_minus,'o',ms=10,label=station)
#    pdb.set_trace()    
    
jun8 = 2456086.500   # In JD
cr.plt.plot([jun8,jun8],[-3e-6,3e-6],'r',label='June 8')
cr.plt.plot([t0,t0],[-3e-6,3e-6],'b--',label='Pulsar obs Jan 25')
#cr.plt.plot([t0,t0],[-3e-6,3e-6],'w--')
cr.plt.legend(loc=3)
cr.plt.ylabel("Delta Delta_t to station CS002 [sec]")
cr.plt.xlabel("Obs Time [JD]")
cr.plt.title('Station Delays wrt CS002: '+antenna_set)
