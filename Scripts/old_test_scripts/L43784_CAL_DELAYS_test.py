''' L43784
'''


from pycrtools import *
import Beam_Tools as bt
from pycrtools.tasks import fitbaseline
beam_suffix='*pol0_sample*'
filenames=glob.glob(beam_suffix)
beams=open(filenames)
RFI1  = range(3605,3625)
RFI2  = range(3540,3562)
RFI3  = range(3470,3495)
RFI4  = range(3410,3430)
RFI5 = range(3350,3367)
RFI6 = range(3290,3305)
RFI7 = range(3210,3250)
RFI8 = range(3165,3185)
RFI9 = range(3100,3130)
RFI10 = range(3086,3089)
RFI11 = range(3045,3065)
RFI12 = range(2980,3000)
RFI13 = range(2925,2943)
RFI14 = range(2870,2890)
RFI15 = range(2820,2835)
RFI16 = range(2804,2808)
RFI17 = range(2764,2778)
RFI18 = range(2710,2724)
beams['RFI_CHANNELS']=RFI1
beams['RFI_CHANNELS']=RFI2
beams['RFI_CHANNELS']=RFI3
beams['RFI_CHANNELS']=RFI4
beams['RFI_CHANNELS']=RFI5
beams['RFI_CHANNELS']=RFI6
beams['RFI_CHANNELS']=RFI7
beams['RFI_CHANNELS']=RFI8
beams['RFI_CHANNELS']=RFI9
beams['RFI_CHANNELS']=RFI10
beams['RFI_CHANNELS']=RFI11
beams['RFI_CHANNELS']=RFI12
beams['RFI_CHANNELS']=RFI13
beams['RFI_CHANNELS']=RFI14
beams['RFI_CHANNELS']=RFI15
beams['RFI_CHANNELS']=RFI16
beams['RFI_CHANNELS']=RFI17
beams['RFI_CHANNELS']=RFI18
freq_max = 144.8
freq_min = 132.8
beams['DM']=26.76
pulse_time_range=[0.592,0.598]
ST_DELAYS = bt.ccBeams(beams,freq_range=[freq_min,freq_max],time_range=pulse_time_range,verbose=1)
freq_range2=[151,169.5]
pulse_time_range2 = [.021,.0255]
ST_DELAYS2 = bt.ccBeams(beams,freq_range=freq_range2,time_range=pulse_time_range2,verbose=1)
plt.ion()
plt.plot(ST_DELAYS2.vec())
ST_DELAYS2
plt.clf()
plt.plot(ST_DELAYS2.vec())
plt.plot(ST_DELAYS.vec())
plt.clf()
plt.plot(ST_DELAYS2.vec(),label='pulse top')
plt.plot(ST_DELAYS.vec(),label='pulse bottom')
plt.legend()
plt.ylabel("Time Delay [s]")
beams
beams['STATION_NAME']
sorted(beams.keys())
beams['BEAM_ANTENNA_SET']
beams['ANTENNA_SET']
plt.xticks(range(12), [beams['STATION_NAME'][i]+beams['ANTENNA_SET'][i] for i in range(12)], rotation='vertical')
plt.title("L43784, superterp")
plt.title("L43784, superterp (masked bottom pulse)")
