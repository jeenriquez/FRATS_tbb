#! /usr/bin/env python
'''Pipeline 1a'''




#------------------------------



#Find file names for given observation.

LOFARDATA=os.environ["LOFARDATA"].rstrip('/')+'/'
filenames = glob.glob(LOFARDATA)

beams_dir = 'beams.sub_station/'


#------------------------------
#Find starting time and save the calibration array.
#Use MultiTBBData to open all the files.

f = cr.open([LOFARDATA+fname for fname in filenames])

#f.calcTimeStart()...need to develop this.








#------------------------------
#Post - Beamforming


#f['DM'] = dm
