import numpy as np



#moon_astero_dist= distance(abs(decimal_RA(moon2[:,1:4]) - decimal_RA(RA_DEC[:,1:4])),abs(decimal_DEC(moon2[:,4:7]) - decimal_DEC(RA_DEC[:,4:7])))



#------------------------------------------
def decimal_RA(RA):
# Function to calculate the RA decimal value from degres/arcmin/arcsec .

    import sys

    if len(RA[0,:]) < 3 :
        print 'Dimmentions dont match.'
        sys.exit()

    RA_deci = RA[:,0]*15. + RA[:,1]/60.	+ RA[:,2]/(60.*60.)
    return RA_deci


#------------------------------------------
def decimal_DEC(DEC):
# Function to calculate the DEC decimal value from degres/arcmin/arcsec .

    import sys

    if len(DEC[0,:]) < 3 :
        print 'Dimmentions dont match.'
        sys.exit()

    DEC_deci = DEC[:,0] + DEC[:,1]/60.	+ DEC[:,2]/(60.*60.)
    return DEC_deci


#------------------------------------------
def distance(x,y):
# Function to calculate distance between elements of two vectors.

    from math import sqrt

    z = x*0.
    for i in xrange(len(x)):
        z[i] = sqrt(x[i]**2. + y[i]**2.)

    return z



