

from pycrtools import *
import Beam_Tools as bt
from pycrtools.tasks import fitbaseline
beam_suffix='*CS002*pol0_sample*'
filenames=glob.glob(beam_suffix)
beams=open(filenames)
beams['DM']=26.76

TAB,dynspec,cleandynspec = bt.addBeams(beams,dyncalc=True,clean=True)
zoom_dynspec = bt.cutDynspec(dynspec,time_range=[0.0,0.2])
flat_dynspec = bt.flatDynspec(zoom_dynspec,axis='y',verbose=True)


#----

#-----------------------------------

flat_numpy_dynspec=hArray(float,flat_dynspec,flat_dynspec)
flat_numpy_dynspec =flat_numpy_dynspec.toNumpy()
x_numpy_values = flat_dynspec.par.xvalues.toNumpy()
numpy_coeffs = np.polyfit(x_numpy_values,flat_numpy_dynspec,10)
numpy_fit = np.polyval(numpy_coeffs,x_numpy_values)

flat_fit = hArray(numpy_fit)
flat_fit.sqrt()
flat_fit=1/flat_fit


for block in range(12480):
   TAB2[block] = TAB[block]*flat_fit


zoom_dynspec2 = bt.cutDynspec(dynspec2,time_range=[0.0,0.2])
flat_dynspec2 = bt.flatDynspec(zoom_dynspec2,axis='y',verbose=True)



#----------------------#
#Talked with Pim;
#   Then use 1/sqrt(fit) since the amplitude is X**2+Y**2 (at least for now)
#   I should also convolved the pulsar signal with something (like gausian?), to removed the ringing.
#   The real frequency gain calibration is done by assuming the galaxy in the case of the LBAs, and the spectrum of a stron source for the HBAs.

#----------------------#
stop
#-----------------------------------





uni_dynspec=hArray(float,[8193],flat_dynspec/flat_dynspec)
uni_var = hArray(float,[1,5])
uni_coeff = hArray(float,[1,5])

uni_xpowers=hArray(float,[1,256,5])
powers=hArray(int,[1,5],range(5))
ndata=-1



hLinearFit(uni_coeff, uni_var, xvec, yvec, ndata)
hLinearFitPolynomialX(uni_xpowers, xvec, powers)


#self.xpowers[..., [0]:self.nselected_bins].linearfitpolynomialx(self.clean_bins_x[..., [0]:self.nselected_bins], self.powers[...])
#self.chisquare = self.coeffs[...].linearfit(self.covariance[...], self.xpowers[...], self.clean_bins_y[...], self.nselected_bins)  # self.weights[...],

#self.fitted_bins_y[...].polynomial(self.fitted_bins_x, self.coeffs[...], self.powers[...])



#-----------------------------------
