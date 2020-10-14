import numpy as np
import camb
from matplotlib import pyplot as plt
import time
import scipy
from scipy.stats import chisquare


def get_spectrum(pars,lmax=1199):
    #print('pars are ',pars)
    H0=pars[0]
    ombh2=pars[1]
    omch2=pars[2]
    tau=pars[3]
    As=pars[4]
    ns=pars[5]
    pars=camb.CAMBparams()
    pars.set_cosmology(H0=H0,ombh2=ombh2,omch2=omch2,mnu=0.06,omk=0,tau=tau)
    pars.InitPower.set_params(As=As,ns=ns,r=0)
    pars.set_for_lmax(lmax,lens_potential_accuracy=0)
    results=camb.get_results(pars)
    powers=results.get_cmb_power_spectra(pars,CMB_unit='muK')
    cmb=powers['total']
    tt=cmb[:,0]    #you could return the full power spectrum here if you wanted to do say EE
    return tt


plt.ion()

pars=np.asarray([65,0.02,0.1,0.05,2e-9,0.96])
wmap=np.loadtxt('wmap_tt_spectrum_9yr_v5.txt')

plt.clf();
#plt.errorbar(wmap[:,0],wmap[:,1],wmap[:,2],fmt='*')
#plt.plot(wmap[:,0],wmap[:,1],'.')

cmb=get_spectrum(pars)

#plt.plot(cmb)
cmb1 = cmb[2:]
plt.plot(cmb1)

rms=np.std(wmap[:,0]-cmb1)
N=rms**2
chisq=np.sum((wmap[:,0]-cmb1)**2)/N**2
print("RMS scatter about fit is ",rms)

print(chisq)
