import numpy as np
import camb
from matplotlib import pyplot as plt

#from scipy.stats import chisquare


def get_spectrum(pars,lmax=1199):
    #print('pars are ',pars)
    H0=pars[0]
    ombh2=pars[1]
    omch2=pars[2]
    tau=0.05
    As=pars[3]
    ns=pars[4]
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
plt.clf()



wmap=np.loadtxt('wmap_tt_spectrum_9yr_v5.txt')

pars=np.asarray([65,0.02,0.1,2e-9,0.96])
mod = get_spectrum(pars)
model = mod[2:]
y = wmap[:,0]

print (len(model))
print(len(y))





def num_deriv(fun,x,pars,dpar):
    #calculate numerical derivatives of 
    #a function for use in e.g. Newton's method or LM
    derivs=np.zeros([len(x),len(pars)])
    for i in range(len(pars)):
        pars2=pars.copy()
        pars2[i]=pars2[i]+dpar[i]
        f_right=fun(pars2)
        f_right=f_right[2:]
        pars2[i]=pars[i]-dpar[i]
        f_left=fun(pars2)
        f_left=f_left[2:]
        derivs[:,i]=(f_right-f_left)/(2*dpar[i])
    return derivs


noise = 0.1

Ninv=np.eye(len(model))/noise**2
#pars=np.asarray([65,0.02,0.1,2e-9,0.96])
dpar=np.asarray([0.65,0.0002,0.001,2e-11,0.0096])
for i in range(10):
    mod = get_spectrum(pars)
    model = mod[2:]
    derivs=num_deriv(get_spectrum,model,pars,dpar)
    resid=y-model
    lhs=derivs.T@Ninv@derivs
    rhs=derivs.T@Ninv@resid
    lhs_inv=np.linalg.inv(lhs)
    step=lhs_inv@rhs
    pars=pars+step
    print(pars)

par_errs=np.sqrt(np.diag(np.linalg.inv(lhs)))
print('final parameters are ',pars,' with errors ',par_errs)



























