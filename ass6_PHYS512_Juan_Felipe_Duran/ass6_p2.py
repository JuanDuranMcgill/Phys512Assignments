

import numpy as np
import matplotlib.pyplot as plt
import time 





#This function gives Lorentzian shape
def lorentz_dev(n):
    x=np.random.rand(n)
    return np.tan(np.pi*(x-0.5))



#This function has the accept ratio of Gaussia/Lorentzian
def gauss_from_lorentz(x,cons):
    #gauss over lorentzian
    accept_prob=cons*np.exp(-0.5*x**2)/(1/(1+x**2))
    
    accept=np.random.rand(len(accept_prob))<accept_prob
    return x[accept]
    

n=1000000

#Calling to get a lorenztian distribution

t1=time.time()
y=lorentz_dev(n)
z=gauss_from_lorentz(y,1/1.23)
t2=time.time()



print("accept ratio was", len(z)/len(y))



#Print histogram
a,b=np.histogram(z,200)
bb=0.5*(b[:-1]+b[1:])
plt.bar(bb,a, width= 1/bb.max())

#Print red outline of Gaussian
Normalize=np.max(a)
pred=np.exp(-0.5*bb**2)
pred=pred*Normalize
plt.plot(bb,pred,'r')

print("In terms of efficiency, the generator took", t2-t1,"seconds")


