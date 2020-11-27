
import numpy as np
import matplotlib.pyplot as plt
import time



#U FROM 0 TO 1, V FROM -1 TO 1
n=1000000
u=np.random.rand(n)
v=(2*np.random.rand(n)-1)*0.85

#Here is the ratio
rat=v/u

t1=time.time()
#Here we accept u< sqrt(p(v/u))
accept=u<np.sqrt(np.exp(-0.5*rat**2))
#Output Gauss
mygauss=rat[accept]
t2=time.time()


print("accept fraction is ", np.mean(accept))

mygauss.sort()



#Plot histogram
a,b=np.histogram(mygauss,200)
bb=0.5*(b[:-1]+b[1:])
plt.bar(bb,a, width= 1/bb.max())


#Plot Gaussian red outline
Normalize=np.max(a)
pred=np.exp(-0.5*bb**2)
pred=pred*Normalize
plt.plot(bb,pred,'r')

print("In terms of efficiency, the generator took", t2-t1,"seconds")
