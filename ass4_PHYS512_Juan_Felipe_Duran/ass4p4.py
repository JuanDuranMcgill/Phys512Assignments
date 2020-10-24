import numpy as np
from matplotlib import pyplot as plt

x=np.arange(1000)
tau=50
nhit=20
x_hit=np.asarray(np.floor(len(x)*np.random.rand(nhit)),dtype='int')
y_hit=np.random.rand(nhit)**2

f=0.0*x
for i in range(nhit):
    f[x_hit[i]]=f[x_hit[i]]+y_hit[i]
g=np.exp(-1.0*x/tau)


def paddedConvolution(a1,a2):
    
    a1=np.pad(a1,[0,len(a1)],'constant')
    a2=np.pad(a2,[0,len(a2)],'constant')
    conv=np.fft.irfft(np.fft.rfft(a1)*np.fft.rfft(a2))
    
    return conv

    
plt.plot(paddedConvolution(f,g))
