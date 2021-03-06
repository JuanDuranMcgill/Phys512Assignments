import numpy as np
from matplotlib import pyplot as plt
N=1024
x=np.arange(N)
k=15.4
y=np.sin(2*np.pi*x*k/N)
yft=np.fft.fft(y)



#Removing leakage from sin
window = 0.5-0.5*np.cos((2*np.pi*x)/N)


yft2=np.fft.rfft(y*window)


plt.plot(np.abs(yft2),'.')
