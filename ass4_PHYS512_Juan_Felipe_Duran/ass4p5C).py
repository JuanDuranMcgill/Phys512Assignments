import numpy as np
from matplotlib import pyplot as plt
N=1024
x=np.arange(N)
k=15.4
y=np.sin(2*np.pi*x*k/N)
yft=np.fft.fft(y)



#Plotting sin with leakage
plt.plot(np.abs(yft),'.')