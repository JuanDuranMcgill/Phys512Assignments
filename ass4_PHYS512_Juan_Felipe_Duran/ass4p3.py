import numpy
from matplotlib import pyplot as plt


x = numpy.arange(-10,10,0.1)
y = numpy.exp(-0.5*x**2/(1.5**2))


def convShift(array,amount):
    N =array.size
    kvec = numpy.arange(N)
    yft = numpy.fft.fft(array)
    J = numpy.complex(0,1)
    dx = amount
    yft_new=yft*numpy.exp(-2*numpy.pi*J*kvec*dx/N)
    y_new=numpy.real(numpy.fft.ifft(yft_new))
    
    return y_new



def corr(f,g):
    ft1 = numpy.fft.fft(f)
    ft2 = numpy.fft.fft(g)
    return numpy.real(numpy.fft.ifft(ft1)*numpy.conj(numpy.fft.ifft(ft2)))

    
plt.plot(x,y)
#plt.plot(x,corr(y,y))
plt.plot(x,convShift(y,25.0), 'r')
plt.plot(x,corr(y,convShift(y,25.0)))

"""
How does the correlation depend on the shift?
-The correlation only happens when functions overlap each other
so the more shifter away a function is, the less correlation there will be. 