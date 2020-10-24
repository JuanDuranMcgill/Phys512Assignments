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
    
plt.plot(x,y)
plt.plot(x,convShift(y,100.0), 'r')