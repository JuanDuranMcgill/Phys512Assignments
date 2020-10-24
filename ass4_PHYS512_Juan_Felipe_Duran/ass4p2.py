import numpy
from matplotlib import pyplot as plt


x = numpy.arange(-10,10,0.1)
y = numpy.exp(-0.5*x**2/(1.5**2))


#f*g = ift(dft(f)*conj(dft(g))).

def corr(f,g):
    ft1 = numpy.fft.fft(f)
    ft2 = numpy.fft.fft(g)
    return numpy.real(numpy.fft.ifft(ft1)*numpy.conj(numpy.fft.ifft(ft2)))

plt.plot(x,y)

plt.plot(x,corr(y,y))