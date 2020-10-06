# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 13:29:24 2020

@author: Utilisateur
H



Here is my attempt at problem 2. I fit a spline to the data and then compared
it to a linear fit in order to get an estimate of the error. I got the max 
error that way. 
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate



#load file
data = np.loadtxt("lakeshore.txt")


#load x and y values, make temperature y axis and voltage on x axis
y=data[:,0]
x=data[:,1]

y=y[::-1]
x=x[::-1]
ye = y[::2]
xe = x[::2]

#When interpolating use this linspace
xx = np.linspace(x[0],x[-1],2001)

#Interpolating part to get y values
spln=interpolate.splrep(x,y)
yy=interpolate.splev(xx,spln)

#Interpolating with spline that takes every other element
ye = y[::2]
xe = x[::2]
splne=interpolate.splrep(xe,ye)
yye=interpolate.splev(xx,splne)




#Plotting plot
plt.subplot(2, 1, 1)
plt.plot(xx,yy, color='green')
plt.ylabel('Temperature(K)')


#Plotting error estimate by comparing spline with spline of everyother element
plt.subplot(2, 1, 2)
plt.plot(xx,yy-yye)
plt.ylabel('Error')
plt.xlabel('Voltage(V)')

#Printing max error
print('the max error is ' , np.max(np.abs(yy-yye)))


#Function that returns temperature for input voltage
def VoltageToTemp(x):
    y = interpolate.splev(x,spln)
    return y

#Testing the function with a print
voltageInput = 0.8
print('the temp of ', voltageInput,'V is ', VoltageToTemp(voltageInput), 'K')





