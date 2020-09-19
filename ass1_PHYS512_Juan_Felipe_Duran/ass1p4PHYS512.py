# -*- coding: utf-8 -*-
"""

Problem 4: There was a singularity in the problem, but both integrators seemed
to be able to handle it properly
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import constants as c
from scipy import integrate
import math 



#Define some constants
e = c.epsilon_0
pi = math.pi
K_e = 1/(4*pi*e)
R = 10
dens = 1



#Define the in class integrator
def integrate_step(fun,x1,x2,tol):
    print('integrating from ',x1,' to ',x2)
    x=np.linspace(x1,x2,5)
    y=fun(x)
    area1=(x2-x1)*(y[0]+4*y[2]+y[4])/6
    area2=(x2-x1)*( y[0]+4*y[1]+2*y[2]+4*y[3]+y[4])/12
    myerr=np.abs(area1-area2)
    if myerr<tol:
        return area2
    else:
        xm=0.5*(x1+x2)
        a1=integrate_step(fun,x1,xm,tol/2)
        a2=integrate_step(fun,xm,x2,tol/2)
        return a1+a2



xx=np.linspace(0,50,50)



def F(x):
    res = np.zeros_like(x)
    #This loop prints the output of the integral for different values of z
    for i,val in enumerate(x):
        z=val
        #This is the function of the integration of the electric field of a spherical shell
        def f(x):    
            return K_e*(2*pi*R**2*dens)*((z-R*np.cos(x))*np.sin(x))/((R**2+z**2-2*R*z*np.cos(x))**(3/2))
        #Call quad
        y,err = integrate.quad(f,0,pi)
        res[i]=y
    return res

def Fclass(x):
    res = np.zeros_like(x)
    #This loop prints the output of the integral for different values of z
    for i,val in enumerate(x):
        z=val
        #This is the function of the integration of the electric field of a spherical shell
        def f(x):    
            return K_e*(2*pi*R**2*dens)*((z-R*np.cos(x))*np.sin(x))/((R**2+z**2-2*R*z*np.cos(x))**(3/2))
        #Call quad
        y = integrate_step(f,0,pi,10)
        res[i]=y
    return res

#Plot quad
plt.plot(xx,F(xx))

#Plot class integrator in stars
plt.plot(xx,Fclass(xx), '*')








