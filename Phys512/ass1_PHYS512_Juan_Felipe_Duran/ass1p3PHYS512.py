# -*- coding: utf-8 -*-

"""
Problem 3:
    The error of the rational function for the Lorentzian should be small. this 
    is because Lorentzian is a rational function. However, when we set it to
    n=4 and m=5, the rational does not behave as well as expected.  
    When we use pinv instead of pinv, we actually calculate the pseudo inverse.
    This will output a pseudo inverse, even for matrices that do not have inverses.
    As a consequence, the values of p and q are very much smaller which allows the function
    to behave better even at higher inputs of n and m. 
"""


import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import math 
#definition the rational function stuff
def rat_eval(p,q,x):
    top=0
    for i in range(len(p)):
        top=top+p[i]*x**i
    bot=1
    for i in range(len(q)):
        bot=bot+q[i]*x**(i+1)
    return top/bot

def rat_fit(x,y,n,m):
    assert(len(x)==n+m-1)
    assert(len(y)==len(x))
    mat=np.zeros([n+m-1,n+m-1])
    for i in range(n):
        mat[:,i]=x**i
    for i in range(1,m):
        mat[:,i-1+n]=-y*x**i
    pars=np.dot(np.linalg.inv(mat),y)
    p=pars[:n]
    q=pars[n:]
    return p,q

#setting values for n and m
n=3
m=4

#Setting up the Cos
pi = math.pi
x = np.linspace(-pi/2,pi/2,n+m-1)
y =np.cos(x)
xx = np.linspace(x[0],x[-1],1000)
y_true = np.cos(xx)

#Setting up the rational function
p,q =rat_fit(x,y,n,m)
print('values for p and q are',p,q)
pred = rat_eval(p,q,xx)

#Setting up the polyfit
fitp =np.polyfit(x,y,n+m-1)
pred_poly = np.polyval(fitp,xx)

#Setting up the spline
spln=interpolate.splrep(x,y)
pred_spln=interpolate.splev(xx,spln)


#Printing all three
plt.plot(x,y, '*')
plt.plot(xx,pred)
plt.plot(xx,pred_poly)
plt.plot(xx,pred_spln)

#Printing error estimates
print('rat agrees to',np.std(y_true-pred)*100,'%')
print('poly agrees to',np.std(y_true-pred_poly)*100,'%')
print('spln agrees to',np.std(y_true-pred_spln)*100,'%')