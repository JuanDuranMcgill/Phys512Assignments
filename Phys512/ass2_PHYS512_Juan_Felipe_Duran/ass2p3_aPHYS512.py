  

"""
For this problem I would use implicit with radau because it deals with big intervals 
and small intervals

"""

import numpy as np
from scipy import integrate
import time

def fun(x,y,half_life=[1e9,24.10/365,6.70/(24*365),245500,75380,1600,
                       3.8235/365,3.10/(60*24*365),26.8/(60*24*365),
                       19.9/(60*24*365),
                       (164.3*10e-6)/(60*60*24*365),
                       22.3,5.015,
                       138.376/365]):
    
    dydx=np.zeros(len(half_life)+len(half_life)-1)
    dydx[0]=-y[0]/half_life[0]
    dydx[1]=y[0]/half_life[0]-y[1]/half_life[1]
    dydx[2]=y[1]/half_life[1]
    dydx[3]=y[1]/half_life[1]-y[2]/half_life[2]
    dydx[4]=y[2]/half_life[2]
    dydx[5]=y[2]/half_life[2]-y[3]/half_life[3]
    dydx[6]=y[3]/half_life[3]
    dydx[7]=y[3]/half_life[3]-y[4]/half_life[4]
    dydx[8]=y[4]/half_life[4]
    dydx[9]=y[4]/half_life[4]-y[5]/half_life[5]
    dydx[10]=y[5]/half_life[5]
    dydx[11]=y[5]/half_life[5]-y[6]/half_life[6]
    dydx[12]=y[6]/half_life[6]
    dydx[13]=y[6]/half_life[6]-y[7]/half_life[7]
    dydx[14]=y[7]/half_life[7]
    dydx[15]=y[7]/half_life[7]-y[8]/half_life[8]
    dydx[16]=y[8]/half_life[8]
    dydx[17]=y[8]/half_life[8]-y[9]/half_life[9]
    dydx[18]=y[9]/half_life[9]
    dydx[19]=y[9]/half_life[9]-y[10]/half_life[10]
    dydx[20]=y[10]/half_life[10]
    dydx[21]=y[10]/half_life[10]-y[11]/half_life[11]
    dydx[22]=y[11]/half_life[11]
    dydx[23]=y[11]/half_life[11]-y[12]/half_life[12]
    dydx[24]=y[12]/half_life[12]
    dydx[25]=y[12]/half_life[12]-y[13]/half_life[13]
    dydx[26]=y[13]/half_life[13]
    return dydx


y0=np.asarray([1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]) 
x0=0
x1=1
t1=time.time();
ans_rk4=integrate.solve_ivp(fun,[x0,x1],y0);
t2=time.time();
print('took ',ans_rk4.nfev,' evaluations and ',t2-t1,' seconds to solve with RK4.')
t1=time.time()
ans_stiff=integrate.solve_ivp(fun,[x0,x1],y0,method='Radau')
t2=time.time()
print('took ',ans_stiff.nfev,' evaluations and ',t2-t1,' seconds to solve implicitly')
print('final values were ',ans_rk4.y[0,-1],' and ',ans_stiff.y[0,-1])
