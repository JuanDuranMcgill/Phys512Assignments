# -*- coding: utf-8 -*-
"""
a) Pb206 is stable while U238 decays. Therefore,  Their ratio is roughly c/e^-d,
It's simply going to be a function that increases exponentially as  U238 decays.

"""

import numpy as np
from matplotlib import pyplot as plt

#PART A)
x = np.linspace(0,20000,1000)

t_half=4458

#U238:
u=(1/2)*np.exp(-x/t_half)
#Pb206:
p = 1

#ratio between Pb and U:
u1 = p/u
###############################################################################
#part B)

'''
b) For this part the ratio is roughly e^-c/e^-d and since the one above is smaller,
this is simply going to look like a linear negative slope
'''
x2= np.linspace(0,20000,1000)

t_half2=245500
t_half3=75380


#U234
uu = (1/2)*np.exp(-x/t_half2)


#THO230
t = (1/2)*np.exp(-x/t_half3)

#ratio between THO230 and U234:
u2 = t/uu


plt.figure(1)
plt.title('Pb206/U238 and THO230/U234')

plt.subplot(211)
plt.plot(x, u1)
plt.subplot(212)
plt.plot(x2,u2)
