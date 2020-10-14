'''
ass2p2

I was able to get a fit with error less than 10^-6
with order 10, # coefficients = 8. 


'''


import numpy as np
from matplotlib import pyplot as plt






#Here I call chebyshev "Hey chebyshev come right here!!"
#Here is the linspacee
x = np.linspace(-1,1,200)
#Here is the function, I shifted x for (x+3)/4 to get values between 0.5 and 1
y=np.log2((x+3)/4)

#Here I fit it like Legendre 
fitLegendre=np.polynomial.legendre.legfit(x,y,10)
predLegendre=np.polynomial.legendre.legval(x,fitLegendre)





#Shifting x to get values btween 0.5 and 1
x = (x+3)/4

#Plotting the predicted values and the real values side by side 

plt.plot(x,y)
plt.plot(x,predLegendre)

#Printing errors

print('Legendre rms error is ',np.sqrt(np.mean((predLegendre-y)**2)),' with max error ',np.max(np.abs(predLegendre-y)))
