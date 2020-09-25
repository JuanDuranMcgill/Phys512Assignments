'''
ass2p2

I was able to get a fit with error less than 10^-6
with order 10, # coefficients = 8. 


'''


import numpy as np
from matplotlib import pyplot as plt


def cheb_mat_uniform(nx,ord): # This function inputs number of points and order
                             #It outputs the Chebyshev matrix and x
    x=np.linspace(-1,1,nx)
    mat=np.zeros([nx,ord+1])
    mat[:,0]=1.0
    if ord>0:
        mat[:,1]=x
    if ord>1:
        for i in range(1,ord):
            mat[:,i+1]=2*x*mat[:,i]-mat[:,i-1]
    return mat,x


#Here I call chebyshev "Hey chebyshev come right here!!"
n=200
ord=10
mat,x=cheb_mat_uniform(n,ord)


#Here is the function, I shifted x for (x+3)/4 to get values between 0.5 and 1
y=np.log2((x+3)/4)

#Manipulating the matrix with y to get the fit
lhs=np.dot(mat.transpose(),mat)
rhs=np.dot(mat.transpose(),y)
fitp=np.dot(np.linalg.inv(lhs),rhs)


#Get the prediction with given number of coefficients
ncoeff=8
pred=np.dot(mat[:,:ncoeff],fitp[:ncoeff])


#Shifting x to get values btween 0.5 and 1
x = (x+3)/4

#Plotting the predicted values and the real values side by side 
plt.plot(x,pred)
plt.plot(x,y)



print('rms error is ',np.sqrt(np.mean((pred-y)**2)),' with max error ',np.max(np.abs(pred-y)))
