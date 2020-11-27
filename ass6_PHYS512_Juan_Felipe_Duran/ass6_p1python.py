

import numpy as np
import matplotlib.pyplot as plt



#Functions to populate array of rnandom numbers

def get_rands_nb(vals):
    n=len(vals)
    
    for i in range(n):
        vals[i]=1e8*np.random.rand()
    return vals

def get_rands(n):
    vec=np.empty(n,dtype='int32')
    get_rands_nb(vec)
    return vec


#Samples
n=30000

#Array of random numbers
vec=get_rands(n*3)

#Shape to 3 columns
vv=np.reshape(vec,[n,3])

x=vv[:,0]
y=vv[:,1]
z=vv[:,2]

#Plot in 3d

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(x, y, z, marker='.')

ax.view_init(50, 50)

print("For the python random numbers, I was not able to do find planes. Also, I was not able to perform this with my local machine")