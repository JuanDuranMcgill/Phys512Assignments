
import numpy as np
#from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.pyplot as plt



#Open rand_points.txt to get values from C 
f=open("rand_points.txt","r")
lines=f.readlines()
resultx=[]
resulty=[]
resultz=[]
for x in lines:
    resultx.append(int(x.split(' ')[0]))
    resulty.append(int(x.split(' ')[1]))
    resultz.append(int(x.split(' ')[2]))
f.close()




points=[]
counter = 0

for r in resultx:
    row = []

    row.append(resultx[counter])
    row.append(resulty[counter])
    row.append(resultz[counter])
    
    counter += 1
    points.append(row)
    


#Plot in 3d
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(resultx, resulty, resultz, marker='.'
           )

ax.view_init(0, 50)
plt.ion()

print("By inspection, I was able to count 30 planes, looking at the left part of the figure.")