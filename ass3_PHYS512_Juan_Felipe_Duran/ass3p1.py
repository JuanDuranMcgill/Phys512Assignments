import sys
import csv
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from scipy.optimize import curve_fit

# Read CSV
csvFileName = "dish_zenith.csv"
csvData = []
with open(csvFileName, 'r') as csvFile:
    csvReader = csv.reader(csvFile, delimiter=',')
    for csvRow in csvReader:
        csvData.append(csvRow)

# Get X, Y, Z

csvData = np.array(csvData)
csvData = csvData.astype(np.float)
X, Y, Z = csvData[:,0], csvData[:,1], csvData[:,2]

# Plot X,Y,Z
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_trisurf(X, Y, Z, color='white', edgecolors='grey', alpha=0.5)
ax.scatter(X, Y, Z, c='red')
plt.show()


def func(X, a, x0, y0, z0):
    x,y = X
    return a*((x-x0)**2+(y-y0)**2)+z0

params,hey = curve_fit(func,(X,Y),Z)
perr = np.sqrt(np.diag(hey))
print("parameters are", params)
print("errors are", perr)
    