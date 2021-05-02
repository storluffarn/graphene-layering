
import numpy as np
from scipy.signal import argrelextrema

traw = np.loadtxt(open("time.csv", "rb"), delimiter=",", skiprows=1)
xraw = np.loadtxt(open("xout.csv", "rb"), delimiter=",", skiprows=1)
qraw = np.loadtxt(open("qout.csv", "rb"), delimiter=",", skiprows=1)
fraw = np.loadtxt(open("tomout.csv", "rb"), delimiter=",", skiprows=1)
times = traw[:,0]
xpos = xraw[:,0]
xvel = xraw[:,1]
xacc = xraw[:,2]
qpos = qraw[:,0]
qvel = qraw[:,1]
qacc = qraw[:,2]
forces = fraw[:,2]

maxima = argrelextrema(forces, np.greater)
minima = argrelextrema(forces, np.less)

print (maxima)
print (minima)

print (forces[minima[0][3]])
print (times[minima[0][3]])
print (xpos[minima[0][3]])
print (xvel[minima[0][3]])
print (xacc[minima[0][3]])
print (qpos[minima[0][3]])
print (qvel[minima[0][3]])
print (qacc[minima[0][3]])

