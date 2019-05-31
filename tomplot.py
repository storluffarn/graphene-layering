
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

traw = np.loadtxt(open("time.csv", "rb"), delimiter=",", skiprows=1)
xraw = np.loadtxt(open("xout.csv", "rb"), delimiter=",", skiprows=1)
qraw = np.loadtxt(open("qout.csv", "rb"), delimiter=",", skiprows=1)
fraw = np.loadtxt(open("tomout.csv", "rb"), delimiter=",", skiprows=1)

nscale = 1e9

tdata = nscale * traw[:,0]
sdata = nscale * traw[:,1]
xdata = nscale * xraw[:,0]
qdata = nscale * qraw[:,0]
fdata = nscale * fraw[:,2]

fig0, ax = plt.subplots()
ax.plot(tdata,sdata)
fig0.savefig("out0.png")

fig1, ax = plt.subplots()
ax.plot(tdata,xdata, color='#ff7f0e')
ax.set(xlabel='t (nm)', ylabel='x (nm)')
ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
ax.xaxis.grid(True, which='minor',linestyle='dotted')
ax.yaxis.grid(True, linestyle='dotted')
fig1.savefig("out1.png")

fig2, ax = plt.subplots()
ax.plot(tdata,qdata)
ax.set(xlabel='t (nm)', ylabel='q (nm)')
fig2.savefig("out2.png")

fig3, ax = plt.subplots()
ax.plot(tdata,fdata)
ax.set(xlabel='t (nm)', ylabel='F (nN)')
fig3.savefig("out3.png")

