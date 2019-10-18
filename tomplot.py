
import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

kb = 1.38064852e-23
xmass = 1e-23
xdamp = 1.875e+13
T = 250
tstep = 1.5e-14
alltime  = 1.2e-08

traw = np.loadtxt(open("time.csv", "rb"), delimiter=",", skiprows=1)
xraw = np.loadtxt(open("xout.csv", "rb"), delimiter=",", skiprows=1)
qraw = np.loadtxt(open("qout.csv", "rb"), delimiter=",", skiprows=1)
fraw = np.loadtxt(open("tomout.csv", "rb"), delimiter=",", skiprows=1)
sraw = np.loadtxt(open("slips.csv", "rb"), delimiter=",", skiprows=0)
avgraw = np.loadtxt(open("avgs.csv", "rb"), delimiter=",", skiprows=1)

nscale = 1.0e9  # WARNING this is *only* for plotting, danger!

tdata = nscale * traw[:,0]
sdata = nscale * traw[:,1]
xdata = nscale * xraw[:,0]
qdata = nscale * qraw[:,0]
fdata = nscale * fraw[:,2]
slips = nscale * sraw[:,1]
slipst = nscale * sraw[:,0]
avgs  = nscale * avgraw[:,1]
avgst  = nscale * avgraw[:,0]

xaccdata = xraw[:,2]

fig0, ax = plt.subplots()
ax.plot(tdata,sdata)
fig0.savefig("plots/plot_vt.png")

fig1, ax = plt.subplots()
ax.plot(tdata,xdata)
ax.set(xlabel='t (nm)', ylabel='x (nm)')
ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
ax.xaxis.grid(True, which='minor',linestyle='dotted')
ax.yaxis.grid(True, linestyle='dotted')
fig1.savefig("plots/plot_x.png")

fig2, ax = plt.subplots()
ax.plot(tdata,qdata)
ax.set(xlabel='t (nm)', ylabel='q (nm)')
fig2.savefig("plots/plot_q.png")

fig3, ax = plt.subplots()
ax.plot(tdata,fdata)
ax.set(xlabel='t (nm)', ylabel='x (nm)')
ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
ax.xaxis.grid(True, which='minor',linestyle='dotted')
ax.yaxis.grid(True, linestyle='dotted')
ax.set(xlabel='t (nm)', ylabel='F (nN)')
fig3.savefig("plots/plot_f.png")

# for debugging

fig4, ax = plt.subplots()

nbins = 100
n,f1,patches = ax.hist(xdata,nbins,edgecolor='black')

fig4.savefig("plots/histogram.png")

# mean square displacement

#xdata /= nscale 
#mean = 0
#
#for x in xdata :
#    mean += x**2
#
#mean /= len(xdata)
#
#print("simulated msd: {}".format(mean))
#
## analytical
#
#diff = kb*T/(xmass*xdamp)
#
#msd = 2*diff*alltime
#
#c = 2*xdamp*kb*T/xmass
#
#msd2 = c/(2*xdamp)*(1-math.exp(-xdamp*alltime)**2)/xdamp**2 + c*alltime/xdamp**2 - c*(1-math.exp(-xdamp*alltime)) / xdamp**3
#
#print ("calculated msd: {}".format(msd))
# print (msd2)

fig5, ax = plt.subplots()

nbins = 100
n,f1,patches = ax.hist(xaccdata,nbins,edgecolor='black')

ax.set(xlabel='t (ns)', ylabel='acc (m/s^2)')
fig5.savefig("plots/plot_xacc.png")

#xaccstd = np.std(xaccdata)
#
#xaccstdtrue = np.sqrt(2*xmass*kb*T*xdamp/tstep)
#
#print ("std of acc: {}".format(xaccstd))
#print ("true std of acc: {}".format(xaccstdtrue))

fig6, ax = plt.subplots()
ax.plot(tdata,fdata)
ax.plot(slipst,slips,'or')
ax.set(xlabel='t (nm)', ylabel='x (nm)')
ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
ax.xaxis.grid(True, which='minor',linestyle='dotted')
ax.yaxis.grid(True, linestyle='dotted')
ax.set(xlabel='t (nm)', ylabel='F (nN)')
fig6.savefig("plots/plot_fs.png")

fig7, ax = plt.subplots()
ax.plot(tdata,fdata)
ax.plot(avgst,avgs)
ax.set(xlabel='t (nm)', ylabel='x (nm)')
ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
ax.xaxis.grid(True, which='minor',linestyle='dotted')
ax.yaxis.grid(True, linestyle='dotted')
ax.set(xlabel='t (nm)', ylabel='F (nN)')
fig7.savefig("plots/plot_favg.png")

fig8, ax = plt.subplots()
nbins = 250
n,f1,patches = ax.hist(slipst,nbins,edgecolor='black')
ax.set(xlabel='t (ns)', ylabel='count')
fig8.savefig("plots/slip_hist.png")


