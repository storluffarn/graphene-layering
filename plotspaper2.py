
import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.colors as colors

kb = 1.38064852e-23
xmass = 1e-23
xdamp = 1.875e+13
T = 300

modif = 0.5
tstep = modif * 3e-14	
steps = 1.0/modif * 10.0e5	
alltime = steps*tstep 

traw = np.loadtxt(open("time.csv", "rb"), delimiter=",", skiprows=1)
xraw = np.loadtxt(open("xout.csv", "rb"), delimiter=",", skiprows=1)
qraw = np.loadtxt(open("qout.csv", "rb"), delimiter=",", skiprows=1)
fraw = np.loadtxt(open("tomout.csv", "rb"), delimiter=",", skiprows=1)
sraw = np.loadtxt(open("slips.csv", "rb"), delimiter=",", skiprows=1)
avgraw = np.loadtxt(open("avgs.csv", "rb"), delimiter=",", skiprows=1)
#sangraw = np.loadtxt(open("sangslips.csv", "rb"), delimiter=",", skiprows=0)
#
params = {}
with open("params.dat") as f:
    for line in f:
        (key,val) = line.split()
        params[key] = float(val)

nscale = 1.0e9  # WARNING this is *only* for plotting, danger!

tdata = nscale * traw[:,0]
sdata = nscale * traw[:,1]
xdata = nscale * xraw[:,0]
qdata = nscale * qraw[:,0]
fdata = nscale * fraw[:,2]
supdata = nscale * traw[:,1]
#kindata = fraw[:,0]
avgs  = nscale * avgraw[:,1]
avgst  = nscale * avgraw[:,0]
#sangx = nscale * sangraw[:,0]
#sangy = nscale * sangraw[:,1]
#xaccdata = xraw[:,2]
slipst = nscale * sraw[:,0]
slipsf = nscale * sraw[:,1]
#
#print (sraw)
#print (len(sraw))

## use this for slip pos dist
#if len(sraw) == 0 :     # ugly errpr handling
#    slipst = [0]
#    slipsf = [0]
#    slipspos = [0]
#elif len(sraw) == 3 :   # even uglier, should be it rows = 1
#    slipst = nscale * sraw[0]
#    slipsf = nscale * sraw[1]
#    slipspos = nscale * sraw[2]
#else :
#    slipst = nscale * sraw[:,0]
#    slipsf = nscale * sraw[:,1]
#    slipspos = nscale * sraw[:,2]

# use this for relax time dist dist
#if len(sraw) == 0 :     # ugly errpr handling
#    slipst = [0]
#elif len(sraw) == 1 :   # even uglier, should be it rows = 1
#    slipst = nscale * sraw[0]
#else :
#    slipst = nscale * sraw

fig0, ax = plt.subplots()
ax.plot(tdata,sdata)
fig0.savefig("plots/plot_vt.png")

fig1, ax = plt.subplots()
ax.plot(tdata,xdata)
ax.set(xlabel='t (nm)', ylabel='x (nm)')
#ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
#ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
#ax.xaxis.grid(True, which='minor',linestyle='dotted')
#ax.xaxis.grid(True, which='major',linestyle='dotted')
#ax.yaxis.grid(True, linestyle='dotted')
fig1.savefig("plots/plot_x.png")

fig2, ax = plt.subplots()
ax.plot(tdata,qdata)
ax.set(xlabel='t (nm)', ylabel='q (nm)')
fig2.savefig("plots/plot_q.png")

fig3, ax = plt.subplots()
ax.plot(tdata,fdata)
#ax.plot(supdata,fdata)
ax.set(xlabel='t (nm)', ylabel='x (nm)')
#ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
#ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
#ax.xaxis.grid(True, which='minor',linestyle='dotted')
#ax.xaxis.grid(True, which='major',linestyle='dotted')
#ax.yaxis.grid(True, linestyle='dotted')
ax.set(xlabel='t (nm)', ylabel='F (nN)')
fig3.savefig("plots/plot_f.png")

## for debugging
#
#fig4, ax = plt.subplots()
#
#nbins = 100
#n,f1,patches = ax.hist(slipspos,nbins,edgecolor='black')
#
#fig4.savefig("plots/histogram.png")
#
## mean square displacement
#
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
## for brownian motion
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
#print (msd2)
#
#fig5, ax = plt.subplots()
#
#nbins = 100
#n,f1,patches = ax.hist(xaccdata,nbins,edgecolor='black')
#
#ax.set(xlabel='t (ns)', ylabel='acc (m/s^2)')
#fig5.savefig("plots/plot_xacc.png")
#
#xaccstd = np.std(xaccdata)
#
#xaccstdtrue = np.sqrt(2*xmass*kb*T*xdamp/tstep)
#
#print ("std of acc: {}".format(xaccstd))
#print ("true std of acc: {}".format(xaccstdtrue))
#
## for an harmonic oscilator
#
#kinmean = np.mean(kindata)
#
#print ('mean kinetic energy: {}'.format(kinmean))
#print ('true mean kinetic energy: {}'.format(0.5*kb*T))

fig6, ax = plt.subplots()
ax.plot(tdata,fdata)
ax.plot(slipst,slipsf,'or')
ax.set(xlabel='t (nm)', ylabel='x (nm)')
#ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
#ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
#ax.xaxis.grid(True, which='minor',linestyle='dotted')
#ax.yaxis.grid(True, linestyle='dotted')
ax.set(xlabel='t (nm)', ylabel='F (nN)')
fig6.savefig("plots/plot_fs.png")
#
fig7, ax = plt.subplots()
ax.plot(tdata,fdata)
ax.plot(avgst,avgs)
ax.set(xlabel='t (nm)', ylabel='x (nm)')
#ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
#ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
#ax.xaxis.grid(True, which='minor',linestyle='dotted')
#ax.yaxis.grid(True, linestyle='dotted')
#ax.set(xlabel='t (nm)', ylabel='F (nN)')
fig7.savefig("plots/plot_favg.png")

#fig8, ax = plt.subplots()
#nbins = 250
#n,f1,patches = ax.hist(slipst,nbins,edgecolor='black')
#ax.set(xlabel='t (ns)', ylabel='count')
#fig8.savefig("plots/slip_hist.png")

#fig9, ax = plt.subplots()
#
#xmax = 1.3e-9
#qmax = 0.6e-9
#pixels = 1000
#toevs = 6.242e+18
#
#xlist, ylist = np.mgrid[0:xmax:pixels*1j, 0:qmax:pixels*1j]
#
#spring = params["spring"]
#supvel = params["supvel"]
#latcona = params["latcona"]
#latconb = params["latconb"]
#barr1 = params["barr1"]
#barr2 = params["barr2"]
#kappa1 = params["kappa1"]
#kappa2 = params["kappa2"]
#nu2 = params["nu2"]
#nu4 = params["nu4"]
#
#def pot(x,q) : 
#    out = (barr1 + kappa1*(x-q)**2)*(1-np.cos(2*np.pi*q/latcona)) + (barr2 + kappa2*(x-q)**2)*(1-np.cos(2*np.pi*x/latconb)) 
#    return out
#
#zi = pot(xlist,ylist)
#
#im = ax.pcolormesh(nscale*xlist,nscale*ylist,toevs*zi.reshape(xlist.shape),cmap='YlOrBr',norm=colors.LogNorm(vmin=0.025,vmax=4),rasterized=True)
#ax.plot(xdata,xdata-qdata,'bo',markersize=0.1)
#ax.set(xlabel='x (nm)', ylabel='x-q (nm)')
#
#cbar = fig9.colorbar(im)
#cbar.ax.set_ylabel("$V_\mathrm{tip-sheet} + V_\mathrm{tip-substrate}$")
#
#fig9.savefig("plots/potland.png")
#
#fig10, ax = plt.subplots()
#ax.plot(sangx,sangy)
#fig10.savefig("plots/sang.png")
