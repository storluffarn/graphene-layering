
import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.colors as colors
import time

kb = 1.38064852e-23

modif = 0.5
tstep = modif * 3e-14	
steps = 1.0/modif * 10.0e5	
alltime = steps*tstep 

traw = np.loadtxt(open("time.csv", "rb"), delimiter=",", skiprows=1)
xraw = np.loadtxt(open("xout.csv", "rb"), delimiter=",", skiprows=1)
qraw = np.loadtxt(open("qout.csv", "rb"), delimiter=",", skiprows=1)
fraw = np.loadtxt(open("tomout.csv", "rb"), delimiter=",", skiprows=1)
sraw = np.loadtxt(open("slips.csv", "rb"), delimiter=",", skiprows=0)
storaw = np.loadtxt(open("sliptos.csv", "rb"), delimiter=",", skiprows=0)
sfraw = np.loadtxt(open("slipsf.csv", "rb"), delimiter=",", skiprows=0)
sqraw = np.loadtxt(open("slipsq.csv", "rb"), delimiter=",", skiprows=0)
#sraw = np.loadtxt(open("slipsdetailed.csv", "rb"), delimiter=",", skiprows=1)
avgraw = np.loadtxt(open("avgs.csv", "rb"), delimiter=",", skiprows=1)
#sangraw = np.loadtxt(open("sangslips.csv", "rb"), delimiter=",", skiprows=0)
sangsimraw = np.loadtxt(open("sangcompsim.csv", "rb"), delimiter=",", skiprows=0)
sangtheoraw = np.loadtxt(open("sangcomptheo.csv", "rb"), delimiter=",", skiprows=0)
sangtheohighraw = np.loadtxt(open("sangcomptheohigh.csv", "rb"), delimiter=",", skiprows=0)


params = {}
with open("params.dat") as f:
    for line in f:
        (key,val) = line.split()
        params[key] = float(val)

spring = params["spring"]
supvel = params["supvel"]
latcona = params["latcona"]
latconb = params["latconb"]
barr1 = params["barr1"]
barr2 = params["barr2"]
kappa1 = params["kappa1"]
kappa2 = params["kappa2"]
nu2 = params["nu2"]
nu4 = params["nu4"]

nscale = 1.0e9  # WARNING this is *only* for plotting, danger!

tdata = nscale * traw[:,0]
#sdata = nscale * traw[:,1]
xdata = nscale * xraw[:,0]
qdata = nscale * qraw[:,0]
fdata = nscale * fraw[:,2]
supdata = nscale * traw[:,1]
#kindata = fraw[:,0]
avgs  = nscale * avgraw[:,1]
avgst  = nscale * avgraw[:,0]
qavgs  = nscale * avgraw[:,3]
#sangx = nscale * sangraw[:,0]
#sangy = nscale * sangraw[:,1]
#xaccdata = xraw[:,2]
slipst = nscale * sraw[:,0]
#slipsft = nscale * fraw[:,0]
slipsff = nscale * sraw[:,1]
#slipsqt = nscale * qraw[:,0]
slipsqq = nscale * sraw[:,2]
sliptost = nscale * storaw[:,0]
slipsftof = nscale * storaw[:,1]
slipsqtoq = nscale * storaw[:,2]
sangsimtt = sangsimraw[:,0]
sangsimvv = sangsimraw[:,2]
sangsimtf = nscale * sangsimraw[:,1]
sangsimvf = nscale * sangsimraw[:,3]
sangtheott = sangtheoraw[:,0]
sangtheovv = sangtheoraw[:,2]
sangtheotf = nscale * sangtheoraw[:,1]
sangtheovf = nscale * sangtheoraw[:,3]
sangtheohightt = sangtheohighraw[:,0]
sangtheohighvv = sangtheohighraw[:,2]
sangtheohightf = nscale * sangtheohighraw[:,1]
sangtheohighvf = nscale * sangtheohighraw[:,3]
#tmp = len(simf)
#simf /= simf[round(tmp/2)]
#simv /= simv[round(tmp/2)]
#sangf /= sangf[round(tmp/2)]
#sangv /= sangv[round(tmp/2)]

#tdata = tdata[:round(len(tdata)*0.5)]
#fdata = fdata[:round(len(fdata)*0.5)]
#qdata = qdata[:round(len(qdata)*0.5)]
#avgs = avgs[:round(len(avgs)*0.5)]
#avgst = avgst[:round(len(avgst)*0.5)]
#qavgs = qavgs[:round(len(qavgs)*0.5)]

#
##print (sraw)
##print (len(sraw))
#
### use this for slip pos dist
#if len(qraw) == 0 :     # ugly errpr handling
#    slipsqt = [0]
#    slipsqq = [0]
#elif len(qraw) == 3 :   # even uglier, should be it rows = 1
#    slipsqt = nscale * qraw[0]
#    slipsqq = nscale * qraw[1]
#else :
#    slipsqt = nscale * qraw[:,0]
#    slipsqq = nscale * qraw[:,1]
#
## use this for relax time dist dist
##if len(sraw) == 0 :     # ugly errpr handling
##    slipst = [0]
##elif len(sraw) == 1 :   # even uglier, should be it rows = 1
##    slipst = nscale * sraw[0]
##else :
##    slipst = nscale * sraw
#
### support movement
#
#fig0, ax = plt.subplots()
#ax.plot(tdata,sdata)
#fig0.savefig("plots/plot_vt.png")
#

### standard plots

#fig1, ax = plt.subplots()
#ax.plot(tdata,xdata)
#ax.set(xlabel='t (nm)', ylabel='x (nm)')
#ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
#ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
#ax.xaxis.grid(True, which='minor',linestyle='dotted')
#ax.xaxis.grid(True, which='major',linestyle='dotted')
#ax.yaxis.grid(True, linestyle='dotted')
#fig1.savefig("plots/plot_x.png")
#
#fig2, ax = plt.subplots()
#ax.plot(tdata,qdata)
##ax.plot(supdata,qdata)
#ax.set(xlabel='t (nm)', ylabel='q (nm)')
#fig2.savefig("plots/plot_q.png")
#
#fig3, ax = plt.subplots()
#ax.plot(tdata,fdata)
##ax.plot(supdata,fdata)
#ax.set(xlabel='t (nm)', ylabel='x (nm)')
#ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
#ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
#ax.xaxis.grid(True, which='minor',linestyle='dotted')
#ax.xaxis.grid(True, which='major',linestyle='dotted')
#ax.yaxis.grid(True, linestyle='dotted')
#ax.set(xlabel='t (nm)', ylabel='F (nN)')
#fig3.savefig("plots/plot_f.png")
#
### for debugging
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
#
### Slipping points in noysy signal
#
#fig6, ax = plt.subplots()
#ax.plot(tdata,fdata)
#ax.plot(slipsft,slipsff,'or')
##ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
##ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
##ax.xaxis.grid(True, which='minor',linestyle='dotted')
##ax.yaxis.grid(True, linestyle='dotted')
#ax.set(xlabel='t (nm)', ylabel='F (nN)')
#fig6.savefig("plots/plot_fs.png")
##
#fig6b, ax = plt.subplots()
#ax.plot(tdata,qdata)
#ax.plot(slipsqt,slipsqq,'or')
##ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
##ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
##ax.xaxis.grid(True, which='minor',linestyle='dotted')
##ax.yaxis.grid(True, linestyle='dotted')
#ax.set(xlabel='t (nm)', ylabel='q (nm)')
#fig6b.savefig("plots/plot_qs.png")
##
#fig7, ax = plt.subplots()
#f1, = ax.plot(tdata,fdata, label = "Signal")
#f2, = ax.plot(avgst,avgs, "Noise reduction")
#ax.set(xlabel='t (nm)', ylabel='x (nm)')
##ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
##ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
##ax.xaxis.grid(True, which='minor',linestyle='dotted')
##ax.yaxis.grid(True, linestyle='dotted')
##ax.set(xlabel='t (nm)', ylabel='F (nN)')
#fig7.savefig("plots/plot_favg.png")
#
#fig7b, ax = plt.subplots()
#f1, = ax.plot(tdata,qdata, label = "Signal")
#f2, = ax.plot(avgst,qavgs, label = "Noise reduction")
##ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
##ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
##ax.xaxis.grid(True, which='minor',linestyle='dotted')
##ax.yaxis.grid(True, linestyle='dotted')
#ax.set(xlabel='t (nm)', ylabel='q (nN)')
#fig7b.savefig("plots/plot_qavg.png")

### Historgram

#fig8, ax = plt.subplots()
#nbins = 250
#n,f1,patches = ax.hist(slipst,nbins,edgecolor='black')
#ax.set(xlabel='t (ns)', ylabel='count')
#fig8.savefig("plots/slip_hist.png")

### paths in potential landscape

#fig9, ax = plt.subplots()
#
#xmax = 1.3e-9
#qmax = 0.6e-9
#pixels = 1000
#toevs = 6.242e+18
#
#xlist, ylist = np.mgrid[0:xmax:pixels*1j, 0:qmax:pixels*1j]
#
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
#
## f-q with slips

fig11, (ax1,ax2) = plt.subplots(nrows = 2, sharex = True)

f12, = ax1.plot(tdata,qdata, label = "Signal")
f22, = ax1.plot(avgst,qavgs, label = "Noise reduction")
ax1.set(ylabel='q (nm)')
ax1.xaxis.set_major_locator(ticker.MultipleLocator(5.0))
ax1.xaxis.set_minor_locator(ticker.MultipleLocator(2.5))
#ax1.xaxis.grid(True, which='minor', linestyle='dotted')
#ax1.xaxis.grid(True, which='major', linestyle='dotted')
#ax1.yaxis.grid(True, which='major', linestyle='dotted')
ax1.legend(loc='upper right')

for slip in slipst:
    ax1.axvline(x=slip, color = 'gray', linestyle = 'dotted')

f21, = ax2.plot(tdata,fdata, label = "Signal")
f22, = ax2.plot(avgst,avgs, label = "Noise reduction")
ax2.set(xlabel='t (ns)', ylabel='F (nN)')
ax2.xaxis.set_major_locator(ticker.MultipleLocator(5.0))
ax2.xaxis.set_minor_locator(ticker.MultipleLocator(1.0))
#ax2.xaxis.grid(True, which='minor', linestyle='dotted')
#ax2.xaxis.grid(True, which='major', linestyle='dotted')
#ax2.yaxis.grid(True, which='major', linestyle='dotted')
#ax2.legend(loc='upper right')

for slip in slipst :
    ax2.axvline(x=slip, color = 'gray', linestyle = 'dotted')
#for slipto in sliptost:
#    ax2.axvline(x=slipto, color = 'gray', linestyle = 'dotted')

plt.subplots_adjust(hspace=0)

timestamp = time.strftime("%Y%m%d-%H%M%S")

fig11.savefig("plots/fq_slips")
#fig11.savefig("plots/fq_slips{}".format(timestamp))

## stable unstable 

#fig12, (ax2,ax1) = plt.subplots(nrows = 2, sharex = True)
#
#ax1.plot(tdata,fdata, color='#ff7f0e', label = "")
#ax1.set(ylabel='F (nN)')
#ax1.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
#ax1.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
#ax1.xaxis.grid(True, which='minor', linestyle='dotted')
#ax1.xaxis.grid(True, which='major', linestyle='dotted')
#ax1.yaxis.grid(True, which='major', linestyle='dotted')
#ax1.set(xlabel='t (nm)', ylabel='q (nN)')
##ax1.legend(loc='lower right')
#
#f2, = ax2.plot(tdata,qdata, color='#ff7f0e', label = "")
#ax2.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
#ax2.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
#ax2.xaxis.grid(True, which='minor', linestyle='dotted')
#ax2.xaxis.grid(True, which='major', linestyle='dotted')
#ax2.yaxis.grid(True, which='major', linestyle='dotted')
##ax2.legend(loc='lower right')
#
#plt.subplots_adjust(hspace=0)
#
#fig12.savefig("plots/fqthermal.png")
##
#fig13, ax = plt.subplots()
#ax.plot(sangsimtt,sangsimtf,'*',label='simulation')
#ax.plot(sangtheohightt,sangtheohightf,linestyle='dashed',color='#ff7f0e', label='theory high')
#ax.plot(sangtheott,sangtheotf, label='theory low')
#ax.set(xlabel='T (K)', ylabel='F (nN)')
##ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
##ax.xaxis.set_minor_locator(icker.MultipleLocator(0.25))
##ax.xaxis.grid(True, which='minor',linestyle='dotted')
##ax.xaxis.grid(True, which='major',linestyle='dotted')
##ax.yaxis.grid(True, linestyle='dotted#i')
#ax.legend(loc='lower right')
#ax.set_ylim([0,4])
#
#fig13.savefig("plots/sangcompf.png")
#
#fig14, ax = plt.subplots()
#ax.plot(sangsimvv,sangsimvf,'*',label='simulation')
#ax.plot(sangtheohighvv,sangtheohighvf,color='#ff7f0e',linestyle='dashed',label='theory high')
#ax.plot(sangtheovv,sangtheovf,color='#ff7f0e', label='theory low')
#
#ax.set(xlabel='v (m/s)', ylabel='F (nN)')
##ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
##ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
##ax.xaxis.grid(True, which='minor',linestyle='dotted')
##ax.xaxis.grid(True, which='major',linestyle='dotted')
##ax.yaxis.grid(True, linestyle='dotted')
#ax.legend(loc='lower right')
#ax.set_ylim([0,4])
#
#fig14.savefig("plots/sangcompv.png")
#
#
### Landscape plots (remastered from Mathematica)

### density plots

#fig15, ax = plt.subplots()
#
## constants
#
#xmax = 6.5e-9
#qmax = 1e-9
#bins = 1000
#
##xlist, qlist = np.mgrid[0:xmax:bins*1j, 0:qmax:bins*1j]
#xlist, qlist = np.mgrid[-0.25e-9:xmax:bins*1j, -qmax:qmax:bins*1j] ## use this for retrace
#
#nuscale = 0.0;
#
#def pot(x,q) :
#    return (barr1 + kappa1*q**2) * (1 - np.cos(2*np.pi/latcona*(x - q))) + (barr2 + kappa2*q**2) * (1 - np.cos(2*np.pi/latconb*x)) + nuscale * (nu2*q**2 + nu4*q**4)
#
## + n2*np.square(qlist) + n4*np.square(qlist) # not a proper part of the land scape
#
#zi = pot(xlist,qlist)
#
#toev = 6.242e+18
#
#im = ax.pcolormesh(nscale*xlist,nscale*qlist,toev*zi.reshape(xlist.shape),cmap='YlOrBr',rasterized=True)
#ax.plot(xdata,qdata)
#ax.set(xlabel='$x$ (nm)', ylabel='$q$ (nm)')
#
#cbar = fig15.colorbar(im)
#cbar.set_label(label='$V_\\mathrm{tip-sheet} + V_\\mathrm{tip-substrate}$')
#
#fig15.savefig("plots/landscape.png")
#

### retrace plots

#fig16, ax = plt.subplots()
#ax.plot(supdata,fdata)
#ax.set(xlabel='support position (nm)', ylabel='F (nN)')
#ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
#ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
#ax.xaxis.grid(True, which='minor',linestyle='dotted')
#ax.xaxis.grid(True, which='major',linestyle='dotted')
#ax.yaxis.grid(True, linestyle='dotted')
#fig16.savefig("plots/retrace_f.png")
#
#fig17, ax = plt.subplots()
#ax.plot(supdata,xdata)
#ax.set(xlabel='support position (nm)', ylabel='x (nm)')
#ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
#ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
#ax.xaxis.grid(True, which='minor',linestyle='dotted')
#ax.xaxis.grid(True, which='major',linestyle='dotted')
#ax.yaxis.grid(True, linestyle='dotted')
#fig17.savefig("plots/retrace_x.png")
#
#fig18, ax = plt.subplots()
#ax.plot(supdata,qdata)
#ax.set(xlabel='support position (nm)', ylabel='q (nm)')
#ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
#ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
#ax.xaxis.grid(True, which='minor',linestyle='dotted')
#ax.xaxis.grid(True, which='major',linestyle='dotted')
#ax.yaxis.grid(True, linestyle='dotted')
#fig18.savefig("plots/retrace_q.png")



