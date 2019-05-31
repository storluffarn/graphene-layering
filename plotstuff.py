
# simple plotting test program

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker

font = {'size':12}
matplotlib.rc('font', **font)

# import csv
xraw = np.loadtxt(open("stdx.csv", "rb"), delimiter=",", skiprows=0)
qraw = np.loadtxt(open("stdq.csv", "rb"), delimiter=",", skiprows=0)
fraw = np.loadtxt(open("stdfric.csv", "rb"), delimiter=",", skiprows=0)
icraw = np.loadtxt(open("incomf.csv", "rb"), delimiter=",", skiprows=0)
v1raw = np.loadtxt(open("v1f.csv", "rb"), delimiter=",", skiprows=0)
v2raw = np.loadtxt(open("v2f.csv", "rb"), delimiter=",", skiprows=0)
icrawx = np.loadtxt(open("incomx.csv", "rb"), delimiter=",", skiprows=0)
v1rawx = np.loadtxt(open("v1x.csv", "rb"), delimiter=",", skiprows=0)
v2rawx = np.loadtxt(open("v2x.csv", "rb"), delimiter=",", skiprows=0)
v21rawx = np.loadtxt(open("v21x.csv", "rb"), delimiter=",", skiprows=0)
v22rawx = np.loadtxt(open("v22x.csv", "rb"), delimiter=",", skiprows=0)
icrawq = np.loadtxt(open("incomq.csv", "rb"), delimiter=",", skiprows=0)
v1rawq = np.loadtxt(open("v1q.csv", "rb"), delimiter=",", skiprows=0)
v2rawq = np.loadtxt(open("v2q.csv", "rb"), delimiter=",", skiprows=0)
g1rawx = np.loadtxt(open("g1.1xr.csv", "rb"), delimiter=",", skiprows=0)
g1rawq = np.loadtxt(open("g1.1qr.csv", "rb"), delimiter=",", skiprows=0)
g1rawf = np.loadtxt(open("g1.1fr.csv", "rb"), delimiter=",", skiprows=0)
g2rawt = np.loadtxt(open("g2.1xr.csv", "rb"), delimiter=",", skiprows=0)
g2rawq = np.loadtxt(open("g2.1qr.csv", "rb"), delimiter=",", skiprows=0)
g2rawf = np.loadtxt(open("g2.1fr.csv", "rb"), delimiter=",", skiprows=0)
v21rawq = np.loadtxt(open("v21q.csv", "rb"), delimiter=",", skiprows=0)
v22rawq = np.loadtxt(open("v22q.csv", "rb"), delimiter=",", skiprows=0)
v21rawf = np.loadtxt(open("v21f.csv", "rb"), delimiter=",", skiprows=0)
v22rawf = np.loadtxt(open("v22f.csv", "rb"), delimiter=",", skiprows=0)
irrrawx = np.loadtxt(open("irrx.csv", "rb"), delimiter=",", skiprows=0)
irrrawq = np.loadtxt(open("irrq.csv", "rb"), delimiter=",", skiprows=0)
irrrawf = np.loadtxt(open("irrf.csv", "rb"), delimiter=",", skiprows=0)
fnu1raw = np.loadtxt(open("nunu.dat", "rb"), delimiter=",", skiprows=0)
fnu2raw = np.loadtxt(open("nufric.dat", "rb"), delimiter=",", skiprows=0)
fnu3raw = np.loadtxt(open("fricnuan.csv", "rb"), delimiter=",", skiprows=0)
thrmrawx = np.loadtxt(open("thrx.csv", "rb"), delimiter=",", skiprows=0)
thrmrawq = np.loadtxt(open("thrq.csv", "rb"), delimiter=",", skiprows=0)
thrmrawf = np.loadtxt(open("thrf.csv", "rb"), delimiter=",", skiprows=0)

# slice list and scale it to nm
xscale = 1e9
qscale = 1e9
tscale = 1e9
fscale = 1e9
nscale = 1.0/(7**-4*14e20)

xdata = xscale * xraw[:,1]
qdata = qscale * qraw[:,1]
fdata = fscale * fraw[:,1]
tdata = tscale * xraw[:,0] # multiple data have the have the same time steps
icdata = fscale * icraw[:,1]
v1data = fscale * v1raw[:,1]
v2data = fscale * v2raw[:,1]
icxdata = fscale * icrawx[:,1]
v1xdata = fscale * v1rawx[:,1]
v2xdata = fscale * v2rawx[:,1]
icqdata = fscale * icrawq[:,1]
v1qdata = fscale * v1rawq[:,1]
v2qdata = fscale * v2rawq[:,1]
g1xdata = fscale * g1rawx[:,1]
g1qdata = fscale * g1rawq[:,1]
g1fdata = fscale * g1rawf[:,1]
g2tdata = fscale * g2rawt[:,0]
g2qdata = fscale * g2rawq[:,1]
g2fdata = fscale * g2rawf[:,1]
v21xdata = fscale * v21rawx[:,1]
v21qdata = fscale * v21rawq[:,1]
v21fdata = fscale * v21rawf[:,1]
v22xdata = fscale * v22rawx[:,1]
v22qdata = fscale * v22rawq[:,1]
v22fdata = fscale * v22rawf[:,1]
irrxdata = fscale * irrrawx[:,1]
irrqdata = fscale * irrrawq[:,1]
irrfdata = fscale * irrrawf[:,1]
tdata2 = tscale * v1raw[:,0] # v plots have longer time
tdata3 = tscale * v2raw[:,0] # v2 has half the number of time steps
nunudata = nscale * fnu1raw
nufrdata = fscale * fnu2raw
annudata = nscale * fnu3raw[:,0]
anfrdata = fscale * fnu3raw[:,1]
thrmxdata = xscale * thrmrawx[:,1]
thrmqdata = qscale * thrmrawq[:,1]
thrmfdata = fscale * thrmrawf[:,1]

print("imports done!")

### simple two stacked figures

fig1, (ax1, ax2) = plt.subplots(2,1,sharex=True)

# plotting stuff
ax1.plot(tdata,xdata,color='#ff7f0e')
ax1.set(ylabel='x (nm)')
ax1.grid(linestyle='dotted')

ax2.plot(tdata,qdata,color='#2ca02e')
ax2.set(xlabel='d (nm)', ylabel='x (nm)')
ax2.grid(linestyle='dotted')

# saving stuff
fig1.savefig("test.png")

print("figure 1 done!")

### two stacked figures, no spacing and lines for lattice periods

fig2, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)

# because the data is in nm
f1, = ax1.plot(tdata, fdata, color='#ff7f0e')
f2, = ax2.plot(tdata, qdata, color='#ff7f0e')

ax1.set(ylabel='F (nN)')
ax1.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
ax1.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
ax1.xaxis.grid(True, which='minor', linestyle='dotted')
ax1.yaxis.grid(True, linestyle='dotted')

#ax1.annotate('A', (0, 1), xytext=(4, -4), xycoords='axes fraction', textcoords='offset points', fontweight='bold', backgroundcolor='w', ha='left', va='top')

ax2.tick_params(top = True)
ax2.set(xlabel='vt (nm)', ylabel='q (nm)')
ax2.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
ax2.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
ax2.xaxis.grid(True, which='minor', linestyle='dotted')
ax2.yaxis.grid(True, linestyle='dotted')
#ax2.legend((f1, f2), ('x(t)', 'q(t)'), loc='lower right')

yticks = ax2.yaxis.get_major_ticks()
yticks[-1].set_visible(False)
plt.subplots_adjust(hspace=0)

# saving stuff
fig2.savefig("test2.png")
print("figure 2 done!")

### density plots

fig3, ax = plt.subplots()

# constants

oos = 1.0/7.0 #ugly constant

un = 4.0e-20
kn = oos**2 * 3.0
en = 15e12
k  = 2.0
v  = 1.0
u1 = un
u2 = 0.5*un
k1 = kn
k2 = 0.5*kn
b  = 1.0
g  = 1.0
n2 = oos**2 * 18.75
n4 = oos**4 * 14e-20
a  = 2.5e-10
mx = 1e-23
mq = oos * 0.25e-22
ex = 1.25*en
eq = oos * 20*en
h  = 1e-15
tmax = 5e6

xmax = 4e-9
qmax = 1e-9
bins = 1000

xlist, qlist = np.mgrid[0:xmax:bins*1j, 0:qmax:bins*1j]

def pot(x,q) :
    tmp = (u1 + k1*q**2) * (1 - np.cos(2*np.pi/a*(x - b*q))) + (u2 + k2*qlist**2) * (1 - np.cos(g*2*np.pi/a*x))
    return tmp
# + n2*np.square(qlist) + n4*np.square(qlist) # not a proper part of the land scape

zi = pot(xlist,qlist)

toev = 6.242e+18

im = ax.pcolormesh(xscale*xlist,qscale*qlist,toev*zi.reshape(xlist.shape),cmap='YlOrBr',rasterized=True)
ax.plot(xdata,qdata)
ax.set(xlabel='x (nm)', ylabel='q (nm)')

cbar = fig3.colorbar(im)
cbar.ax.set_ylabel("$V_\mathrm{tip-sheet} + V_\mathrm{tip-substrate}$")

fig3.savefig("test3.png")
print("figure 3 done!")

### tryin go get it all into one plot

nx = 2
ny = 2

fig4, f4_axs = plt.subplots(ncols=nx, nrows=ny, figsize=(7,3), sharex=True)

plt.subplots_adjust(hspace=0)

gs = f4_axs[0, 1].get_gridspec()

# remove the underlying axes
#for i in range (ny) : 
#    for ax in f7_axs[0:, i] :
#        ax.remove()

f4_axs[0,1].remove()
f4_axs[1,1].remove()

#f7_axs[0,0].remove()

axbig = fig4.add_subplot(gs[0:, -1])

fig4.tight_layout() # need to comment out to make hspace = 0 to work

fig4.savefig("test4.png")
print("figure 4 done!")

### Analytical plot

fig5, ax = plt.subplots()

# plotting stuff

# experimental data

# sims
dong = np.array([1.5, 0.9, 0.6, 0.5]) / 1.5 # graphene
li = np.array([8.0, 6.0, 4.8, 4.6]) / 8.0 # graphene


#exps from "Frictional Characteristics of Atomically Thin Sheets"
carp1 = [1.0, 0.80, 0.55, 0.45] # graphene
carp2 = [1.0, 0.82, 0.75] 
carp3 = [1.0, 0.78] 
carp4 = [1.0, 0.5] # MoS2
carp5 = [1.0, 0.9, 0.8, 0.55]
carp6 = [1.0, 0.75]
carp7 = [1.0, 0.83, 0.67, 0.60] # NbSe2

x2 = [1,2]
x22 = [1,3]
x3 = [1,2,3]
x32 = [1,3,4]
x4 = [1,2,3,4]

# normalizing

anscale = fscale * 2.06121e-9 # this is f(nu_4)

nufrdata = nufrdata / anscale
anfrdata = anfrdata / anscale

f1, = ax.plot(nunudata,nufrdata,'g^')
f2, = ax.plot(annudata,anfrdata,color='#ff7f0e')

f3, = ax.plot(x4,carp1,'bo') # graphene
ax.plot(x3,carp2,'bo')
ax.plot(x2,carp3,'bo')
f4, = ax.plot(x22,carp4,'bs') # MoS2
ax.plot(x4,carp5,'bs')
ax.plot(x22,carp6,'bs')
f5, = ax.plot(x4,carp7,'bD') # NbSe2

f6, = ax.plot(x4,dong,'mP')
f7, = ax.plot(x4,li,'mX')

ax.set_xlim(0.0,4.5)
ax.set_ylim(0.0,1.75)
ax.set(xlabel='$\\nu_4 / \\nu_4^0$ or \n number of layers', ylabel='F/F($\\nu_4^0$) or \n friction relative to one layer' )
ax.legend((f1, f2, f3, f4, f5, f6, f7), ('model sim.', 'analytic', 'graphene [22]', 'MoS2 [22]', 'NbSe2 [22]', 'graphene [6]', 'grahpene [24]'), loc='upper right')

ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.5))
ax.xaxis.grid(linestyle='dotted',which='minor')
ax.yaxis.grid(linestyle='dotted')

# saving stuff
fig5.tight_layout()
fig5.savefig("test5.png")
print("figure 5 done!")

### incom plot

fig6, ax = plt.subplots()

# plotting stuff
ax.plot(tdata,icdata,color='#ff7f0e')
ax.set(xlabel='d (nm)', ylabel='F (nN)')
ax.grid(linestyle='dotted')

## in fig2 (std fric) style

fig6a, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)

# because the data is in nm
f1, = ax1.plot(tdata, icdata, color='#ff7f0e')
f2, = ax2.plot(tdata, icqdata, color='#ff7f0e')

ax1.set(ylabel='F (nN)')
ax1.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
ax1.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
ax1.xaxis.grid(True, which='minor',linestyle='dotted')
ax1.yaxis.grid(True, linestyle='dotted')

ax2.tick_params(top = True)
ax2.set(xlabel='vt (nm)', ylabel='q (nm)')
ax2.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
ax2.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
ax2.xaxis.grid(True, which='minor', linestyle='dotted')
ax2.yaxis.grid(True, linestyle='dotted')
#ax2.legend((f21, f22), ('v = 1 m/s', 'v = 2 m/s'), loc='lower right')

yticks = ax2.yaxis.get_major_ticks()
yticks[-1].set_visible(False)
plt.subplots_adjust(hspace=0)

## potland

fig6b, ax = plt.subplots()

g = (1 + np.sqrt(5))/2.0

zi = pot(xlist,qlist)

toev = 6.242e+18

im = ax.pcolormesh(xscale*xlist,qscale*qlist,toev*zi.reshape(xlist.shape),cmap='YlOrBr',rasterized=True)
ax.plot(icxdata,icqdata)
ax.set(xlabel='x (nm)', ylabel='q (nm)')

cbar = fig6b.colorbar(im)
cbar.ax.set_ylabel("$V_\mathrm{tip-sheet} + V_\mathrm{tip-substrate}$")

## moire stuff

fig6a2, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)

# because the data is in nm
f10, = ax1.plot(g2tdata, g2fdata, color='#2ca02e')
f11, = ax1.plot(g2tdata, g1fdata, color='#ff7f0e')
f20, = ax2.plot(g2tdata, g2qdata, color='#2ca02e')
f21, = ax2.plot(g2tdata, g1qdata, color='#ff7f0e')

ax1.set(ylabel='F (nN)')
ax1.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
#ax1.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
ax1.xaxis.grid(True, which='major',linestyle='dotted')
ax1.yaxis.grid(True, linestyle='dotted')

ax2.tick_params(top = True)
ax2.set(xlabel='vt (nm)', ylabel='q (nm)')
ax2.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
#ax2.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
ax2.xaxis.grid(True, which='major', linestyle='dotted')
ax2.yaxis.grid(True, linestyle='dotted')
ax2.legend((f21, f20), ('deformable sheet','rigid sheet'), loc='lower right')

yticks = ax2.yaxis.get_major_ticks()
yticks[-1].set_visible(False)
plt.subplots_adjust(hspace=0)


fig6b2, ax = plt.subplots()

xmax2 = 12.5e-9
qmax2 = 0.9e-9
bins2 = 1000

xlist2, qlist2 = np.mgrid[0:xmax2:bins2*1j, 0:qmax2:bins2*1j]

g = 0.95

zi = pot(xlist2,qlist2)

toev = 6.242e+18

im = ax.pcolormesh(xscale*xlist2,qscale*qlist2,toev*zi.reshape(xlist2.shape),cmap='YlOrBr',rasterized=True)
ax.plot(g1xdata,g1qdata)
ax.set(xlabel='x (nm)', ylabel='q (nm)')

cbar = fig6b2.colorbar(im)
cbar.ax.set_ylabel("$V_\mathrm{tip-sheet} + V_\mathrm{tip-substrate}$")


# saving stuff
fig6.savefig("test6.png")
fig6a.savefig("test6a.png")
fig6b.savefig("test6b.png")
fig6a2.savefig("test6a2.png")
fig6b2.savefig("test6b2.png")

print("figure 6 done!")

### velocity dep

fig7, ax = plt.subplots()

# plotting stuff

f11, = ax.plot(tdata2, v1data, color='#ff7f0e')
f12, = ax.plot(tdata3, v2data, color='#2ca02e')
f21, = ax.plot(tdata2, v21fdata, color='#ff7f0e', dashes=[3,2])
f22, = ax.plot(tdata3, v22fdata, color='#2ca02e', dashes=[3,2])

ax.set(xlabel='vt (nm)', ylabel='F (nN)' )
ax.legend((f21, f22, f11, f12), ('v = 1 m/s, std params', 'v = 2 m/s, std params', 'v = 1 m/s, mod params', 'v = 2 m/s, mod params'), loc='center right')

ax.grid(linestyle='dotted') 

## in fig2 (std fric) style

fig7a, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)

# because the data is in nm
f11, = ax1.plot(tdata2, v1data, color='#ff7f0e')
f12, = ax1.plot(tdata3, v2data, color='#2ca02e')
f21, = ax2.plot(tdata2, v1qdata, color='#ff7f0e')
f22, = ax2.plot(tdata3, v2qdata, color='#2ca02e')

ax1.set(ylabel='F (nN)')
ax1.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
ax1.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
#ax1.xaxis.grid(True, which='minor')
ax1.grid(True, linestyle='dotted')

ax2.tick_params(top = True)
ax2.set(xlabel='vt (nm)', ylabel='q (nm)')
ax2.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
ax2.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
ax2.grid(True, linestyle='dotted')
ax2.legend((f21, f22), ('v = 1 m/s', 'v = 2 m/s'), loc='lower right')

yticks = ax2.yaxis.get_major_ticks()
yticks[-1].set_visible(False)
plt.subplots_adjust(hspace=0)

## potland

fig7b, ax = plt.subplots()

#un = 1.3e-20
#kn = 0.2
#u1 = un 
#k1 = kn
#u2 = 0.0*un
#k2 = 0.0*kn 
#b = 1.0 
#g = 1.0
#n2 = 0*18.75
#n4 = 8e18 
#a = 2.5e-10

xmax = 6.5e-9
qmax = 1.1e-9
bins = 1000

xlist, qlist = np.mgrid[0:xmax:bins*1j, 0:qmax:bins*1j]

def pot(x,q) :
    tmp = (u1 + k1*q**2) * (1 - np.cos(2*np.pi/a*(x - b*q))) 
    return tmp
# + n2*np.square(qlist) + n4*np.square(qlist) # not a proper part of the land scape

zi = pot(xlist,qlist)

toev = 6.242e+18

im = ax.pcolormesh(xscale*xlist,qscale*qlist,toev*zi.reshape(xlist.shape),cmap='YlOrBr',rasterized=True)
f1, = ax.plot(v1xdata[:-7000],v1qdata[:-7000])
f2, = ax.plot(v2xdata[:-3500],v2qdata[:-3500],color='purple')
f3, = ax.plot(v21xdata,v21qdata,dashes=[3,2],color='#1f77b4')
f4, = ax.plot(v22xdata,v22qdata,dashes=[3,2],color='purple')
ax.legend((f3, f4, f1, f2), ('v = 1 m/s, std params', 'v = 2 m/s, std params', 'v = 1 m/s, mod params', 'v = 2 m/s, mod params'), loc='center right')
ax.set(xlabel='x (nm)', ylabel='q (nm)')
cbar = fig7b.colorbar(im)
cbar.ax.set_ylabel("$V_\mathrm{tip-sheet}$")

## only mod pram

fig7b2, ax = plt.subplots()

xmax = 6.5e-9
qmax = 0.3e-9
bins = 1000

xlist, qlist = np.mgrid[0:xmax:bins*1j, 0:qmax:bins*1j]

def pot(x,q) :
    tmp = (u1 + k1*q**2) * (1 - np.cos(2*np.pi/a*(x - b*q))) 
    return tmp
# + n2*np.square(qlist) + n4*np.square(qlist) # not a proper part of the land scape

zi = pot(xlist,qlist)

toev = 6.242e+18

im = ax.pcolormesh(xscale*xlist,qscale*qlist,toev*zi.reshape(xlist.shape),cmap='YlOrBr',rasterized=True)
f1, = ax.plot(v1xdata[:-7000],v1qdata[:-7000])
f2, = ax.plot(v2xdata[:-3500],v2qdata[:-3500],color='purple')
ax.legend((f1, f2), ('v = 1 m/s', 'v = 2 m/s'), loc='lower right')
ax.set(xlabel='x (nm)', ylabel='q (nm)')
cbar = fig7b2.colorbar(im)
cbar.ax.set_ylabel("$V_\mathrm{tip-sheet}$")

## without new parameters

fig7c, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)

# because the data is in nm
f11, = ax1.plot(tdata2, v1data, color='#ff7f0e')
f12, = ax1.plot(tdata3, v2data, color='#2ca02e')
f21, = ax2.plot(tdata2, v21fdata, color='#ff7f0e')
f22, = ax2.plot(tdata3, v22fdata, color='#2ca02e')

ax1.set(ylabel='F (nN)')
ax1.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
ax1.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
#ax1.xaxis.grid(True, which='minor')
ax1.grid(True, linestyle='dotted')

ax2.tick_params(top = True)
ax2.set(xlabel='vt (nm)', ylabel='F (nN)')
ax2.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
ax2.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
ax2.grid(True, linestyle='dotted')
ax2.legend((f21, f22), ('v = 1 m/s', 'v = 2 m/s'), loc='lower right')

yticks = ax2.yaxis.get_major_ticks()
yticks[-1].set_visible(False)
plt.subplots_adjust(hspace=0)

# saving stuff
fig7.savefig("test7.png")
fig7a.savefig("test7a.png")
fig7b.savefig("test7b.png")
fig7c.savefig("test7c.png")
fig7b2.savefig("test7b2.png")
print("figure 7 done!")

### Sketch

fig8, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)

xmax = 12*np.pi 
qmax = 6*np.pi 
bins = 1000

xlist, qlist = np.mgrid[0:xmax:bins*1j, 0:qmax:bins*1j]

def pot1(x,q) :
    tmp = np.cos(x-1.0*q) - 0.5*np.cos(x)
    return tmp

def pot2(x,q) :
    tmp = np.cos(x-1.0*q) - 0.5*np.cos((1+np.sqrt(5)/2.0)*x)
    return tmp

z1 = pot1(xlist,qlist)
z2 = pot2(xlist,qlist)

# plasma or inferno?
ax1.contour(xlist,qlist,z1.reshape(xlist.shape),colors='black')
ax1.contourf(xlist,qlist,z1.reshape(xlist.shape),cmap='plasma',alpha=0.75)
ax2.pcolormesh(xlist,qlist,z2.reshape(xlist.shape),cmap='inferno',rasterized=True)

#ax1.set_aspect(0.5)
#ax2.set_aspect(0.5)

ax1.tick_params(labelleft = False)
ax2.tick_params(labelleft = False)
ax2.tick_params(labelbottom = False)

# ax.set(xlabel='x (nm)', ylabel='q (nm)')

fig8.savefig("test8.png")
print("figure 8 done!")

### irregular substrate

fig9, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)

# because the data is in nm
f1, = ax1.plot(tdata, irrfdata, color='#ff7f0e')
f2, = ax2.plot(tdata, irrqdata, color='#ff7f0e')

ax1.set(ylabel='F (nN)')
ax1.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
ax1.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
ax1.xaxis.grid(True, which='minor', linestyle='dotted')
ax1.yaxis.grid(True, which='major', linestyle='dotted')

#ax1.annotate('A', (0, 1), xytext=(4, -4), xycoords='axes fraction', textcoords='offset points', fontweight='bold', backgroundcolor='w', ha='left', va='top')

ax2.tick_params(top = True)
ax2.set(xlabel='vt (nm)', ylabel='q (nm)')
ax2.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
ax2.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
ax2.xaxis.grid(True, which='minor', linestyle='dotted')
ax2.yaxis.grid(True, which='major', linestyle='dotted')
#ax2.legend((f1, f2), ('x(t)', 'q(t)'), loc='lower right')

yticks = ax2.yaxis.get_major_ticks()
yticks[-1].set_visible(False)
plt.subplots_adjust(hspace=0)

# saving stuff
fig9.savefig("test9.png")
print("figure 9 done!")

### 

fig10, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)

# because the data is in nm
f1, = ax1.plot(tdata, thrmfdata, color='#ff7f0e', marker = 'o', markersize=1)
f2, = ax2.plot(tdata, thrmqdata, color='#ff7f0e', marker = 'o', markersize=1)

ax1.set(ylabel='F (nN)')
ax1.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
ax1.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
ax1.xaxis.grid(True, which='minor', linestyle='dotted')
ax1.yaxis.grid(True, which='major', linestyle='dotted')

#ax1.annotate('A', (0, 1), xytext=(4, -4), xycoords='axes fraction', textcoords='offset points', fontweight='bold', backgroundcolor='w', ha='left', va='top')

ax2.tick_params(top = True)
ax2.set(xlabel='vt (nm)', ylabel='q (nm)')
ax2.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
ax2.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
ax2.xaxis.grid(True, which='minor', linestyle='dotted')
ax2.yaxis.grid(True, which='major', linestyle='dotted')
#ax2.legend((f1, f2), ('x(t)', 'q(t)'), loc='lower right')

yticks = ax2.yaxis.get_major_ticks()
yticks[-1].set_visible(False)
plt.subplots_adjust(hspace=0)

# saving stuff
fig10.savefig("test10.png")
print("figure 10 done!")

### save all figs as vectors

fig1.savefig("std1.pdf", dpi = 600)
fig2.savefig("std2.pdf", dpi = 600)
fig3.savefig("pot.pdf", dpi = 600)
fig5.savefig("fricnu.pdf", dpi = 600)
fig6.savefig("incom.pdf", dpi = 600)
fig6a.savefig("incomfq.pdf",dpi = 600)
fig6b.savefig("incompot.pdf",dpi = 600)
fig6a2.savefig("incomfq2.pdf",dpi = 600)
fig6b2.savefig("incompot2.pdf",dpi = 600)
fig7.savefig("vdep.pdf", dpi = 600)
fig7a.savefig("vdepfq.pdf", dpi = 600)
fig7c.savefig("vdep2.pdf", dpi = 600)
fig7b.savefig("vdepot.pdf", dpi = 600)
fig7b2.savefig("vdepot2.pdf", dpi = 600)
fig8.savefig("sketch.pdf", dpi = 600)
fig8.savefig("sketch.pdf", dpi = 600)
fig9.savefig("irrsub.pdf", dpi = 600)
fig10.savefig("thermal.pdf", dpi = 600)

print("exports done!")




