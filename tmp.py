
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

### Trying things simple

fig3 = plt.figure(constrained_layout=True,figsize=(7,3))
gs = fig3.add_gridspec(2, 4)
f3_ax1 = fig3.add_subplot(gs[0:, 2:])
#f3_ax1.set_title('gs[0, :]')
f3_ax2 = fig3.add_subplot(gs[0, :-2])
#f3_ax2.set_title('gs[1, :-1]')
f3_ax3 = fig3.add_subplot(gs[1, :-2])
#f3_ax3.set_title('gs[1:, -1]')

fig3.savefig("tmp.png")

### Trying things more hard

nx = 2
ny = 2

fig7, f7_axs = plt.subplots(ncols=nx, nrows=ny, figsize=(7,3), sharex=True)

plt.subplots_adjust(hspace=0)

gs = f7_axs[0, 1].get_gridspec()

# remove the underlying axes
#for i in range (ny) : 
#    for ax in f7_axs[0:, i] :
#        ax.remove()

f7_axs[0,1].remove()
f7_axs[1,1].remove()

#f7_axs[0,0].remove()

axbig = fig7.add_subplot(gs[0:, -1])

fig7.tight_layout() # need to comment out to make hspace = 0 to work

fig7.savefig("tmp2.png")


