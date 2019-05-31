
# imports
from os import system
import numpy as np

#preliminaries

gridsize = 25

#init = 3.5E1    # for q^2
init = 7**-4*14E20   # for q^4

lower = 0.5
upper = 4.0

#log grid
#logup = log(init*lower)
#loglo = log(init*upper)
#diff = logup - loglo
#step = diff / float(gridsize)

#for (uint l = 0; l <= gridsize; l++)
#{
#    double knot = exp(step*l + logmin);
#    grid.push_back(knot);
#}
#
#grids.push_back(grid);

diff = upper - lower
step = diff / float (gridsize)

curr = lower

system("> nunuout.dat")
system("> nufricout.dat")

nus = np.array([]);

for k in range (0, gridsize+1) :
    nu = curr*init
    nus = np.append(nus,nu)

    system("./nu_script.wls {}".format(nu))
    
    print("current nu is {}".format(nu))
    
    curr += step

# convert to sensible format

#with open("nunuout.dat") as infile, open ("nunu.dat",'w') as outfile : 
#    for line in infile :
#        outfile.write(line.replace('*^','e'))
#
#with open("nufricout.dat") as infile, open ("nufric.dat",'w') as outfile : 
#    for line in infile :
#        outfile.write(line.replace('*^','e'))

np.savetxt("nu_out.dat",nus)



