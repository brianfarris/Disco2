import math
from numpy import *
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import sys
import h5py
import numpy as np
import readChkpt as rc

if(len(sys.argv) < 2):
    print("\nGive me a checkpoint (.h5) file.\n")
    sys.exit()

filename = sys.argv[1]

print("Loading Data...")
t,r,phi,z,rho,P,vr,vp,vz,dV = rc.readChkpt(filename)

N = len(r)

print("Loaded.\nT: {0:.12g}".format(t))

xmin = 2.0
xmax = 4.0
ymin = -2.0
ymax = 2.0

#convert to cartesian
x = r * np.cos(phi)
y = r * np.sin(phi)
vx = vr * np.cos(phi) - r*vp * np.sin(phi)
vy = vr * np.sin(phi) + r*vp * np.cos(phi)

num_colours = 1000

triang = tri.Triangulation(x,y)
plt.figure()
plt.subplot(221)
plt.tricontourf(triang, rho, num_colours, cmap=plt.cm.afmhot_r)
plt.colorbar()
plt.xlabel(r"$x$")
plt.ylabel(r"$y$")
plt.subplot(222)
plt.tricontourf(triang, P, num_colours, cmap=plt.cm.afmhot_r)
plt.colorbar()
plt.xlabel(r"$x$")
plt.ylabel(r"$y$")
plt.subplot(223)
plt.tricontourf(triang, vr, num_colours, cmap=plt.cm.afmhot_r)
plt.colorbar()
plt.xlabel(r"$x$")
plt.ylabel(r"$y$")
plt.subplot(224)
plt.tricontourf(triang, vp, num_colours, cmap=plt.cm.afmhot_r)
plt.colorbar()
plt.xlabel(r"$x$")
plt.ylabel(r"$y$")
plt.savefig("disc.png")

plt.show()

