from numpy import *
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import sys
import h5py
import numpy as np
import discopy as dp

equat = False
if(len(sys.argv) < 3):
    print("\nGive me a checkpoint (.h5) file and a parameter to plot.\n")
    sys.exit()

filename = sys.argv[1]
diagnum = int(sys.argv[2])

if(len(sys.argv) > 3):
    if(sys.argv[3] == 'e'):
        equat = True

dp = rc.readCheckpoint(filename)
t = dat[0]
r = dat[1]
phi = dat[2]
z = dat[3]

#open hdf5 file
#f = h5py.File(filename,'r')
#Data = f['Data']

#read in coords
#phi = array(Data[:,0])
#r = np.array(Data[:,1])
#z = np.array(Data[:,2])

#Get Data
data = array(dat[diagnum+1])

#close file
#f.close()

#convert to cartesian
x = r * np.cos(phi)
y = r * np.sin(phi)

# Create the Triangulation; no triangles so Delaunay triangulation created.     
if(equat):
    inds = z==z[len(z)/2]
    triang = tri.Triangulation(x[inds], y[inds])

    #set plot range
    spread = data[inds].max() - data[inds].min()
    v = np.linspace(data[inds].min()-0.1*spread,data[inds].max()+0.1*spread, 200, endpoint=True)

#Plot.
plt.figure()
plt.subplot(221)
plt.plot(x, data, 'k.')
plt.xlabel("x")
plt.subplot(222)
plt.plot(y, data, 'k.')
plt.xlabel("y")
plt.subplot(223)
plt.plot(z, data, 'k.')
plt.xlabel("z")
plt.savefig("prim_plot.png")

if(equat):
    plt.figure()
    plt.xlim(x.min(),x.max())
    plt.ylim(y.min(),y.max())
    plt.gca().set_aspect('equal')
    plt.tricontourf(triang, data[inds],v,cmap=plt.cm.afmhot_r)
    plt.colorbar()
    plt.savefig("prim_plot_equat.png")

plt.show()

