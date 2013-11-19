from numpy import *
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import sys
import h5py
import numpy as np

if(len(sys.argv) < 2):
    print("\nGive me a checkpoint (.h5) file.\n")
    sys.exit()

filename = sys.argv[1]

#open hdf5 file
f = h5py.File(filename,'r')
Data = f['Data']

#read in coords
phi = array(Data[:,0])
r = array(Data[:,1])
z = array(Data[:,2])

#read in data
rho = np.array(Data[:,3])
P = np.array(Data[:,4])
vr = np.array(Data[:,5])
vphi = np.array(Data[:,6])
vz = np.array(Data[:,7])

#close file
f.close()

#convert to cartesian
x = r * np.cos(phi)
y = r * np.sin(phi)
vx = vr * cos(phi) - r * vphi * sin(phi) 
vy = vr * sin(phi) + r * vphi * cos(phi) 
equat = z==z[len(z)/2]

#triangulate!
triang = tri.Triangulation(x[equat], y[equat]);
spread = np.linspace(rho[equat].min(), rho[equat].max(), 200, endpoint=True)

# plot.                                                                               
plt.figure()
plt.subplot(221)
plt.plot(x, rho, 'k.')
plt.xlabel(r'$x$')
plt.ylabel(r'$\rho_0$')
plt.subplot(222)
plt.plot(x, P, 'k.')
plt.xlabel(r'$x$')
plt.ylabel(r'$P$')
plt.subplot(223)
plt.plot(x, vx, 'k.')
plt.xlabel(r'$x$')
plt.ylabel(r'$v^x$')

plt.subplot(224)
plt.xlim(x.min(), x.max())
plt.ylim(y.min(), y.max())
plt.gca().set_aspect('equal')
plt.tricontourf(triang, rho[equat], spread, cmap=plt.cm.afmhot_r)
plt.colorbar()
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')

plt.savefig("prim_plot_x.png")

plt.show()

