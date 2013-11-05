from numpy import *
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import sys
import h5py
import numpy as np

if(len(sys.argv) < 3):
    print("\nGive me a snapshot (.h5) file and a parameter to plot.\n")
    sys.exit()

filename = sys.argv[1]
diagnum = int(sys.argv[2])

#open hdf5 file
f = h5py.File(filename,'r')
Data = f['EQUAT']

#read in radius and phi
Radius_arr = array(Data[:,0])
Phi_arr = array(Data[:,1])

#convert to cartesian
xpoints_arr = np.multiply(Radius_arr,np.cos(Phi_arr)) 
ypoints_arr = np.multiply(Radius_arr,np.sin(Phi_arr))

data_arr = array(Data[:,diagnum])

#close file
f.close()

# Create the Triangulation; no triangles so Delaunay triangulation created.                  
triang = tri.Triangulation(xpoints_arr, ypoints_arr)

#set plot range
spread = data_arr.max() - data_arr.min()
v = np.linspace(data_arr.min()-0.1*spread,data_arr.max()+0.1*spread, 200, endpoint=True)

# plot.                                                                               
plt.figure(1)
plt.xlim(-1.0,1.0)
plt.ylim(-1.0,1.0)
plt.gca().set_aspect('equal')
plt.tricontourf(triang, data_arr,v,cmap=plt.cm.afmhot_r)
plt.colorbar()

plt.show()

