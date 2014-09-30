from numpy import *
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.colors import LogNorm,Normalize
import sys
import math
import h5py
import numpy as np

filename = sys.argv[1]

f = h5py.File(filename,'r')
Data = f['Data']

Phi_arr = array(Data[:,0])
Radius_arr = array(Data[:,1])
xpoints_arr = np.multiply(Radius_arr,np.cos(Phi_arr)) + np.random.rand(Radius_arr.size)/500.
ypoints_arr = np.multiply(Radius_arr,np.sin(Phi_arr)) + np.random.rand(Radius_arr.size)/500.


data_arr = array(Data[:,8])
print "done reading data in"
f.close()

min_radius = 0.0

# Create the Triangulation; no triangles so Delaunay triangulation created.                  
triang = tri.Triangulation(xpoints_arr, ypoints_arr)

# Mask off unwanted triangles. 
xmid = xpoints_arr[triang.triangles].mean(axis=1)
ymid = ypoints_arr[triang.triangles].mean(axis=1)
mask = np.where(xmid*xmid + ymid*ymid < min_radius*min_radius, 1, 0)
triang.set_mask(mask)

# pcolor plot.                                                                               
plt.figure(1)

#log scale
#log_datamin=-6.
#log_datamax=0.1
#datamin=np.power(10,log_datamin)
#datamax=np.power(10,log_datamax)
#data_arr[np.where(data_arr>datamax)] = datamax-1.e-10
#data_arr[np.where(data_arr<datamin)] = datamin+1.e-10
#v = np.power(10.,np.linspace(log_datamin,log_datamax,400))

datamin=0.0
datamax=1.1
v = np.linspace(datamin,datamax,200)

plt.xlim(-1.,1.)
plt.ylim(-1.,1.)
plt.gca().set_aspect('equal')
#plt.tricontourf(triang, data_arr,levels=v,norm=LogNorm(vmin=v.min(),vmax=v.max()),cmap=plt.cm.gist_heat_r)
plt.tricontourf(triang, data_arr,levels=v,norm=Normalize(vmin=v.min(),vmax=v.max()),cmap=plt.cm.gist_heat_r)
plt.colorbar(ticks=[0.,0.2,0.4,0.6,0.8,1.0])


plt.savefig('vortex2D.png')


