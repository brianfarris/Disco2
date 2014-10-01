from numpy import *
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import sys
import math
import h5py
import numpy as np

filename = sys.argv[1]

f = h5py.File(filename,'r')
Data = f['EQUAT']

Phi_arr = array(Data[:,1])
Radius_arr = array(Data[:,0])
xpoints_arr = np.multiply(Radius_arr,np.cos(Phi_arr)) 
ypoints_arr = np.multiply(Radius_arr,np.sin(Phi_arr))

vr_arr = array(Data[:,4])
vp_arr = array(Data[:,5])
vy_arr = np.multiply(np.cos(array(Phi_arr)),vp_arr) + np.multiply(np.sin(array(Phi_arr)),vr_arr)
print "done reading data in"
f.close()

# Create the Triangulation; no triangles so Delaunay triangulation created.                  
triang = tri.Triangulation(xpoints_arr, ypoints_arr)

# Mask off unwanted triangles. 
min_radius = 0.00001
xmid = xpoints_arr[triang.triangles].mean(axis=1)
ymid = ypoints_arr[triang.triangles].mean(axis=1)
mask = np.where(xmid*xmid + ymid*ymid < min_radius*min_radius, 1, 0)
triang.set_mask(mask)

datamin = 0.0;
datamax = 0.5;
    
fig,ax = plt.subplots()
v = np.linspace(datamin,datamax, 200, endpoint=True)
ax.set_xlim(-10.,10.)
ax.set_ylim(-10.,10.)
cax = ax.tricontourf(triang, vy_arr,200,cmap=plt.cm.gist_heat_r)
cbar = fig.colorbar(cax)


plt.savefig('shear2D.png')


