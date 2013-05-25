#from pylab import *
#from scipy import *
from numpy import *
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import sys
import math
import h5py
import numpy as np
##from scipy.interpolate import griddata

filename = sys.argv[1]
plottype = sys.argv[2]

f = h5py.File(filename,'r')
Data = f['Data']

Phi_arr = array(Data[:,0])
Radius_arr = array(Data[:,1])
xpoints_arr = np.multiply(Radius_arr,np.cos(Phi_arr)) + np.random.rand(Radius_arr.size)/500.
ypoints_arr = np.multiply(Radius_arr,np.sin(Phi_arr)) + np.random.rand(Radius_arr.size)/500.

data_arr = array(Data[:,8])
print "done reading data in"
f.close()

# Creating a Triangulation without specifying the triangles results in the                   
# Delaunay triangulation of the points.                                                      

min_radius = 0.0

# Create the Triangulation; no triangles so Delaunay triangulation created.                  
triang = tri.Triangulation(xpoints_arr, ypoints_arr)

# Mask off unwanted triangles. 
xmid = xpoints_arr[triang.triangles].mean(axis=1)
ymid = ypoints_arr[triang.triangles].mean(axis=1)
mask = np.where(xmid*xmid + ymid*ymid < min_radius*min_radius, 1, 0)
triang.set_mask(mask)


datamin = -0.0;
datamax = 0.75;
    
# pcolor plot.                                                                               
plt.figure(1)
v = np.linspace(datamin,datamax, 200, endpoint=True)
plt.xlim(-1.,1.)
plt.ylim(-1.,1.)
plt.gca().set_aspect('equal')
plt.tricontourf(triang, data_arr,200,cmap=plt.cm.gist_heat_r)
plt.colorbar()

if (plottype=="screen"):
    plt.show()
else:
    save_filename = filename+'.png'
    plt.savefig(save_filename)

