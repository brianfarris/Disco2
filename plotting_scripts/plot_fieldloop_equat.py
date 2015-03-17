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
diagnum = int(sys.argv[3])

f = h5py.File(filename,'r')
Data = f['EQUAT']

Phi_arr = array(Data[:,1])
Radius_arr = array(Data[:,0])
xpoints_arr = np.multiply(Radius_arr,np.cos(Phi_arr)) 
ypoints_arr = np.multiply(Radius_arr,np.sin(Phi_arr))

#data_arr = np.log10(array(Data[:,3]))
#data_arr = .75/256*np.divide(array(Data[:,15]),np.sqrt(array(Data[:,13])))
#data_arr = 0.75/256*array(Data[:,15])
#data_arr = np.sqrt(np.array(Data[:,13]))
data_arr = np.array(Data[:,diagnum])
print "done reading data in"
f.close()

# Creating a Triangulation without specifying the triangles results in the                   
# Delaunay triangulation of the points.                                                      

min_radius = 0.00001

# Create the Triangulation; no triangles so Delaunay triangulation created.                  
triang = tri.Triangulation(xpoints_arr, ypoints_arr)

# Mask off unwanted triangles. 
xmid = xpoints_arr[triang.triangles].mean(axis=1)
ymid = ypoints_arr[triang.triangles].mean(axis=1)
mask = np.where(xmid*xmid + ymid*ymid < min_radius*min_radius, 1, 0)
triang.set_mask(mask)


datamin = -0.001;
datamax = 0.0;
    
# pcolor plot.                                                                               
plt.figure(1)
v = np.linspace(datamin,datamax, 200, endpoint=True)
#plt.subplot(211)
plt.xlim(-0.75,.75)
plt.ylim(-.75,.75)
plt.gca().set_aspect('equal')
plt.tricontourf(triang, data_arr,200,cmap=plt.cm.gist_heat_r)
plt.colorbar()

#plt.subplot(212)

if (plottype=="screen"):
    plt.show()
else:
    save_filename = filename+'.png'
    plt.savefig(save_filename)

