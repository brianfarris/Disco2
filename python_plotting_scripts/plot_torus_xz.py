import sys
import math 
import h5py
from pylab import *
from scipy import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

filename = sys.argv[1]

#open file
f = h5py.File(filename,'r')
Data = f['/Data']

#setup interpolation points
xi = np.linspace(0.4,2.0,200)
yi = array([0.0])
zi = np.linspace(-1.0,1.0,200)

#read in coordinates
Phi_arr = array(Data[:,0])
Radius_arr = array(Data[:,1])
Z_arr = array(Data[:,2])

#convert to cartesian
xpoints_arr = np.multiply(Radius_arr,np.cos(Phi_arr))
ypoints_arr = np.multiply(Radius_arr,np.sin(Phi_arr))
zpoints_arr = Z_arr
data_arr = array(Data[:,3])

#close file
f.close()

#choose subset of points within some distance of meridional plane
slicewidth = 0.1
merid_slice_index = np.where(np.fabs(ypoints_arr)<slicewidth)
xpoints_arr = xpoints_arr[merid_slice_index]
ypoints_arr = ypoints_arr[merid_slice_index]
zpoints_arr = zpoints_arr[merid_slice_index]
data_arr = data_arr[merid_slice_index]

#interpolate to regular grid
gridfunc = griddata((xpoints_arr,ypoints_arr,zpoints_arr),data_arr,(xi[None,:],yi[:,None],zi[:,None]))

#plotting range
v = np.linspace(0.0,0.6, 200, endpoint=True)

#plot
CS = plt.contourf(xi,zi,gridfunc,v,cmap=plt.cm.afmhot_r,extend='both')
plt.colorbar(CS) # draw colorbar
plt.xlim(0.4,2.0)
plt.ylim(-1.0,1.0)
plt.title(filename)
plt.show()

