from numpy import *
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import sys
import h5py
import numpy as np

filename = sys.argv[1]

#open hdf5 file
f = h5py.File(filename,'r')
Data = f['EQUAT']

#read in radius and phi
Radius_arr = array(Data[:,0])
Phi_arr = array(Data[:,1])

#convert to cartesian
xpoints_arr = np.multiply(Radius_arr,np.cos(Phi_arr)) 
ypoints_arr = np.multiply(Radius_arr,np.sin(Phi_arr))

dens_arr = array(Data[:,2])
vr_arr = array(Data[:,3])
vphi_arr = array(Data[:,4])
pres_arr = array(Data[:,5])

#close file
f.close()

vx_arr = vr_arr * np.cos(Phi_arr) - vphi_arr * np.sin(Phi_arr)
vy_arr = vr_arr * np.sin(Phi_arr) + vphi_arr * np.cos(Phi_arr)

print("Plotting!\n")

plt.figure()
plt.subplot(221)
plt.plot(xpoints_arr, dens_arr, 'k.')
plt.ylabel("Density")
plt.subplot(222)
plt.plot(xpoints_arr, pres_arr, 'k.')
plt.ylabel("Pressure")
plt.subplot(223)
plt.plot(xpoints_arr, vx_arr, 'k.')
plt.ylabel("Velocity X")
plt.subplot(224)
plt.plot(xpoints_arr, vy_arr, 'k.')
plt.ylabel("Velocity Y")

plt.show()

