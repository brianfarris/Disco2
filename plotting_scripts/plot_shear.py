from numpy import *
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import sys
import math
import h5py
import numpy as np
import re

# read filename from command line
filename = sys.argv[1]
print "filename: ",filename

# parse filename to get time
time_string = re.findall("\d+.\d+", filename)
time = float(time_string[0])
print time
t0=0.25

# open file
f = h5py.File(filename,'r')
Data = f['/EQUAT']

# set data arrays
Phi_arr = array(Data[:,1])
Radius_arr = array(Data[:,0])
xpoints_arr = np.multiply(Radius_arr,np.cos(Phi_arr))
ypoints_arr = np.multiply(Radius_arr,np.sin(Phi_arr))
vr_arr = array(Data[:,4])
vp_arr = array(Data[:,5])
vx_arr = np.multiply(-vr_arr,np.cos(Phi_arr))+np.multiply(vp_arr,np.sin(Phi_arr))
vy_arr = np.multiply(vr_arr,np.sin(Phi_arr))+np.multiply(vp_arr,np.cos(Phi_arr))

# we only care about data near the x-axis
yslice_width = 0.1
yslice_pos = 3.0
xpoints_arr = xpoints_arr[np.where(abs(ypoints_arr-yslice_pos)<yslice_width)]
vx_arr = vx_arr[np.where(abs(ypoints_arr-yslice_pos)<yslice_width)]
vy_arr = vy_arr[np.where(abs(ypoints_arr-yslice_pos)<yslice_width)]

print "done reading data in"
f.close()

#get analytic solution
x_analytic = np.array(linspace(0,10,400))
vy_analytic = np.array(linspace(0,10,400))
vy_analytic_0 = np.array(linspace(0,10,400))
nu = 0.05 # hard-coded viscosity
t=t0+time
for i in range(len(x_analytic)):
  vy_analytic[i] = 1./sqrt(2*math.pi*nu*t)*exp(-(x_analytic[i]-3.)*(x_analytic[i]-3.)/(4.*nu*t))
  vy_analytic_0[i] = 1./sqrt(2*math.pi*nu*t0)*exp(-(x_analytic[i]-3.)*(x_analytic[i]-3.)/(4.*nu*t0))

#make actual plot
fig,ax = plt.subplots(figsize=(8,8))
ax.set_xlim(0.5,6.)
ax.plot(x_analytic,vy_analytic_0,linewidth=3.,color='black',label=r'$t/t_0=0.0$')
ax.plot(xpoints_arr[np.argsort(xpoints_arr)],vy_arr[np.argsort(xpoints_arr)],color='blue',linewidth=3.0,label=r'simulation $t/t_0 = '+str(time/t0)+'$')
ax.plot(x_analytic,vy_analytic,linewidth=3.,color='red',label=r'analytic $t/t_0='+str(time/t0)+'$')
ax.set_xlabel(r'$x$',size='x-large')
ax.set_ylabel(r'$v^y$',rotation='horizontal',size='x-large')

# create legend
legend = ax.legend(loc='upper right', shadow=True)

plt.savefig('shear1D.pdf')
