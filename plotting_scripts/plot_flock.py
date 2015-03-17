import sys
import math 
import h5py
from pylab import *
from scipy import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

def makeplot(xi,zi,data,plotsizex,plotsizez,datamin,datamax,filename):
    v = np.linspace(datamin,datamax, 300, endpoint=True)
    CS = plt.contourf(xi,zi,data,v,cmap=plt.cm.BrBG)
    plt.colorbar() # draw colorbar
    plt.xlim(1.5,plotsizex)
    plt.ylim(-plotsizez,plotsizez)
    plt.title(filename)
    save_filename = filename+'.png'
    plt.savefig(save_filename)

def interp_data(xi,yi,zi,Data,index):
    slicewidth = 0.1
    Phi_arr = array(Data[:,0])
    Radius_arr = array(Data[:,1])
    Z_arr = array(Data[:,2])
        
    xpoints_arr = np.multiply(Radius_arr,np.cos(Phi_arr))
    ypoints_arr = np.multiply(Radius_arr,np.sin(Phi_arr))
    zpoints_arr = Z_arr
    
    data_arr = array(Data[:,index])
    merid_slice_index = np.where(np.fabs(ypoints_arr)<slicewidth)
    xpoints_arr = xpoints_arr[merid_slice_index]
    ypoints_arr = ypoints_arr[merid_slice_index]
    zpoints_arr = zpoints_arr[merid_slice_index]
    data_arr = data_arr[merid_slice_index]
    print "done reading data in"
    gridfunc = griddata((xpoints_arr,ypoints_arr,zpoints_arr),data_arr,(xi[None,:],yi[:,None],zi[:,None]))
    return gridfunc

filename = sys.argv[1]
plotsizex = 4.0
plotsizez = 0.5
index = 9

f = h5py.File(filename,'r')
Data = f['Data']

#setup interpolation points
xi = np.linspace(1.5,plotsizex,300)
yi = array([0.0])
zi = np.linspace(-plotsizez,plotsizez,300)
datai = interp_data(xi,yi,zi,Data,index)
print "done interpolating"
f.close()

datamin = -np.fabs(datai).max()
datamax = np.fabs(datai).max();
makeplot(xi,zi,datai,plotsizex,plotsizez,datamin,datamax,filename)
