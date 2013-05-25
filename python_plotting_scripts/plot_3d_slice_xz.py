import sys
import math 
import h5py
from pylab import *
from scipy import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

def torus_analytic_rho(x,y,z,C,lk,Gamma,K):
    R = sqrt(x*x+y*y)
    PoRho = (Gamma-1.0)/(Gamma*K)*(C+1.0/sqrt(R*R+z*z)-1./2.*lk*lk/R/R)
    n = 1./(Gamma-1.0)
    if (PoRho>0):
        rho = (PoRho/K)**n
    else:
        rho=1.0e-5
    return(rho)

def makeplot(xi,zi,data,plotsizex,plotsizez,datamin,datamax,filename):
    origin = 'upper'
    v = np.linspace(datamin,datamax, 300, endpoint=True)
    CS = plt.contourf(xi,zi,data,v,cmap=plt.cm.binary,extend='both')
#    CS.cmap.set_over('black')
#    CS.cmap.set_over('white')
    plt.colorbar(CS) # draw colorbar
    plt.xlim(0.4,plotsizex)
    plt.ylim(-plotsizez,plotsizez)
    plt.title(filename)
    save_filename = filename+'.png'
    plt.savefig(save_filename)

def interp_data(xi,yi,zi,Data,index,datatype,norm,grid0):
    slicewidth = 0.1
    Phi_arr = array(Data[:,0])
    Radius_arr = array(Data[:,1])
    Z_arr = array(Data[:,2])
        
    xpoints_arr = np.multiply(Radius_arr,np.cos(Phi_arr))
    ypoints_arr = np.multiply(Radius_arr,np.sin(Phi_arr))
    zpoints_arr = Z_arr
    if datatype=="vector":
        Vi = array(Data[:,index])
        Vj = array(Data[:,index+1])
        Vk = array(Data[:,index+2])
        data_arr = np.sqrt(np.multiply(Vi,Vi)+np.multiply(Vj,Vj)+np.multiply(Vk,Vk))
    else:
        #data_arr_denom = array(Data[:,index])
        #data_arr_num = array(Data[:,index+1])
        #data_arr = np.divide(data_arr_num,data_arr_denom)
        #data_arr = 0.01*np.power(array(Data[:,index]),5./3.)
        data_arr = array(Data[:,index])
    merid_slice_index = np.where(np.fabs(ypoints_arr)<slicewidth)
    xpoints_arr = xpoints_arr[merid_slice_index]
    ypoints_arr = ypoints_arr[merid_slice_index]
    zpoints_arr = zpoints_arr[merid_slice_index]
    data_arr = data_arr[merid_slice_index]
    rho_arr_analytic = np.zeros(len(data_arr))
    for i in arange(len(data_arr)):
        rho_arr_analytic[i] = torus_analytic_rho(xpoints_arr[i],ypoints_arr[i],zpoints_arr[i],-0.1875,1.58113883008,5./3.,0.01)
    print "done reading data in"
    gridfunc = griddata((xpoints_arr,ypoints_arr,zpoints_arr),data_arr,(xi[None,:],yi[:,None],zi[:,None]))/norm-grid0
    gridfunc_analytic = griddata((xpoints_arr,ypoints_arr,zpoints_arr),rho_arr_analytic,(xi[None,:],yi[:,None],zi[:,None]))/norm-grid0
    #error = np.divide(np.fabs(gridfunc-gridfunc_analytic),(np.fabs(gridfunc)+np.fabs(gridfunc_analytic)))
    error = gridfunc-gridfunc_analytic
    return gridfunc

filename = sys.argv[1]
plotsizex = double(sys.argv[2])
plotsizez = double(sys.argv[3])
norm = double(sys.argv[4])
grid0 = double (sys.argv[5])
datatype = sys.argv[6]
index = int(sys.argv[7])
logscale = int(sys.argv[8])

f = h5py.File(filename,'r')
Data = f['/Data']

#print np.size(Data)
#for i in arange(np.size(Data)):
#    print Data[i,2],Data[i,3]

#setup interpolation points
xi = np.linspace(0.4,plotsizex,300)
yi = array([0.0])
zi = np.linspace(-plotsizez,plotsizez,300)
datai = interp_data(xi,yi,zi,Data,index,datatype,norm,grid0)
print "done interpolating"
f.close()

datamin = 0.0
datamax = 0.6;
if (logscale==1):
    datai=np.log10(np.fabs(datai))
    datamin = math.log10(datamin)
    datamax = math.log10(datamax)
makeplot(xi,zi,datai,plotsizex,plotsizez,datamin,datamax,filename)
