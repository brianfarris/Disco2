import math
from itertools import imap, izip
import h5py as h5
import numpy as np

def readCheckpoint(filename):

    f = h5.File(filename, "r")
    Data = f['Data'][...]
    t = f['T'][0]
    f.close()

    #Standard way of removing duplicates (thnks to dfm and dh)
    Data = np.array(dict(izip(imap(hash, imap(tuple,Data[:,0:3])), Data)).values())
    
    #That step messed up the data, sort it by z, r, phi
    sort_inds = np.lexsort((Data[:,0], Data[:,1], Data[:,2]))
    Data = Data[sort_inds]

    #read in coords
    piph = np.array(Data[:,0])
    r = np.array(Data[:,1])
    z = np.array(Data[:,2])

    #read in data
    rho = np.array(Data[:,3])
    P = np.array(Data[:,4])
    vr = np.array(Data[:,5])
    vp = np.array(Data[:,6])
    vz = np.array(Data[:,7])
    try:
        q = np.array(Data[:,8])
    except:
        q = np.zeros(r.shape)

    #calculate cell volumes
    z_vals = np.unique(z)
    r_vals = np.unique(r)

    phi = np.ones(len(r))
    dV = np.ones(len(r))
    dr = np.ones(len(r))
    dphi = np.ones(len(r))

    #dz
    if len(z_vals) > 1:
        dV[:] = (z_vals[len(z_vals)-1] - z_vals[0]) / (len(z_vals) - 1.0)

    #dr
    dr_vals = np.ones(r_vals.shape)
    Rm = 0.5*(r_vals[-2]+r_vals[-1])
    Rp = 1.5*r_vals[-1] - 0.5*r_vals[-2]
    for i in range(len(r_vals)-1, -1, -1):
        Rm = 2*r_vals[i] - Rp
        dr_vals[i] = Rp - Rm
        Rp = Rm

    for i in range(len(r_vals)):
        dV[r==r_vals[i]] *= dr_vals[i]
        dr[r==r_vals[i]] = dr_vals[i]

    #dphi
    for i in range(len(r_vals)):
        inds = (r==r_vals[i])
        my_phi = piph[inds]
        dp = 2*math.pi/len(my_phi)
        dphi[inds] = dp
        phi[inds] = piph[inds] - 0.5*dp
        dV[inds] *= r[inds] * dphi[inds]

    return t, r, phi, z, rho, P, vr, vp, vz, dV, q, piph

