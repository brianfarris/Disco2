import h5py as h5
import numpy as np

def readChkpt(filename):

    f = h5.File(filename)
    Data = f['Data']

    #read in coords
    phi = np.array(Data[:,0])
    r = np.array(Data[:,1])
    z = np.array(Data[:,2])

    #read in data
    t = f['T'][0]
    rho = np.array(Data[:,3])
    P = np.array(Data[:,4])
    vr = np.array(Data[:,5])
    vp = np.array(Data[:,6])
    vz = np.array(Data[:,7])

    #close file
    f.close()

    #calculate cell volumes
    dV = np.ones(len(r))
    z_vals = np.unique(z)
    r_vals = np.unique(r)

    #dz
    if len(z_vals) > 1:
        dV[:] = (z_vals[len(z_vals)-1] - z_vals[0]) / (len(z_vals) - 1.0)

    #dr
    Rp = 0.0
    for i in range(len(r_vals)):
        Rm = Rp
        R = r_vals[i]
        Rp = 2.0*R - Rm
        dV[r==R] *= 0.5*(Rp*Rp - Rm*Rm)

    #dphi
    for i in range(len(r_vals)):
        inds = (r==r_vals[i])
        dphi = np.zeros(len(phi[inds]))
        dphi[1:] = phi[inds][1:] - phi[inds][:-1]
        dphi[0] = phi[inds][0] - phi[inds][-1]
        dphi[dphi<0] += 2*np.pi
        dV[inds] *= dphi

    return t, r, phi, z, rho, P, vr, vp, vz, dV

