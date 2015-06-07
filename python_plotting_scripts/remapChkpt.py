import sys
import math
import numpy as np
import scipy.optimize as opt
import h5py as h5
import readChkpt as rc
import readParfile as rp

class Grid:
    # A class to contain a Disco grid.

    nr_noghost = 0
    nz_noghost = 0
    nr_tot = 0
    nz_tot = 0
    ng_rmin = 0
    ng_rmax = 0
    ng_zmin = 0
    ng_zmax = 0
    np = None

    rFaces = None
    zFaces = None
    pFaces = None

    prim = None
    grad = None

    def __init__(self, pars):
        self._setNRZ(pars)
        self._setFacePosRZ(pars)
        self._setNP(pars)
        self._setFacePosP(pars)

    def loadChkpt(self, filename):
        dat = rc.readChkpt(filename)
        r = dat[1]
        z = dat[3]
        piph = dat[11]

        if self._checkSame(r, z, piph):
            print("Looks like the same grid! Copying...")
            self._loadChkptIdentical(dat)
        else:
            print("Looks like different grids! Remapping...")
            self._loadChkptRemap(dat)

    def _setNRZ(self, pars):
        # Same logic as in sim_alloc_arr() in sim_alloc_arr.c

        nr = pars["NumR"]
        nz = pars["NumZ"]
        ng = pars["ng"]

        self.nr_noghost = nr
        self.nz_noghost = nz
        self.ng = ng

        if pars["NoInnerBC"] == 0:
            self.nr_tot = 2*ng + nr
            self.ng_rmin = ng
            self.ng_rmax = ng
        else:
            self.nr_tot = ng + nr + 1
            self.ng_rmin = 1
            self.ng_rmax = ng

        if nz == 1:
            self.nz_tot = 1
            self.ng_zmin = 0
            self.ng_zmax = 0
        else:
            self.nz_tot = nz + 2*ng
            self.ng_zmin = ng
            self.ng_zmax = ng

    def _setFacePosRZ(self, pars):
        # This is meant to be an exact copy of the set_rz() and set_N_p() 
        # functions in sim_set.c

        Rmin = pars["R_Min"]
        Rmax = pars["R_Max"]
        Zmin = pars["Z_Min"]
        Zmax = pars["Z_Max"]

        Rscale = pars["RLogScale"]
        Zscale = pars["ZLogScale"]
        sigma = pars["HiResSigma"]
        R0 = pars["HiResR0"]
        fac = pars["HiResFac"]

        r2_max = Rmax + fac*sigma*math.sqrt(math.pi/4.0) \
                            * (math.erf((Rmax-R0)/sigma) + 1.0)
        r2_min = Rmin + fac*sigma*math.sqrt(math.pi/4.0) \
                            * (math.erf((Rmin-R0)/sigma) + 1.0)
        r1_max = Rscale * math.log(1.0 + r2_max/Rscale)
        r1_min = Rscale * math.log(1.0 + r2_min/Rscale)
        z1_max = np.sign(Zmax) * Zscale * math.log(1.0+math.fabs(Zmax)/Zscale)
        z1_min = np.sign(Zmin) * Zscale * math.log(1.0+math.fabs(Zmin)/Zscale)

        deltaR = (r1_max - r1_min) / self.nr_noghost
        deltaZ = (z1_max - z1_min) / self.nz_noghost
        
        r1 = r1_min + deltaR * (np.arange(self.nr_tot+1) - self.ng_rmin)
        z1 = z1_min + deltaZ * (np.arange(self.nz_tot+1) - self.ng_zmin)

        r2 = Rscale * (np.exp(r1/Rscale)-1.0)

        self.rFaces = np.zeros(self.nr_tot+1)
        self.zFaces = np.zeros(self.nz_tot+1)

        rTol = 1.0e-12
        fTol = 2*np.finfo(float).eps

        for i, r2val in enumerate(r2):
            self.rFaces[i] = opt.brentq(_r_func, 0.0, 2*Rmax, maxiter=1000,
                                            args=(r2val,fac,sigma,R0))

        if pars["NoInnerBC"] != 0:
            self.rFaces[0] = 0.0

        self.zFaces = np.sign(z1) * Zscale * (np.exp(np.fabs(z1)/Zscale)-1.0)

    def _setNP(self, pars):
        self.np = np.zeros((self.nz_tot,self.nr_tot), dtype=np.int64)
        NP = pars["NP_CONST"]
        aspect = pars["aspect"]

        if NP > 0:
            self.np[:,:] = NP
        else:
            r = self.rFaces[1:]
            dr = self.rFaces[1:]-self.rFaces[:-1]
            self.np[:,:] = [int(x) for x in 2*math.pi*(1.0 + (r/dr-1.0)/aspect)]

    def _setFacePosP(self, pars):

        phi0 = np.random.rand(self.nz_tot,self.nr_tot)

        self.pFaces = []

        for k in xrange(self.nz_tot):
            slice = []
            for i, npi in enumerate(self.np[k]):
                dphi = 2*math.pi / npi
                piph = dphi * (phi0[k,i] + np.arange(npi))
                piph[piph >= 2*math.pi] -= 2*math.pi
                slice.append(piph)
            self.pFaces.append(slice)

    def _checkSame(self, R, Z, Piph):
        # If the grid data (R, Z, Piph) are identical
        # in structure to the current grid return True, 
        # otherwise return False.
        # Assumes R,Z are cell-centered, not face-centered.

        z_vals = np.unique(Z)
        r_vals = np.unique(R)

        if len(z_vals) != self.nz_tot:
            return False
        if len(r_vals) != self.nr_tot:
            return False

        for k,z in enumerate(z_vals):

            if z < self.zFaces[k] or z > self.zFaces[k+1]:
                return False

            zind = (Z==z)
            Rz = R[zind]
            Zz = Z[zind]
            Piphz = Piph[zind]

            for i,r in enumerate(r_vals):
            
                if r < self.rFaces[i] or r > self.rFaces[i+1]:
                    return False

                rind = (R==r)
                if len(Piphz[rind]) != self.np[k][i]:
                    return False

        return True

    def _loadChkptIdentical(self, dat):

        if self.prim is None:
            self.prim = []

        R = dat[1]
        Z = dat[3]
        Piph = dat[11]
        prim = np.array([dat[4], dat[5], dat[6], dat[7], dat[8], dat[10]])
        prim = prim.T

        Rvals = np.unique(R)
        Zvals = np.unique(Z)

        for k, z in enumerate(Zvals):
            
            slice = []
            zind = z==Z
            
            for i, r in enumerate(Rvals):
                rind = r==R
                self.pFaces[k][i] = Piph[zind][rind]

                #annulus = np.zeros((self.np[k][i], prim.shape[0]))
                annulus = prim[zind][rind]
                slice.append(annulus)
            
            self.prim.append(slice)

    def _loadChkptRemap(self, dat):
        print("Not Same!")
        return

def _r_func(r, r2, fac, sigma, r0):
    #print r, r2, fac, sigma, r0
    f = r - r2 + fac*sigma*math.sqrt(0.25*math.pi)*(math.erf((r-r0)/sigma)+1.0)
    #print f
    return f


if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("I need a parfile to nom.")
        sys.exit()

    elif len(sys.argv) == 2:
        pars = rp.readParfile(sys.argv[1])
        g = Grid(pars)

        print g.rFaces
        print g.zFaces
        print g.pFaces

    elif len(sys.argv) > 2:
        
        pars = rp.readParfile(sys.argv[1])
        g = Grid(pars)
        g.loadChkpt(sys.argv[2])



