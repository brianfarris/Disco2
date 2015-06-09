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
    nq = 0

    _pars = None

    def __init__(self, pars=None, chkpt=None, archive=None):

        if archive is not None:
            self.loadArchive(archive)
        
        elif pars is not None:
            self.loadPars(pars)

            if chkpt is not None:
                self.loadChkpt(chkpt)

    def loadPars(self, pars):
        self._pars = pars
        self._setNRZ()
        self._setFacePosRZ()
        self._setNP()
        self._setFacePosP()

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

    def loadArchive(self, filename):
        f = h5.File(filename, "r")

        parGroup = f["parameters"]
        gridGroup = f["grid"]
        dataGroup = f["data"]

        self._pars = dict()
        for key in parGroup.attrs:
            self._pars[key] = parGroup.attrs[key]

        self.rFaces = gridGroup['r_faces'][...]
        self.zFaces = gridGroup['z_faces'][...]
        self.np = gridGroup['np'][...]
        piph_arr = gridGroup['p_faces'][...][2,:]

        self.nr_tot = gridGroup.attrs['nr_tot']
        self.nz_tot = gridGroup.attrs['nz_tot']
        self.ng_rmin = gridGroup.attrs['ng_rmin']
        self.ng_rmax = gridGroup.attrs['ng_rmax']
        self.ng_zmin = gridGroup.attrs['ng_zmin']
        self.ng_zmax = gridGroup.attrs['ng_zmax']
        self.nq = dataGroup.attrs['nq']

        j = 0
        self.pFaces = []
        for k in xrange(self.nz_tot):
            slice = []
            for i in xrange(self.nr_tot):
                np = self.np[k,i]
                slice.append(piph_arr[j:j+np])
                j += np
            self.pFaces.append(slice)

        if self.nq > 0:
            j = 0
            self.prim = []
            prim = dataGroup["prim"][...][2:,:].T
            for k in xrange(self.nz_tot):
                slice = []
                for i in xrange(self.nr_tot):
                    np = self.np[k,i]
                    slice.append(prim[j:j+np,:])
                    j += self.np[k,i]
                self.prim.append(slice)

        f.close()

    def saveArchive(self, filename):
        f = h5.File(filename, "w")

        parGroup = f.create_group("parameters")
        gridGroup = f.create_group("grid")
        dataGroup = f.create_group("data")

        for key in self._pars:
            parGroup.attrs[key] = self._pars[key]

        gridGroup.create_dataset("r_faces", data=self.rFaces)
        gridGroup.create_dataset("z_faces", data=self.zFaces)
        gridGroup.create_dataset("np", data=self.np)
        piphDSet = gridGroup.create_dataset("p_faces", (3, self.np.sum()))

        gridGroup.attrs["ng_rmin"] = self.ng_rmin
        gridGroup.attrs["ng_rmax"] = self.ng_rmax
        gridGroup.attrs["ng_zmin"] = self.ng_zmin
        gridGroup.attrs["ng_zmax"] = self.ng_zmax
        gridGroup.attrs["nr_tot"] = self.nr_tot
        gridGroup.attrs["nz_tot"] = self.nz_tot

        dataGroup.attrs['nq'] = self.nq

        j = 0
        for k in xrange(self.nz_tot):
            for i in xrange(self.nr_tot):
                np = self.np[k,i]
                piphDSet[0,j:j+np] = k
                piphDSet[1,j:j+np] = i
                piphDSet[2,j:j+np] = self.pFaces[k][i]

                j += self.np[k,i]

        if self.nq > 0:
            j = 0
            primDSet = dataGroup.create_dataset("prim", (2+self.nq, self.np.sum()))
            for k in xrange(self.nz_tot):
                for i in xrange(self.nr_tot):
                    np = self.np[k,i]
                    primDSet[0,j:j+np] = k
                    primDSet[1,j:j+np] = i
                    for q in xrange(self.nq):
                        primDSet[2+q, j:j+np] = self.prim[k][i][:,q]

                    j += self.np[k,i]

        f.close()

    def _setNRZ(self):
        # Same logic as in sim_alloc_arr() in sim_alloc_arr.c

        nr = self._pars["NumR"]
        nz = self._pars["NumZ"]
        ng = self._pars["ng"]

        self.nr_noghost = nr
        self.nz_noghost = nz
        self.ng = ng

        if self._pars["NoInnerBC"] == 0:
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

    def _setFacePosRZ(self):
        # This is meant to be an exact copy of the set_rz() and set_N_p() 
        # functions in sim_set.c

        Rmin = self._pars["R_Min"]
        Rmax = self._pars["R_Max"]
        Zmin = self._pars["Z_Min"]
        Zmax = self._pars["Z_Max"]

        Rscale = self._pars["RLogScale"]
        Zscale = self._pars["ZLogScale"]
        sigma = self._pars["HiResSigma"]
        R0 = self._pars["HiResR0"]
        fac = self._pars["HiResFac"]

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

        if self._pars["NoInnerBC"] != 0:
            self.rFaces[0] = 0.0

        self.zFaces = np.sign(z1) * Zscale * (np.exp(np.fabs(z1)/Zscale)-1.0)

    def _setNP(self):
        self.np = np.zeros((self.nz_tot,self.nr_tot), dtype=np.int64)
        NP = self._pars["NP_CONST"]
        aspect = self._pars["aspect"]

        if NP > 0:
            self.np[:,:] = NP
        else:
            r = self.rFaces[1:]
            dr = self.rFaces[1:]-self.rFaces[:-1]
            self.np[:,:] = [int(x) for x in 2*math.pi*(1.0 + (r/dr-1.0)/aspect)]

    def _setFacePosP(self):

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

        self.nq = 6

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

    def _calcGrad(self):
        
        # Initialize grad array to zero.  
        grad = []
        for k in xrange(self.nz_tot):
            slice = []
            for i in xrange(self.nr_tot):
                np = self.np[k,i]
                slice.append(np.zeros((np,self.nq,3), dtype=np.float64))
            grad.append(slice)

        # Phi-direction
        for k in xrange(self.nz_tot):
            for i in xrange(self.nr_tot):
                np = self.np[k,i]
                right = np.arange(np) + 1
                right[np-1] = 0
                left = np.arange(np) - 1
                left[0] = np

                dphi = self.pFaces - self.pFaces[left]
                dphi[dphi>2*np.pi] -= 2*np.pi
                dphi[dphi<0] += 2*np.pi

                dR = 0.5*(dphi + dphi[right])
                dL = 0.5*(dphi + dphi[left])

                fR = self.prim[k][i][right,:]
                fC = self.prim[k][i][:,:]
                fL = self.prim[k][i][left,:]

                grad[k][i][:,:,1] = (dL*(fR-fC)/dR + dR*(fC-fL)/dL) / (dL+dR)

        # R-direction
        for k in xrange(self.nz_tot):
            dz = self.zFaces[k+1] - self.zFaces[k]
            for i,r in enumerate(self.rFaces):
                if i == 0 or i == self.nr_tot:
                    continue

                dr = 0.5 * (self.rFaces[i+1] - self.rFaces[i-1])
                r = self.rFaces[i]

                npL = self.np[k][i-1]
                npR = self.np[k][i]
                pFacesL = self.pFaces[k][i-1]
                pFacesR = self.pFaces[k][i-1]
                jL0 = np.argmin(pFacesL)
                jR0 = np.argmin(pFacesR)

                primL = self.prim[k][i-1]
                primR = self.prim[k][i]

                jL = jL0
                jR = jR0

                done = False

                while not done:

                    phiL1, phiL2, phiR1, phiR2 = fixPhi(pFacesL[jL-1], 
                                                        pFacesL[jL],
                                                        pFacesR[jR-1],
                                                        pFacesR[jR])
                    phiL = 0.5 * (phiL1 + phiL2)
                    phiR = 0.5 * (phiR1 + phiR2)

                    phi1 = max(phiL1, phiR1)
                    phi2 = min(phiL2, phiR2)
                    
                    phi = 0.5 * (phi1+phi2)
                    dphi = phi2-phi1

                    fL = primL[jL] + (phi-phiL) * grad[k][i-1][jL,:,1]
                    fR = primR[jR] + (phi-phiR) * grad[k][i  ][jR,:,1]

                    df = (fR - fL) / dr
                    dA = r * dphi * dz

                    grad[k][i-1][jL,:,0] += df*dA
                    grad[k][i  ][jR,:,0] += df*dA

                    if phiL2 < phiR2:
                        jL += 1
                        if jL == npL:
                            jL = 0
                    else:
                        jR += 1
                        if jR == npR:
                            jR = 0

                    if jL == jL0 and jR == jR0:
                        done = True

        for k in xrange(self.nz_tot):
            dz = self.zFaces[k+1] - self.zFaces[k]
            for i in xrange(self.nr_tot):
                r1 = self.rFaces[i]
                r2 = self.rFaces[i+1]
                for j in xrange(self.np[k,i]):
                    dphi = self.pFaces[k][i][j] - self.pFaces[k][i][j-1]
                    if dphi < 0:
                        dphi += 2*np.pi
                    dA = (r1+r2)*dphi*dz
                    grad[k][i][j,:,0] /= dA

        # Z - Direction


def _fixPhi(phiL1, phiL2, phiR1, phiR2):

    phiL = 0.5 * (phiL1 + phiL2)
    phiR = 0.5 * (phiR1 + phiR2)

    if phiL1 > phiL2 and phiR1 > phiR2:
        phiL1 -= 2*np.pi
        phiR1 -= 2*np.pi
    elif phiL1 > phiL2:
        if phiR < np.pi:
            phiL1 -= 2*np.pi
        else:
            phiL2 += 2*np.pi
    elif phiR1 > phiR2:
        if phiL < np.pi:
            phiR1 -= 2*np.pi
        else:
            phiR2 += 2*np.pi

    return phiL1, phiL2, phiR1, phiR2




def _r_func(r, r2, fac, sigma, r0):
    #print r, r2, fac, sigma, r0
    f = r - r2 + fac*sigma*math.sqrt(0.25*math.pi)*(math.erf((r-r0)/sigma)+1.0)
    #print f
    return f


if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Gimme some files.")
        sys.exit()

    if "archive" in sys.argv[1]:
        g = Grid(archive=sys.argv[1])

        print g.rFaces
        print g.zFaces
        print g.pFaces
        print g.prim

    else:
        if len(sys.argv) == 2:
            pars = rp.readParfile(sys.argv[1])
            g = Grid(pars)

            print g.rFaces
            print g.zFaces
            print g.pFaces

        elif len(sys.argv) > 2:
            
            pars = rp.readParfile(sys.argv[1])
            g = Grid(pars)
            g.loadChkpt(sys.argv[2])

    g.saveArchive("archive.h5")



