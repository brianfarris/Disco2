import sys
import numpy as np
import discoGrid as dg
import readParfile as rp

def remapGrid(g0, g1):
    # Maps g0's data onto g1, taking into account differing grid structures.

    if g0.nz_tot != 1 or g1.nz_tot != 1:
        _remapGrid2D(g0, g1)
    else:
        _remapGrid3D(g0, g1)

def _remapGrid2D(g0, g1):
    # Instantiate g1's prim with values linearly interpolated from g0, assuming
    # both grids only have a single 'z' zone.

    if g0.grad is None:
        g0.plm()

    prim = []
    nq = g0.nq

    # Initialize all prims to 0
    for i in xrange(g1.nr_tot):
        prim.append(np.zeros((g1.np[0,i],nq)))

    # Just do it. (TM)

    done = False
    i0 = 0
    i1 = 0

    # Set prim with integral of g0 grid over volume
    while not done:
        print("Interpolating {0:d} and {1:d}".format(i0, i1))
        r0L = g0.rFaces[i0]
        r0R = g0.rFaces[i0+1]
        r1L = g1.rFaces[i1]
        r1R = g1.rFaces[i1+1]

        rL = max(r0L, r1L)
        rR = min(r0R, r1R)

        if rR > rL:

            r0 = 0.5 * (r0L + r0R)
            r = 0.5*(rL + rR)
            dr = rR - rL

            np0 = g0.np[0][i0]
            np1 = g1.np[0][i1]
            pFaces0 = g0.pFaces[0][i0]
            pFaces1 = g1.pFaces[0][i1]

            prim0 = g0.prim[0][i0]
            grad0 = g0.grad[0][i0]

            j0_0 = np.argmin(pFaces0)
            j1_0 = np.argmin(pFaces1)

            j0 = j0_0
            j1 = j1_0

            donePhi = False

            while not donePhi:

                phi0L, phi0R, phi1L, phi1R = dg.fixPhi(pFaces0[j0-1],
                                                        pFaces0[j0],
                                                        pFaces1[j1-1],
                                                        pFaces1[j1])
                phiL = max(phi0L, phi1L)
                phiR = min(phi0R, phi1R)

                phi0 = 0.5*(phi0L + phi0R)
                phi = 0.5*(phiL + phiR)
                dphi = phiR-phiL

                if dphi < 0:
                    print("AAAAAAHHHHHHH!")

                q = prim0[j0,:] + (r-r0)*grad0[j0,:,0] \
                                + (phi-phi0)*grad0[j0,:,1]
                dA = r * dr * dphi

                prim[i1][j1,:] += q*dA

                if phi0R < phi1R:
                    j0 += 1
                    if j0 == np0:
                        j0 = 0
                else:
                    j1 += 1
                    if j1 == np1:
                        j1 = 0

                if j0 == j0_0 and j1 == j1_0:
                    donePhi = True

        if r0R < r1R:
            i0 += 1
        else:
            i1 += 1

        if i1 == g1.nr_tot:
            done = True

    # Calculate densities
    for i in xrange(g1.nr_tot):
        r1 = g1.rFaces[i]
        r2 = g1.rFaces[i+1]
        r = 0.5*(r1+r2)
        dr = r2-r1
        for j in xrange(g1.np[0,i]):
            phi1 = g1.pFaces[0][i][j-1]
            phi2 = g1.pFaces[0][i][j]
            dphi = phi2-phi1
            if dphi < 0:
                dphi += 2*np.pi

            prim[i][j,:] /= r * dr * dphi

    # Done!
    g1.prim = [prim]
    g1.nq = nq
    g1.T = g0.T

def _remapGrid3D(g0, g1):
    # Instantiate g1's prim with values linearly interpolated from g0, assuming
    # both grids only have a single 'z' zone.

    if g0.grad is None:
        g0.plm()

    prim = []
    nq = g0.nq

    # Initialize all prims to 0
    for k in xrange(g1.nz_tot):
        slice = []
        for i in xrange(g1.nr_tot):
            slice.append(np.zeros((g1.np[k,i],nq)))
        prim.append(slice)

    # Just do it. (TM)

    # Set prim with integral of g0 grid over volume
    doneZ = False
    k0 = 0
    k1 = 0

    while not doneZ:
        z0L = g0.zFaces[k0]
        z0R = g0.zFaces[k0+1]
        z1L = g1.zFaces[k1]
        z1R = g1.zFaces[k1+1]

        zL = max(z0L, z1L)
        zR = min(z0R, z1R)

        if zR > zL:

            z0 = 0.5*(z0L+z0R)
            z = 0.5*(zL+zR)
            dz = zR - zL

            doneR = False
            i0 = 0
            i1 = 0

            while not doneR:
                print("Interpolating {0:d} and {1:d}".format(i0, i1))
                r0L = g0.rFaces[i0]
                r0R = g0.rFaces[i0+1]
                r1L = g1.rFaces[i1]
                r1R = g1.rFaces[i1+1]

                rL = max(r0L, r1L)
                rR = min(r0R, r1R)

                if rR > rL:

                    r0 = 0.5 * (r0L + r0R)
                    r = 0.5*(rL + rR)
                    dr = rR - rL

                    np0 = g0.np[k0][i0]
                    np1 = g1.np[k1][i1]
                    pFaces0 = g0.pFaces[k0][i0]
                    pFaces1 = g1.pFaces[k1][i1]

                    prim0 = g0.prim[k0][i0]
                    grad0 = g0.grad[k1][i0]

                    j0_0 = np.argmin(pFaces0)
                    j1_0 = np.argmin(pFaces1)

                    j0 = j0_0
                    j1 = j1_0

                    donePhi = False

                    while not donePhi:

                        phi0L, phi0R, phi1L, phi1R = dg.fixPhi(pFaces0[j0-1],
                                                                pFaces0[j0],
                                                                pFaces1[j1-1],
                                                                pFaces1[j1])
                        phiL = max(phi0L, phi1L)
                        phiR = min(phi0R, phi1R)

                        phi0 = 0.5*(phi0L + phi0R)
                        phi = 0.5*(phiL + phiR)
                        dphi = phiR-phiL

                        if dphi < 0:
                            print("AAAAAAHHHHHHH!")

                        q = prim0[j0,:] + (r-r0)*grad0[j0,:,0] \
                                        + (phi-phi0)*grad0[j0,:,1] \
                                        + (z-z0)*grad0[j0,:,2]
                        dV = r * dr * dphi * dz

                        prim[k1][i1][j1,:] += q*dV

                        if phi0R < phi1R:
                            j0 += 1
                            if j0 == np0:
                                j0 = 0
                        else:
                            j1 += 1
                            if j1 == np1:
                                j1 = 0

                        if j0 == j0_0 and j1 == j1_0:
                            donePhi = True

                if r0R < r1R:
                    i0 += 1
                else:
                    i1 += 1

                if i1 == g1.nr_tot:
                    doneR = True

        if z0R < z1R:
            k0 += 1
        else:
            k1 += 1

        if k1 == g1.nz_tot:
            doneZ = True

    # Calculate densities
    for k in xrange(g1.nz_tot):
        z1 = g1.rFaces[k]
        z2 = g1.rFaces[k+1]
        z = 0.5*(z1+z2)
        dz = z2-z1
        for i in xrange(g1.nr_tot):
            r1 = g1.rFaces[i]
            r2 = g1.rFaces[i+1]
            r = 0.5*(r1+r2)
            dr = r2-r1
            for j in xrange(g1.np[k,i]):
                phi1 = g1.pFaces[k][i][j-1]
                phi2 = g1.pFaces[k][i][j]
                dphi = phi2-phi1
                if dphi < 0:
                    dphi += 2*np.pi

                prim[k][i][j,:] /= r * dr * dphi * dz

    # Done!
    g1.prim = prim
    g1.nq = nq
    g1.T = g0.T


if __name__ == "__main__":

    if len(sys.argv) != 3 and len(sys.argv) != 4:
        print("Remap a Disco2 archive onto a new parfile specified grid.")
        print("usage: python remapArchive.py <archive> <parfile> [<num_processes>]")
        sys.exit()

    archiveName = sys.argv[1]
    parfileName = sys.argv[2]

    if len(sys.argv) == 4:
        numProc = int(sys.argv[3])
    else:
        numProc = 1

    pars = rp.readParfile(parfileName)

    root = ".".join(parfileName.split(".")[:-1])

    print("Making grids...")
    g0 = dg.Grid(archive=archiveName)
    g1 = dg.Grid(pars=pars)

    if g0.nz_tot == 1 and g1.nz_tot == 1:
        print("Remapping...")
        remapGrid(g0, g1)

    print("Saving new grid...")
    g1.saveArchive(root + "_archive" + ".h5")
    g1.saveCheckpoint(root + "_checkpoint" + ".h5", numProc)

