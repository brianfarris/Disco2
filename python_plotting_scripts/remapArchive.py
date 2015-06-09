import sys
import numpy as np
import discoGrid as dg
import readParfile as rp

def remapGrid2D(g0, g1):
    # Instantiate g1's prim with values linearly interpolated from g0.

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
        remapGrid2D(g0, g1)

    print("Saving new grid...")
    g1.saveArchive(root + "_archive" + ".h5")
    g1.saveChkpt(root + "_checkpoint" + ".h5", numProc)

