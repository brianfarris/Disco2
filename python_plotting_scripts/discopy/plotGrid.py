import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as pat
import matplotlib.collections as coll
import matplotlib.cm as cm
import discoGrid as dg
import readParfile as rp

def plotGrid(g, ax, k=0, N=1000):

    PHI = np.linspace(0, 2*math.pi, N)

    lines = []
    
    phiLines = np.zeros((g.nr_tot+1, N, 2))
    phiLines[:,:,0] = g.rFaces[:,None] * np.cos(PHI)[None,:]
    phiLines[:,:,1] = g.rFaces[:,None] * np.sin(PHI)[None,:]

    rLines = np.zeros((g.np[k,:].sum(),2,2))

    j = 0
    for i in xrange(g.nr_tot):
        r1 = g.rFaces[i]
        r2 = g.rFaces[i+1]

        for phi in g.pFaces[k][i]:
            cp = math.cos(phi)
            sp = math.sin(phi)
            rLines[j,0,0] = r1*cp
            rLines[j,0,1] = r1*sp
            rLines[j,1,0] = r2*cp
            rLines[j,1,1] = r2*sp
            j += 1

    rColl = coll.LineCollection(rLines, colors=[(0,0,0,1)], linewidths=0.5)
    phiColl = coll.LineCollection(phiLines, colors=[(0,0,0,1)], linewidths=0.5)
    ax.add_collection(rColl)
    ax.add_collection(phiColl)


def plotGridData(g, ax, k=0, q=0, cmap=cm.spring):
    
    patches = []
    vals = np.zeros(g.np[k].sum(), dtype=np.float64)

    d2r = 180/math.pi

    v = 0
    for i in xrange(g.nr_tot):
        r1 = g.rFaces[i]
        r2 = g.rFaces[i+1]
        dr = r2-r1

        vals[v : v+g.np[k,i]] = g.prim[k][i][:,q]
        v += g.np[k,i]

        for j in xrange(g.np[k,i]):
            p1 = d2r * g.pFaces[k][i][j-1]
            p2 = d2r * g.pFaces[k][i][j]
            wedge = pat.Wedge((0,0), r2, p1, p2, width=dr)
            patches.append(wedge)

    p = coll.PatchCollection(patches, cmap=cmap, linewidths=0)
    p.set_array(vals)
    ax.add_collection(p)
    plt.colorbar(p)

def sizePlot(g, ax, rmax=None):

    rlim = g.rFaces[-1]

    if rmax is not None:
        rlim = rmax

    ax.set_xlim(-rlim, rlim)
    ax.set_ylim(-rlim, rlim)
    ax.set_aspect("equal")

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Need a parfile")
        sys.exit()

    if "archive" in sys.argv[1]:
        g = dg.Grid(archive=sys.argv[1])
        

#        plotGrid(g, ax)

        for q in xrange(g.nq):
            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)
            plotGridData(g, ax, q=q, cmap=cm.afmhot)
            sizePlot(g, ax)
            fig.savefig("grid_{0:d}.png".format(q))
            fig.savefig("grid_{0:d}.pdf".format(q))
            plt.close()

    else:
        pars = rp.readParfile(sys.argv[1])
        print("Making grid...")
        g = dg.Grid(pars)

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

        if len(sys.argv) == 2:
            print("Plotting grid...")
            plotGrid(g, ax)
            sizePlot(g, ax)

        elif len(sys.argv) > 2:
            print("Loading Checkpoint...")
            g.loadChkpt(sys.argv[2])
            print("Plotting data...")
            plotGridData(g, ax, q=0)
            sizePlot(g, ax)

    plt.show()

