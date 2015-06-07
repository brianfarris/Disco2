import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as pat
import matplotlib.collections as coll
import matplotlib.cm as cm
import remapChkpt as rc
import readParfile as rp

def plotGrid(g, ax, k=0):

    PHI = np.linspace(0, 2*math.pi, 1000)

    for i in xrange(g.nr_tot):
        r1 = g.rFaces[i]
        r2 = g.rFaces[i+1]
        ax.plot(r1*np.cos(PHI), r1*np.sin(PHI), 'k')
        r = np.linspace(r1, r2, 2)
        for phi in g.pFaces[k][i]:
            ax.plot(r*math.cos(phi), r*math.sin(phi), 'k')
    ax.plot(g.rFaces[-1]*np.cos(PHI), g.rFaces[-1]*np.sin(PHI), 'k')


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

    pars = rp.readParfile(sys.argv[1])
    print("Making grid...")
    g = rc.Grid(pars)


    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    if len(sys.argv) == 2:
        print("Plotting grid...")
        plotGrid(g, ax)
        sizePlot(g, ax)

    elif len(sys.argv) > 2:
        print("Loading Checkpoint...")
        g.loadChkpt(sys.argv[2])
        print("Plotting grid...")
        #plotGrid(g, ax)
        print("Plotting data...")
        plotGridData(g, ax, q=0)
        sizePlot(g, ax)

    plt.show()

