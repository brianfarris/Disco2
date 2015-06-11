import math
from numpy import *
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib.collections as coll
import matplotlib.colors as clrs
import matplotlib.patches as pat
import matplotlib.ticker as tkr
import sys
import numpy as np
import discopy as dp

M = 1.0
a = 0.0
gridscale = 'linear'
datscale = 'log'
A = a*M
poscmap = plt.cm.afmhot
divcmap = plt.cm.RdBu

BINSEP = 200.0
ROCHE = True
XROCHE = None
YROCHE = None
ZROCHE = None

def plot_equat_single(fig, ax, mesh, dat, gridscale="linear", gridbounds=None,
                    datscale="linear", datbounds=None, label=None, 
                    normal=False, **kwargs):
    N = 400

    #Set data bounds & scale.
    if datbounds is None:
        datbounds = np.array([dat.min(), dat.max()])

    if datscale == "log":
        v = np.logspace(math.floor(math.log10(datbounds[0])), math.ceil(math.log10(datbounds[1])), 
                            base=10.0, num=N)
        norm = clrs.LogNorm()
        locator = tkr.LogLocator()
        formatter = tkr.LogFormatter()
    else:
        if normal:
            datm = max(abs(datbounds[0]), abs(datbounds[1]))
            v = np.linspace(-datm, datm, N)
        else:
            v = np.linspace(datbounds[0], datbounds[1], N)
        norm = clrs.Normalize()
        locator = tkr.AutoLocator()
        formatter = tkr.ScalarFormatter()

    #Plot Disc Image
    C = ax.tricontourf(mesh, dat, v, norm=norm, **kwargs)
    colorbar = fig.colorbar(C, format=formatter, ticks=locator)

    #Plot Roche Equipotential
    if ROCHE:
        plot_roche(ax, mesh, gridbounds)

    #Patches to highlight horizona and ergosphere
    if M > 0.0:
        ergo = pat.Circle((0,0), 2*M)
        horizon = pat.Wedge((0,0), M*(1.0+math.sqrt(1.0-a*a)), 0, 360, 
                            2*M*math.sqrt(1.0-a*a))
        patches = coll.PatchCollection([ergo,horizon], cmap=plt.cm.Greys,alpha=0.4)
        colors = np.array([0.1,0.3])
        patches.set_array(colors)
        ax.add_collection(patches)

    #Formatting
    ax.set_aspect('equal')

    ax.set_xscale(gridscale)
    ax.set_yscale(gridscale)
    if gridbounds is not None:
        ax.set_xlim(gridbounds[0])
        ax.set_ylim(gridbounds[1])

    #Labels
    if label is not None:
        colorbar.set_label(label)
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$y$')

def plot_roche(ax, mesh, gridbounds):

    N = 256

    global XROCHE
    global YROCHE
    global ZROCHE

    if XROCHE == None or YROCHE is None or ZROCHE is None:
        print("Building Roche Lobe")
        if gridbounds is None:
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
        else:
            xlim = gridbounds[0]
            ylim = gridbounds[1]
        x = np.linspace(xlim[0], xlim[1], N)
        y = np.linspace(ylim[0], ylim[1], N)

        XROCHE,YROCHE = np.meshgrid(x, y)
        ZROCHE = -M/np.sqrt((XROCHE)*(XROCHE)+YROCHE*YROCHE) \
                -M/np.sqrt((XROCHE+BINSEP)*(XROCHE+BINSEP)+YROCHE*YROCHE) \
                -M*((XROCHE+0.5*BINSEP)*(XROCHE+0.5*BINSEP)+YROCHE*YROCHE) \
                    / (BINSEP*BINSEP*BINSEP)

    lvl = -4*M/BINSEP
    ax.contour(XROCHE, YROCHE, ZROCHE, levels=[lvl], colors='m')
    #ax.contour(XROCHE, YROCHE, ZROCHE, levels=[1.4*lvl,1.35*lvl,1.3*lvl,1.25*lvl,1.2*lvl,1.15*lvl,1.1*lvl,1.05*lvl,lvl,0.95*lvl,0.9*lvl,0.85*lvl,0.8*lvl,0.75*lvl,0.7*lvl,0.65*lvl,0.6*lvl], colors='m')
    #ax.contour(XROCHE, YROCHE, ZROCHE)


def make_plot(mesh, dat, gridscale="linear", gridbounds=None, 
                datscale="linear", datbounds=None, label=None, 
                title=None, filename=None, normal=False, **kwargs):

    fig = plt.figure(figsize=(12,9))
    ax = fig.add_subplot(1,1,1)
    plot_equat_single(fig, ax, mesh, dat, gridscale, gridbounds,
                    datscale, datbounds, label, normal, **kwargs)
    if title is not None:
        ax.set_title(title)
    if filename is not None:
        print("Saving {0:s}...".format(filename))
        fig.savefig(filename)
    plt.close()

def plot_all(filename, gridscale='linear', plot=True, bounds=None):

    print("Reading {0:s}".format(filename))

    dat = dp.readCheckpoint(filename)
    t = dat[0]
    r = dat[1]
    phi = dat[2]
    z = dat[3]
    rho = dat[4]
    T = dat[5]
    vr = dat[6]
    vp = dat[7]
    dV = dat[9]
    q = dat[10]

    inds = z==z[len(z)/2]
    r = r[inds]
    phi = phi[inds]
    z = z[inds]
    rho = rho[inds]
    T = T[inds]
    vr = vr[inds]
    vp = vp[inds]
    dV = dV[inds]
    q = q[inds]

    if bounds is None:
        bounds = []
        bounds.append([rho[rho==rho].min(), rho[rho==rho].max()])
        bounds.append([T[T==T].min(), T[T==T].max()])
        bounds.append([vr[vr==vr].min(), vr[vr==vr].max()])
        bounds.append([vp[vp==vp].min(), vp[vp==vp].max()])
        bounds.append([q[q==q].min(), q[q==q].max()])
        bounds = np.array(bounds)

    if plot:

        print("Plotting t = {0:g}".format(t))

        x = r*cos(phi)
        y = r*sin(phi)
        mesh = tri.Triangulation(x, y)

        gridbounds = np.array([[-r.max(),r.max()],[-r.max(),r.max()]])

        outpath = filename.split("/")[:-1]
        chckname = filename.split("/")[-1]

        title = "t = {0:.3g}".format(t)

        rhoname = "plot_disc_equat_{0}_{1}.png".format(
                    "_".join(chckname.split(".")[0].split("_")[1:]), "rho")
        Tname = "plot_disc_equat_{0}_{1}.png".format(
                    "_".join(chckname.split(".")[0].split("_")[1:]), "T")
        vrname = "plot_disc_equat_{0}_{1}.png".format(
                    "_".join(chckname.split(".")[0].split("_")[1:]), "vr")
        vpname = "plot_disc_equat_{0}_{1}.png".format(
                    "_".join(chckname.split(".")[0].split("_")[1:]), "vp")
        qname = "plot_disc_equat_{0}_{1}.png".format(
                    "_".join(chckname.split(".")[0].split("_")[1:]), "q")
        #Plot.

        #Rho
        make_plot(mesh, rho, gridscale="linear", gridbounds=gridbounds,
                datscale=datscale, datbounds=bounds[0],
                label=r'$\rho_0$', title=title, filename=rhoname, cmap=poscmap)

        #T
        make_plot(mesh, T, gridscale="linear", gridbounds=gridbounds,
                datscale="log", datbounds=bounds[1],
                label=r'$T$', title=title, filename=Tname, cmap=poscmap)

        
        #Vr
        make_plot(mesh, vr, gridscale="linear", gridbounds=gridbounds,
                datscale="linear", datbounds=bounds[2],
                label=r'$v^r$', title=title, filename=vrname, normal=True,
                cmap=divcmap)

        #Vp
        make_plot(mesh, vp, gridscale="linear", gridbounds=gridbounds,
                datscale="linear", datbounds=bounds[3],
                label=r'$v^\phi$', title=title, filename=vpname, normal=True,
                cmap=divcmap)
        
        #q
        make_plot(mesh, q, gridscale="linear", gridbounds=gridbounds,
                datscale="linear", datbounds=bounds[4],
                label=r'$q$', title=title, filename=qname, cmap=poscmap)

    return bounds


if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("\nGive me a checkpoint (.h5) file.\n")
        sys.exit()

    elif len(sys.argv) == 2:
        filename = sys.argv[1]
        plot_all(filename, gridscale=gridscale)
        plt.show()

    else:
        bounds = np.zeros((5,2))
        bounds[:,0] = np.inf
        bounds[:,1] = -np.inf
        for filename in sys.argv[1:]:
            b = plot_all(filename, gridscale=gridscale, plot=False)
            lower = b[:,0]<bounds[:,0]
            upper = b[:,1]>bounds[:,1]
            bounds[lower,0] = b[lower,0]
            bounds[upper,1] = b[upper,1]

        for filename in sys.argv[1:]:
            plot_all(filename, gridscale=gridscale, plot=True, bounds=bounds)

