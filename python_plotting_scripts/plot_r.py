import math
from numpy import *
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import sys
import numpy as np
import readChkpt as rc

GAM = 1.3333333333

def plot_r_profile_single(r, f, sca, ylabel, R=None, F=None):

    if sca == 'log' and (f < 0).all():
        plt.plot(r, -f, 'k+')
        if R != None and F != None:
            plt.plot(R, -F, 'r')
    else:
        plt.plot(r, f, 'k+')
        if R != None and F != None:
            plt.plot(R, F, 'r')
    
    plt.gca().set_xscale(sca)
    plt.gca().set_yscale(sca)
    plt.xlabel(r"$r$")
    plt.ylabel(ylabel)

def plot_r_profile(filename, sca='linear'):
    dat = rc.readChkpt(filename)
    t = dat[0]
    r = dat[1]
    rho = dat[4]
    P = dat[5]
    vr = dat[6]
    vp = dat[7]
    
    i = r.argmax()
    R = np.linspace(r.min(), r.max(), 1000)
    RHO = rho[i] * np.power(R/r[i], -0.5)
    PPP = P[i] * np.power(R/r[i], -1.5)
    URR = vr[i] * np.power(R/r[i], -0.5)
    UPP = vp[i] * np.power(R/r[i], -1.5)

    mach = np.sqrt((vr*vr+r*r*vp*vp)/(GAM*P/rho))
    MACH = np.sqrt((URR*URR+R*R*UPP*UPP)/(GAM*PPP/RHO))

    print("Plotting t = {0:g}".format(t))

    #Plot.
    fig = plt.figure(figsize=(12,9))
    plt.subplot(231)
    plot_r_profile_single(r, rho, sca, r"$\rho_0$", R, RHO)
    plt.subplot(232)
    plot_r_profile_single(r, P, sca, r"$P$", R, PPP)
    plt.subplot(234)
    plot_r_profile_single(r, vr, sca, r"$v^r$", R, URR)
    plt.plot(R, -np.ones(R.shape), 'b')
    plt.plot(R, -(2.0-R)/(2.0+R), 'b')
    plt.subplot(235)
    plot_r_profile_single(r, vp, sca, r"$v^\phi$", R, UPP)
    plt.subplot(236)
    plot_r_profile_single(r, mach, sca, r"$\mathfrak{M}$", R, MACH)
    plt.subplot(233)
    plot_r_profile_single(r, P/rho, sca, r"$T = P/\rho_0$", R, PPP/RHO)

    plt.tight_layout()

    outname = "plot_r_{0}.png".format("_".join( filename.split(".")[0].split("_")[1:] ))

    print("Saving {0:s}...".format(outname))
    plt.savefig(outname)

    return fig

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("\nGive me a checkpoint (.h5) file.\n")
        sys.exit()

    elif len(sys.argv) == 2:
        filename = sys.argv[1]
        fig = plot_r_profile(filename, sca='log')
        plt.show()

    else:
        for filename in sys.argv[1:]:
            fig = plot_r_profile(filename, sca='log')
            plt.close(fig)

