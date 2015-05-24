import math
from numpy import *
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import sys
import numpy as np
import readChkpt as rc

GAM = 5.0/3.0
M = 1.0
scale = 'log'

def plot_r_profile_single(r, f, sca, ylabel, R=None, F=None):

    if sca == 'log' and (f < 0).any():
        plt.plot(r, -f, 'k+')
        if R != None and F != None:
            plt.plot(R, -F, 'r')
    else:
        plt.plot(r, f, 'k+')
        if R != None and F != None:
            plt.plot(R, F, 'r')
    
    plt.gca().set_xscale(sca)
    if (f == 0.0).all():
        plt.gca().set_yscale('linear')
    else:
        plt.gca().set_yscale(sca)
   
    xmin = math.pow(10, math.floor(math.log10(r.min())))
    plt.axvspan(xmin, 2*M, color='lightgrey', alpha=0.5)
    plt.xlim(xmin, r.max())
    
    plt.xlabel(r"$r$")
    plt.ylabel(ylabel)

def plot_r_profile(filename, sca='linear'):
    dat = rc.readChkpt(filename)
    t = dat[0]
    r = dat[1]
    rho = dat[4]
    p = dat[5]
    vr = dat[6]
    vp = dat[7]
    q = dat[10]

    #real = r>2*M
    real = r>-1.0
    r = r[real]
    rho = rho[real]
    p = p[real]
    vr = vr[real]
    vp = vp[real]

    R = np.linspace(r.min(), r.max(), 500)
    
    rhoh = rho + GAM/(GAM-1.0)*p

    #u0 = 1.0/np.sqrt(1.0-2*M/r-vr*vr/(1-2*M/r)-r*r*vp*vp)
    #VR = vr/(1-2*M/r)
    #VP = r*vp/np.sqrt(1-2*M/r)
    #W = u0 * np.sqrt(1-2*M/r)
    u0 = 1.0/np.sqrt(1.0-2*M/r-4*M/r*vr-(1+2*M/r)*vr*vr-r*r*vp*vp)
    VR = (1+2*M/r)*vr + 2*M/r
    VP = np.sqrt(1+2*M/r)*r*vp
    W = u0/np.sqrt(1+2*M/r)

    T = p/rho
    H = np.sqrt(p*r*r*r/(rhoh*M))/u0

    cs = np.sqrt(GAM * p / rho)
    V = np.sqrt(VR*VR+VP*VP)
    mach = V/cs

    Mdot = -2*math.pi* r * rho * u0*vr

    print("Plotting t = {0:g}".format(t))

    #Plot.
    fig = plt.figure(figsize=(12,9))

    plt.subplot(331)
    plot_r_profile_single(r, rho, sca, r"$\Sigma$")
    #ax1 = plt.gca()
    #ax2 = ax1.twinx()
    #ax2.set_ylabel(r"$q$")
    #ax2.plot(r, q, 'r+')

    plt.subplot(332)
    plot_r_profile_single(r, p, sca, r"$\Pi$")
    
    plt.subplot(333)
    plot_r_profile_single(r, T, sca, r"$T = P/\rho_0$")
    
    plt.subplot(334)
    plot_r_profile_single(r, vr, "linear", r"$v^r$")
    plt.plot(R, (1-2*M/R)/(1+2*M/R), 'r--')
    plt.plot(R, -1*np.ones(R.shape), 'r--')
    plt.plot(R, (0.5-2*M/R)/(1+2*M/R), 'b--')
    plt.plot(R, (-0.5-2*M/R)/(1+2*M/R), 'b--')
    plt.plot(R, (-2*M/R)/(1+2*M/R), ls='--', color='lightgrey')
    plt.gca().set_xscale(sca)
    
    plt.subplot(335)
    plot_r_profile_single(r, vp, sca, r"$v^\phi$")
    plt.plot(R, 1.0/(R*np.sqrt(1+2*M/R)), 'r--')
    plt.plot(R, 0.5/(R*np.sqrt(1+2*M/R)), 'b--')

    plt.subplot(336)
    plot_r_profile_single(r, mach, "log", r"$\mathcal{M}$")
    plt.plot(r, W, 'b+')

    plt.subplot(337)
    plot_r_profile_single(r, H, sca, r"$\mathcal{H}$")

    plt.subplot(338)
    plot_r_profile_single(r, Mdot, "linear", r"$\dot{M}$")

    plt.subplot(339)
    plt.plot(r, VR, 'r+')
    plt.plot(r, VP, 'g+')
    plt.plot(r, V, 'b+')
    plot_r_profile_single(r, cs, "linear", r"$c_s$")

    plt.suptitle(r"$T = {0:.3g}$".format(t))

    plt.tight_layout()

    outname = "plot_idisc_{0}.png".format("_".join( filename.split("/")[-1].split(".")[0].split("_")[1:] ))

    print("Saving {0:s}...".format(outname))
    plt.savefig(outname)

    return fig

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("\nGive me a checkpoint (.h5) file.\n")
        sys.exit()

    elif len(sys.argv) == 2:
        filename = sys.argv[1]
        fig = plot_r_profile(filename, sca=scale)
        plt.show()

    else:
        for filename in sys.argv[1:]:
            fig = plot_r_profile(filename, sca=scale)
            plt.close(fig)

