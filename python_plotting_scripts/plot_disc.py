import math
from numpy import *
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import sys
import numpy as np
import readChkpt as rc

GAM = 1.3333333333
M = 1.0e-6
scale = 'linear'

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

    real = r>2*M
    r = r[real]
    rho = rho[real]
    P = P[real]
    vr = vr[real]
    vp = vp[real]
    
    i = r.argmax()
    R = np.logspace(np.log10(r.min()), np.log10(r.max()), 100)
    RHO = rho[i] * np.power(R/r[i], -0.5)
    PPP = P[i] * np.power(R/r[i], -1.5)
    URR = vr[i] * np.power(R/r[i], -0.5)
    UPP = vp[i] * np.power(R/r[i], -1.5)

    rhoh = rho + GAM*P/(GAM-1.0)
    RHOH = RHO + GAM*PPP/(GAM-1.0)

    #u0 = 1.0/np.sqrt(1.0-2*M/r-vr*vr/(1-2*M/r)-r*r*vp*vp)
    #U00 = 1.0/np.sqrt(1.0-2*M/R-URR*URR/(1-2*M/R)-R*R*UPP*UPP)
    u0 = 1.0/np.sqrt(1.0-2*M/r-4*M/r*vr-(1+2*M/r)*vr*vr-r*r*vp*vp)
    U00 = 1.0/np.sqrt(1.0-2*M/R-4*M/R*URR-(1+2*M/R)*URR*URR-R*R*UPP*UPP)

    fac = 1.0/(1.0-2*M/(r-2*M)*vr)
    vr = vr*fac
    vp = vp*fac

    V = u0*vr / np.sqrt(1-2*M/r+u0*u0*vr*vr)
    VVV = U00*URR / np.sqrt(1-2*M/R+U00*U00*URR*URR)

    mach = np.sqrt((vr*vr/(1-2*M/r)+r*r*vp*vp)/(GAM*P/rhoh))
    MACH = np.sqrt((URR*URR/(1-2*M/R)+R*R*UPP*UPP)/(GAM*PPP/RHOH))

    T = P/rho
    TTT = PPP/RHO

    H = 2*np.sqrt(P*r*r*r/(rhoh*M))/u0
    HHH = 2*np.sqrt(PPP*R*R*R/(RHOH*M))/U00

    print("Plotting t = {0:g}".format(t))

    #Plot.
    fig = plt.figure(figsize=(12,9))
    plt.subplot(331)
    plot_r_profile_single(r, rho/H, sca, r"$\rho_0$", R, RHO/HHH)
    plt.subplot(332)
    plot_r_profile_single(r, P/H, sca, r"$P$", R, PPP/HHH)
    plt.subplot(333)
    plot_r_profile_single(r, T, sca, r"$T = P/\rho_0$", R, TTT)
    
    plt.subplot(334)
    plot_r_profile_single(r, V, sca, r"$V$", R, VVV)
#    plt.plot(R, -np.ones(R.shape), 'b')
#    plt.plot(R, -(2.0-R)/(2.0+R), 'b')
    plt.subplot(335)
    plot_r_profile_single(r, vp, sca, r"$v^\phi$", R, UPP)
    plt.subplot(336)
    plot_r_profile_single(r, mach, sca, r"$\mathcal{M}$", R, MACH)

    plt.subplot(337)
    plot_r_profile_single(r, H, sca, r"$\mathcal{H}$", R, HHH)

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
        fig = plot_r_profile(filename, sca=scale)
        plt.show()

    else:
        for filename in sys.argv[1:]:
            fig = plot_r_profile(filename, sca=scale)
            plt.close(fig)

