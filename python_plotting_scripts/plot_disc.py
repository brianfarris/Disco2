import math
from numpy import *
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import sys
import numpy as np
import readChkpt as rc

GAM = 1.3333333333
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
    plt.xlabel(r"$r$")
    plt.ylabel(ylabel)

def plot_r_profile(filename, sca='linear'):
    dat = rc.readChkpt(filename)
    t = dat[0]
    r = dat[1]
    rho = dat[4]
    T = dat[5]
    vr = dat[6]
    vp = dat[7]

    real = r>2*M
    r = r[real]
    rho = rho[real]
    T = T[real]
    vr = vr[real]
    vp = vp[real]
    
    i = r.argmax()
    R = np.linspace(r.min(), r.max(), 100)
    #R = np.logspace(np.log10(r.min()), np.log10(r.max()), 100)
    RHO = rho[i] * np.power(R/r[i], -1.5)
    TTT = T[i] * np.power(R/r[i], -1.0)
    URR = vr[i] * np.power(R/r[i], -0.5)
    UPP = vp[i] * np.power(R/r[i], -1.5)

    rhoh = rho*(1.0 + GAM/(GAM-1.0)*T)
    RHOH = RHO*(1.0 + GAM/(GAM-1.0)*TTT)

    #u0 = 1.0/np.sqrt(1.0-2*M/r-vr*vr/(1-2*M/r)-r*r*vp*vp)
    #U00 = 1.0/np.sqrt(1.0-2*M/R-URR*URR/(1-2*M/R)-R*R*UPP*UPP)
    u0 = 1.0/np.sqrt(1.0-2*M/r-4*M/r*vr-(1+2*M/r)*vr*vr-r*r*vp*vp)
    U00 = 1.0/np.sqrt(1.0-2*M/R-4*M/R*URR-(1+2*M/R)*URR*URR-R*R*UPP*UPP)

    #fac = 1.0/(1.0-2*M/(r-2*M)*vr)
    #vr = vr*fac
    #vp = vp*fac

    #V = u0*vr / np.sqrt(1-2*M/r+u0*u0*vr*vr)
    #VVV = U00*URR / np.sqrt(1-2*M/R+U00*U00*URR*URR)

    #mach = np.sqrt((vr*vr/(1-2*M/r)+r*r*vp*vp)/(GAM*P/rhoh))
    #MACH = np.sqrt((URR*URR/(1-2*M/R)+R*R*UPP*UPP)/(GAM*PPP/RHOH))

    P = rho*T
    PPP = RHO*TTT

    H = np.sqrt(P*r*r*r/(rhoh*M))/u0
    HHH = np.sqrt(PPP*R*R*R/(RHOH*M))/U00

    print("Plotting t = {0:g}".format(t))

    #Plot.
    fig = plt.figure(figsize=(12,9))

    plt.subplot(331)
    plot_r_profile_single(r, rho, sca, r"$\rho_0$", R, RHO)

    plt.subplot(332)
    plot_r_profile_single(r, P, sca, r"$P = \rho_0 T$", R, PPP)
    
    plt.subplot(333)
    plot_r_profile_single(r, T, sca, r"$T$", R, TTT)
    
    plt.subplot(334)
    plot_r_profile_single(r, vr, "linear", r"$v^r$", R, URR)
    plt.plot(R, (1-2*M/R)/(1+2*M/R), 'r--')
    plt.plot(R, -1*np.ones(R.shape), 'r--')
    plt.plot(R, (0.5-2*M/R)/(1+2*M/R), 'b--')
    plt.plot(R, (-0.5-2*M/R)/(1+2*M/R), 'b--')
    plt.gca().set_xscale("log")
    
    plt.subplot(335)
    plot_r_profile_single(r, vp, sca, r"$v^\phi$", R, UPP)

    #plt.subplot(336)
    #plot_r_profile_single(r, mach, sca, r"$\mathcal{M}$", R, MACH)

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

