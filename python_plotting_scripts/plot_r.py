import math
from numpy import *
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import sys
import numpy as np
import readChkpt as rc

GAM = 1.3333333333
M = 1.0

def plot_r_profile_single(r, f, xsca, ysca, ylabel, R=None, F=None):

    if ysca == 'log' and (f < 0).all():
        plt.plot(r, -f, 'k+')
        if R != None and F != None:
            plt.plot(R, -F, 'r')
    else:
        plt.plot(r, f, 'k+')
        if R != None and F != None:
            plt.plot(R, F, 'r')

    plt.gca().set_xscale(xsca)
    plt.gca().set_yscale(ysca)
    plt.xlabel(r"$r$")
    plt.ylabel(ylabel)
    
    if 2*M > r.min():
        plt.axvspan(plt.xlim()[0], 2*M, color='lightgrey', alpha=0.5)
    

def plot_r_profile(filename, sca='linear'):
    dat = rc.readChkpt(filename)
    t = dat[0]
    r = dat[1]
    rho = dat[4]
    P = dat[5]
    vr = dat[6]
    vp = dat[7]

    R = np.linspace(r.min(), r.max(), 200)

    print("Plotting t = {0:g}".format(t))

    #Plot.
    fig = plt.figure(figsize=(12,9))
    plt.subplot(221)
    plot_r_profile_single(r, rho, sca, sca, r"$\rho_0$")
    plt.subplot(222)
    plot_r_profile_single(r, P, sca, sca, r"$P$")
    plt.subplot(223)
    plot_r_profile_single(r, vr, sca, "linear", r"$v^r$")
    plt.subplot(224)
    plot_r_profile_single(r, vp, sca, "linear", r"$v^\phi$")
    plt.plot(R, np.sqrt(M/(R*R*R)), 'r')

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

