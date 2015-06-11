import math
from numpy import *
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import sys
import numpy as np
import discopy as dp

GAM = 1.3333333333
M = 1.0
xsca = "log"
ysca = "log"

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
    

def plot_r_profile(filename, xsca, ysca='linear'):
    dat = dp.readCheckpoint(filename)
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
    plot_r_profile_single(r, rho, xsca, ysca, r"$\rho_0$")
    plt.subplot(222)
    plot_r_profile_single(r, P, xsca, ysca, r"$P$")
    plt.subplot(223)
    plot_r_profile_single(r, vr, xsca, "linear", r"$v^r$")
    plt.subplot(224)
    plot_r_profile_single(r, vp, xsca, "linear", r"$v^\phi$")
    plt.plot(R, np.sqrt(M/(R*R*R)), 'r')

    plt.title(r"t = {0:.3e}".format(t))
    plt.tight_layout()

    outname = "plot_r_{0}.png".format("_".join( filename.split("/")[-1].split(".")[0].split("_")[1:] ))

    print("Saving {0:s}...".format(outname))
    plt.savefig(outname)

    return fig

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("\nGive me a checkpoint (.h5) file.\n")
        sys.exit()

    elif len(sys.argv) == 2:
        filename = sys.argv[1]
        fig = plot_r_profile(filename, xsca=xsca, ysca=ysca)
        plt.show()

    else:
        for filename in sys.argv[1:]:
            fig = plot_r_profile(filename, xsca=xsca, ysca=ysca)
            plt.close(fig)

