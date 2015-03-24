import math
from numpy import *
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import sys
import numpy as np
import readChkpt as rc

GAM = 5.0/3.0
M = 1.0
a = 0.0
scale = 'log'
eos_x1 = 1.0
eos_x2 = 0.0
eos_x3 = 0.0
alpha = 0.0
A = a*M

c = 2.99792458e10
G = 6.6738e-8
h = 6.62606957e-27
kb = 1.3806488e-16
sb = 1.56055371e59
mp = 1.672621777e-24
rg_solar = 1.4766250385e5
M_solar = 1.9884e33

def P_gas(rho, T):
    return rho * T

def P_rad(rho, T):
    return 4.0*sb/(3.0*c) * (T*mp*c*c)**4 / (c*c)

def P_deg(rho, T):
    return 2*np.pi*h*c/3.0 * np.power(3*rho/(8*np.pi*mp),4.0/3.0) / (c*c)

def e_gas(rho, T):
    return T / (GAM-1.0)

def e_rad(rho, T):
    return 4.0*sb * (T*mp*c*c)**4 / (c*rho*c*c)

def e_deg(rho, T):
    return h*c/(4.0*mp) * np.power(3*rho/(8*np.pi*mp),1.0/3.0) / (c*c)

def dPdr(rho, T):
    return T + h*c/(3*mp)*np.power(3*rho/(8*np.pi*mp),1.0/3.0) / (c*c)

def dPdT(rho, T):
    return rho + 16.0*sb/(3.0*c) * (T*mp*c*c)**3 * mp

def dedr(rho, T):
    return -4.0*sb * (T*mp*c*c)**4 / (c*rho*rho*c*c) + 3*h*c/(32*np.pi*mp*mp) * np.power(3*rho/(8*np.pi*mp),-2.0/3.0) / (c*c)

def dedT(rho, T):
    return c*c/(GAM-1.0) + 16.0*sb*(T*mp*c*c)**3 * mp / (c*rho)

def cs2(rho, T):
    P = P_gas(rho,T)+P_rad(rho,T)+P_deg(rho,T)
    e = e_gas(rho,T)+e_rad(rho,T)+P_deg(rho,T)
    enth = 1.0 + e + P/rho
    return (dPdr(rho,T)*dedT(rho,T) - dedr(rho,T)*dPdT(rho,T)
            + P*dPdT(rho,T)/(rho*rho)) / (dedT(rho,T)*enth)

def plot_r_profile_single(r, f, sca, ylabel, bounds=None, R=None, F=None):

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

    plt.axvspan(M*(1.0-math.sqrt(1.0-a*a)), M*(1.0+math.sqrt(1.0-a*a)), 
                    color='grey', alpha=0.5)
    plt.axvspan(plt.xlim()[0], 2*M, color='lightgrey', alpha=0.5)
    plt.xlim(r.min(), r.max())

    #if bounds is not None:
    #    if sca == "log":
    #        upper = 10.0 ** (math.ceil(math.log10(bounds[1]))+0.1)
    #        if bounds[0] > 0:
    #            lower = 10.0 ** (math.floor(math.log10(bounds[0]))-0.1)
    #        else:
    #            lower = 10.0 ** (math.floor(math.log10(f[f>0].min()))-0.1)
    #        plt.ylim(lower, upper)

def plot_r_profile(filename, sca='linear', plot=True, bounds=None):

    print("Reading {0:s}".format(filename))

    dat = rc.readChkpt(filename)
    t = dat[0]
    r = dat[1]
    rho = dat[4]
    T = dat[5]
    vr = dat[6]
    vp = dat[7]

    #real = r>2*M
    real = r > -1
    r = r[real]
    rho = rho[real]
    T = T[real]
    vr = vr[real]
    vp = vp[real]

    inds = np.argsort(r)
    r = r[inds]
    rho = rho[inds]
    T = T[inds]
    vr = vr[inds]
    vp = vp[inds]

    i = r.argmax()
    #R = np.linspace(r.min(), r.max(), 100)
    R = np.logspace(np.log10(r.min()), np.log10(r.max()), 100)
    RHO = rho[i] * np.power(R/r[i], -1.5)
    TTT = T[i] * np.power(R/r[i], -1.0)
    URR = vr[i] * np.power(R/r[i], -0.5)
    UPP = vp[i] * np.power(R/r[i], -1.5)

    #u0 = 1.0/np.sqrt(1.0-2*M/r-vr*vr/(1-2*M/r)-r*r*vp*vp)
    #U00 = 1.0/np.sqrt(1.0-2*M/R-URR*URR/(1-2*M/R)-R*R*UPP*UPP)
    u0 = 1.0/np.sqrt(1.0-2*M/r - 4*M/r*vr + 4*M*A/r*vp - (1+2*M/r)*vr*vr
                    + 2*A*(1+2*M/r)*vr*vp - (r*r+A*A*(1+2*M/r))*vp*vp)
    #U00 = 1.0/np.sqrt(1.0-2*M/R-4*M/R*URR-(1+2*M/R)*URR*URR-R*R*UPP*UPP)

    #fac = 1.0/(1.0-2*M/(r-2*M)*vr)
    #vr = vr*fac
    #vp = vp*fac

    #V = u0*vr / np.sqrt(1-2*M/r+u0*u0*vr*vr)
    #VVV = U00*URR / np.sqrt(1-2*M/R+U00*U00*URR*URR)

    #mach = np.sqrt((vr*vr/(1-2*M/r)+r*r*vp*vp)/(GAM*P/rhoh))
    #MACH = np.sqrt((URR*URR/(1-2*M/R)+R*R*UPP*UPP)/(GAM*PPP/RHOH))

    Pg = eos_x1*P_gas(rho, T)
    Pr = eos_x2*P_rad(rho, T)
    Pd = eos_x3*P_deg(rho, T)
    Ptot = Pg + Pr + Pd
    epstot = eos_x1*e_gas(rho,T) + eos_x2*e_rad(rho,T) + eos_x3*e_deg(rho,T)

    #EPS = eos_x1*e_gas(RHO,TTT) + eos_x2*e_rad(RHO,TTT) + eos_x3*e_deg(RHO,TTT)
    #PPP = eos_x1*P_gas(RHO, TTT) + eos_x2*P_rad(RHO, TTT)

    rhoh = rho + rho*epstot + Ptot
    #RHOH = RHO + RHO*EPS + PPP

    W = u0 / np.sqrt(1+2*M/r)

    H = np.sqrt(Ptot*r*r*r/(rhoh*M))/u0
    #HHH = np.sqrt(PPP*R*R*R/(RHOH*M))/U00

    Mdot = -2*math.pi*r*rho*u0*vr*H * c*rg_solar*rg_solar/M_solar

    cs = np.sqrt(cs2(rho,T))
    VR = np.abs((1.0 + 2.0*M/r)*vr + 2.0*M/r)
    VP = np.sqrt(1.0+2.0*M/r)*r*vp
    V = np.sqrt(VR*VR+VP*VP)

    if bounds is None:
        bounds = []
        bounds.append([rho[rho==rho].min(), rho[rho==rho].max()])
        bounds.append([T[T==T].min(), T[T==T].max()])
        bounds.append([Ptot.min(), Ptot.max()])
        bounds.append([(-vr[vr<0]).min(), (-vr[vr<0]).max()])
        bounds.append([vp[vp>0].min(), vp[vp>0].max()])
        bounds.append([W.min(), W.max()])
        bounds.append([H.min(), H.max()])
        bounds.append([Mdot[Mdot>0].min(), Mdot[Mdot>0].max()])
        bounds.append([min(cs.min(),np.abs(VR).min(),np.abs(VP).min(),V.min()),
                    max(cs.max(),np.abs(VR).max(),np.abs(VP).max(),V.max())])
        bounds = np.array(bounds)

    if plot:

        print("Plotting t = {0:g}".format(t))

        #Plot.
        fig = plt.figure(figsize=(12,9))

        plt.subplot(331)
        #plot_r_profile_single(r, rho, sca, r"$\rho_0$", R, RHO)
        plot_r_profile_single(r, rho, sca, r"$\rho_0$ ($g/cm^3$)", bounds[0])

        plt.subplot(332)
        #plot_r_profile_single(r, T, sca, r"$T$", R, TTT)
        plot_r_profile_single(r, T, sca, r"$T$ ($m_p c^2$)", bounds[1])

        plt.subplot(333)
        plt.plot(r, Pg, 'g+')
        plt.plot(r, Pr, 'b+')
        plt.plot(r, Pd, 'r+')
        #plot_r_profile_single(r, Ptot, sca, r"$P$", R, PPP)
        plot_r_profile_single(r, Ptot, sca, r"$P$ ($g\ c^2/cm^3$)", bounds[2])

        plt.subplot(334)
        #plot_r_profile_single(r, vr, "linear", r"$v^r$", R, URR)
        plt.plot(r, -((1-2*M/r)/(1+2*M/r) + A*vp), 'r--')
        plt.plot(r, -(-1*np.ones(r.shape) + A*vp), 'r--')
        plt.plot(r, -((0.5-2*M/r)/(1+2*M/r) + A*vp), 'b--')
        plt.plot(r, -((-0.5-2*M/r)/(1+2*M/r) + A*vp), 'b--')
        plt.plot(r, -((0.0-2*M/r)/(1+2*M/r) + A*vp), ls='--', color='grey')
        plot_r_profile_single(r, -vr, sca, r"$v^r$", bounds[3])
        #plt.gca().set_xscale(sca)

        plt.subplot(335)
        #plot_r_profile_single(r, vp, sca, r"$v^\phi$", R, UPP)
        plt.plot(R, 1.0/(np.sqrt(1+2.0*M/R)*R), 'r--')
        plt.plot(R, 0.5/(np.sqrt(1+2.0*M/R)*R), 'b--')
        plot_r_profile_single(r, vp, sca, r"$v^\phi$", bounds[4])
        OMK = np.sqrt(M/(R*R*R))
        plt.plot(R, OMK, 'g-')
        plt.plot(R, OMK/(1+A*OMK), 'g--')
        plt.plot(R, OMK/((1+A*OMK)*np.sqrt(1+2*M/R)), 'y-')

        plt.subplot(336)
        #plot_r_profile_single(r, mach, sca, r"$\mathcal{M}$", R, MACH)
        plot_r_profile_single(r, W, sca, r"$W$", bounds[5])
        plt.gca().set_yscale('linear')

        plt.subplot(337)
        #plot_r_profile_single(r, H, sca, r"$\mathcal{H}$", R, HHH)
        plot_r_profile_single(r, H, sca, r"$H$", bounds[6])

        plt.subplot(338)
        plot_r_profile_single(r, Mdot, sca, r"$\dot{M}$ ($M_\odot / s$)",
                                bounds[7])

        plt.subplot(339)
        plt.plot(r, VR, 'r+')
        plt.plot(r, VP, 'g+')
        plt.plot(r, V,  'b+')
        plot_r_profile_single(r, cs, sca, r"$c_s$ ($c$)", bounds[8])

        plt.tight_layout()

        outpath = filename.split("/")[:-1]
        chckname = filename.split("/")[-1]
        outname = "plot_disc_{0}.png".format("_".join(chckname.split(".")[0].split("_")[1:]))

        print("Saving {0:s}...".format(outname))
        plt.savefig(outname)

    else:
        fig = None

    return fig, bounds

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("\nGive me a checkpoint (.h5) file.\n")
        sys.exit()

    elif len(sys.argv) == 2:
        filename = sys.argv[1]
        fig = plot_r_profile(filename, sca=scale)
        plt.show()

    else:
        all_bounds = np.zeros((9,2))
        #all_bounds[:,0] = np.inf
        #all_bounds[:,1] = -np.inf
        #for filename in sys.argv[1:]:
            #fig, bounds = plot_r_profile(filename, sca=scale, plot=False)
            #fig, bounds = plot_r_profile(filename, sca=scale, plot=False)
            #all_bounds[:,0] = np.minimum(all_bounds[:,0], bounds[:,0])
            #all_bounds[:,1] = np.maximum(all_bounds[:,1], bounds[:,1])

        for filename in sys.argv[1:]:
            fig, bounds = plot_r_profile(filename, sca=scale, plot=True,
                                        bounds=all_bounds)
            plt.close(fig)
