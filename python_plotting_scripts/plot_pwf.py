import math
from numpy import *
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import sys
import numpy as np
import discopy as dp

GAM = 5.0/3.0
M = 3.0
scale = 'log'
x1 = 1.0
x2 = 1.0
x3 = 1.0
alpha = 0.1

c = 2.99792458e10
G = 6.6738e-8
h = 6.62606957e-27
kb = 1.3806488e-16
sb = 1.56055371e59
mp = 1.672621777e-24
rg_solar = 1.4766250385e5
M_solar = 1.9884e33

def xnuc(rho, T):
    rho10 = rho * 1.0e-10
    T10 = T * mp*c*c/kb * 1.0e-10
    xn = 30.97*np.power(rho10,-0.75)*np.power(T10,1.125)*np.exp(-6.096/T10)
    xn[xn>1.0] = 1.0
    return xn

def P_gas(rho, T):
    xn = xnuc(rho, T)
    return (0.25+0.75*xn) * rho * T

def P_rad(rho, T):
    return (11.0/12.0) * 4.0*sb/c * (T*mp*c*c)**4 / (c*c)

def P_deg(rho, T):
    return 2*np.pi*h*c/3.0 * np.power(3*rho/(16*np.pi*mp),4.0/3.0) / (c*c)

def e_gas(rho, T):
    xn = xnuc(rho, T)
    return (0.25+0.75*xn)/(GAM-1) * T

def e_rad(rho, T):
    return (11.0/4.0) * 4.0*sb/c * (T*mp*c*c)**4 / (rho*c*c)

def e_deg(rho, T):
    return h*c/(4.0*mp) * np.power(3*rho/(16*np.pi*mp),1.0/3.0) / (c*c)

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

    plt.xlim(r.min(), r.max())
    plt.axvspan(plt.xlim()[0], 2*M*rg_solar, color='grey', alpha=0.5)

    """
    if bounds is not None:
        if sca == "log":
            upper = 10.0 ** (math.ceil(math.log10(bounds[1]))+0.1)
            if bounds[0] > 0:
                lower = 10.0 ** (math.floor(math.log10(bounds[0]))-0.1)
            else:
                lower = 10.0 ** (math.floor(math.log10(f[f>0].min()))-0.1)
            plt.ylim(lower, upper)
    """

def plot_r_profile(filename, sca='linear', plot=True, bounds=None):

    print("Reading {0:s}".format(filename))

    dat = dp.readCheckpoint(filename)
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

    #i = r.argmax()
    #R = np.linspace(r.min(), r.max(), 100)
    #R = np.logspace(np.log10(r.min()), np.log10(r.max()), 100)
    #RHO = rho[i] * np.power(R/r[i], -1.5)
    #TTT = T[i] * np.power(R/r[i], -1.0)
    #URR = vr[i] * np.power(R/r[i], -0.5)
    #UPP = vp[i] * np.power(R/r[i], -1.5)

    #u0 = 1.0/np.sqrt(1.0-2*M/r-vr*vr/(1-2*M/r)-r*r*vp*vp)
    #U00 = 1.0/np.sqrt(1.0-2*M/R-URR*URR/(1-2*M/R)-R*R*UPP*UPP)
    u0 = 1.0/np.sqrt(1.0-2*M/r-4*M/r*vr-(1+2*M/r)*vr*vr-r*r*vp*vp)
    #U00 = 1.0/np.sqrt(1.0-2*M/R-4*M/R*URR-(1+2*M/R)*URR*URR-R*R*UPP*UPP)

    VR = vr*u0/np.sqrt(1.0 - 2.0*M/r + u0*u0*vr*vr) * c
    VP = vp / (1.0 - 2*M*vr/(r-2*M)) * c/rg_solar
    l = r*r*u0*vp * c*rg_solar

    #fac = 1.0/(1.0-2*M/(r-2*M)*vr)
    #vr = vr*fac
    #vp = vp*fac

    #V = u0*vr / np.sqrt(1-2*M/r+u0*u0*vr*vr)
    #VVV = U00*URR / np.sqrt(1-2*M/R+U00*U00*URR*URR)

    #mach = np.sqrt((vr*vr/(1-2*M/r)+r*r*vp*vp)/(GAM*P/rhoh))
    #MACH = np.sqrt((URR*URR/(1-2*M/R)+R*R*UPP*UPP)/(GAM*PPP/RHOH))

    Pg = x1*P_gas(rho, T)
    Pr = x2*P_rad(rho, T)
    Pd = x3*P_deg(rho, T)
    Ptot = Pg + Pr + Pd
    epstot = x1*e_gas(rho,T) + x2*e_rad(rho,T) + x3*e_deg(rho,T)

    PPtot = Ptot * c*c

    #EPS = x1*e_gas(RHO,TTT) + x2*e_rad(RHO,TTT) + x3*e_deg(RHO,TTT)
    #PPP = x1*P_gas(RHO,TTT) + x2*P_rad(RHO,TTT) + eox_x3*P_deg(RHO,TTT)

    rhoh = rho + rho*epstot + Ptot
    #RHOH = RHO + RHO*EPS + PPP

    Tk = T *mp*c*c / kb

    W = 1.0/np.sqrt(1.0 - ((1.0+2.0*M/r)*vr+2*M/r)**2 - (1.0+2.0*M/r)*r*r*vp*vp)
    H = np.sqrt(Ptot*r*r*r/(rhoh*M))/u0
    HH = np.sqrt(Ptot*r*r*r/(rhoh*M))/u0 * rg_solar
    #HHH = np.sqrt(PPP*R*R*R/(RHOH*M))/U00

    Mdot = -2*math.pi*r*rho*u0*vr*H * c*rg_solar*rg_solar/M_solar
    xn = xnuc(rho, T)

    qee = 5.0e33 * np.power(Tk/1.0e11,9.0)
    qeN = 9.0e33 * xn * (rho/1.0e10) * np.power(Tk/1.0e11,6.0)
    q = qee + qeN

    Lee = qee * 2*math.pi*r*H * rg_solar*rg_solar
    LeN = qeN * 2*math.pi*r*H * rg_solar*rg_solar
    L   = q   * 2*math.pi*r*H * rg_solar*rg_solar

    if bounds is None:
        bounds = []
        bounds.append([rho.min(), rho.max()])
        bounds.append([Tk.min(), Tk.max()])
        bounds.append([(-VR[VR<0]).min(), (-VR[VR<0]).max()])
        bounds.append([HH.min(), HH.max()])
        bounds.append([VP[VP>0].min(), VP[VP>0].max()])
        bounds.append([l[l>0].min(), l[l>0].max()])
        bounds.append([q.min(), q.max()])
        bounds.append([Mdot[Mdot>0].min(), Mdot[Mdot>0].max()])
        bounds.append([L.min(), L.max()])
        bounds.append([PPtot.min(), PPtot.max()])
        bounds.append([0.0,1.0])
        bounds = np.array(bounds)

    if plot:

        print("Plotting t = {0:g}".format(t))

        #Plot 1.
        r *= rg_solar

        fig1 = plt.figure(figsize=(9,9))

        plt.subplot(321)
        plot_r_profile_single(r, rho, sca, r"$\rho_0$ $(g/cm^3)$", bounds[0])

        plt.subplot(322)
        plot_r_profile_single(r, Tk, sca, r"$T$ $(K)$", bounds[1])

        plt.subplot(323)
        plot_r_profile_single(r, -VR, sca, r"$V$ $(cm/s)$", bounds[2])

        plt.subplot(324)
        plot_r_profile_single(r, HH, sca, r"$H$ $(cm)$", bounds[3])

        plt.subplot(325)
        #plot_r_profile_single(r, VP, sca, r"$\Omega$ ($1/s$)", bounds[4])
        plot_r_profile_single(r, VP, sca, r"$\Omega$ $(1/s)$", None)
        plt.plot(r, np.power(r/(M*rg_solar),-1.5)*c/(M*rg_solar), color='grey')

        plt.subplot(326)
        plot_r_profile_single(r, l, sca, r"$\ell$ $(cm^2/s)$", bounds[5])

        plt.tight_layout()

        #Plot 2.
        fig2 = plt.figure(figsize=(9,9))

        plt.subplot(321)
        plt.plot(r, qee, 'g+')
        plt.plot(r, qeN, 'b+')
        plot_r_profile_single(r, q, sca, r"$\dot{q}$ $(erg/cm^3s)$", bounds[6])

        plt.subplot(322)
        plot_r_profile_single(r, Mdot, sca, r"$\dot{M}$ $(M_{\odot}/s)$", 
                                bounds[7])

        plt.subplot(323)
        plt.plot(r, Lee, 'g+')
        plt.plot(r, LeN, 'b+')
        plot_r_profile_single(r, L, sca, r"$dL/dr$ $(erg/cm\ s)$", bounds[8])

        plt.subplot(324)
        plt.plot(r, Pg*c*c, 'r+')
        plt.plot(r, Pr*c*c, 'g+')
        plt.plot(r, Pd*c*c, 'b+')
        plot_r_profile_single(r, PPtot, sca, r"$P$ $(erg/cm^3)$", bounds[9])
        
        plt.subplot(325)
        plot_r_profile_single(r, xn, "linear", r"$X_{nuc}$", bounds[10])
        plt.gca().set_xscale("log")

        plt.tight_layout()


        outpath = filename.split("/")[:-1]
        chckname = filename.split("/")[-1]
        outname1 = "plot_pwf1_{0}.png".format("_".join(chckname.split(".")[0].split("_")[1:]))
        outname2 = "plot_pwf2_{0}.png".format("_".join(chckname.split(".")[0].split("_")[1:]))


        print("Saving {0:s} and {1:s}...".format(outname1, outname2))
        fig1.savefig(outname1)
        fig2.savefig(outname2)
        figs = (fig1,fig2)

    else:
        figs = None

    return figs, bounds

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("\nGive me a checkpoint (.h5) file.\n")
        sys.exit()

    elif len(sys.argv) == 2:
        filename = sys.argv[1]
        figs = plot_r_profile(filename, sca=scale)
        plt.show()

    else:
        all_bounds = np.zeros((11,2))
        all_bounds[:,0] = np.inf
        all_bounds[:,1] = -np.inf
        for filename in sys.argv[1:]:
            figs, bounds = plot_r_profile(filename, sca=scale, plot=False)
            all_bounds[:,0] = np.minimum(all_bounds[:,0], bounds[:,0])
            all_bounds[:,1] = np.maximum(all_bounds[:,1], bounds[:,1])

        for filename in sys.argv[1:]:
            figs, bounds = plot_r_profile(filename, sca=scale, plot=True,
                                        bounds=all_bounds)
            plt.close(figs[0])
            plt.close(figs[1])
