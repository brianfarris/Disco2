import matplotlib
from matplotlib import pyplot as plt
import glob
import numpy as np
import sys

filename = sys.argv[1]

r = np.genfromtxt(filename,usecols=0)
rho = np.genfromtxt(filename,usecols=1)
P = np.genfromtxt(filename,usecols=2)
vp = np.genfromtxt(filename,usecols=6)
vr = np.genfromtxt(filename,usecols=3)

alpha = 0.1
r0 = 5.
d = 4
PoRho_r1 = 0.01
rho_init = np.power(r,-3./5.)*np.exp(-np.power(r/r0,-d))
rho_ss = np.power(r,-3./5.)
P_init = PoRho_r1*np.power(r,-3./2.)*np.exp(-np.power(r/r0,-d))
P_ss = PoRho_r1*np.power(r,-3./2.)
vp_ss = -0.75*PoRho_r1*np.power(r,-2./5.)
r_drP_o_P = -1.5 + d*np.power(r/r0,-d);
P_o_rho_init = PoRho_r1*np.divide(np.power(r,-1.5),np.power(r,-3./5.))
vp_init = np.zeros(r.size)
vr_init = np.zeros(r.size)
vr_ss = -1.5*alpha*PoRho_r1*np.power(r,-2./5.)

fig,ax = plt.subplots(4,1,figsize=(4,6))
ax[0].plot(r,rho_init,linewidth=2.)
ax[0].plot(r,rho_ss,linewidth=2.)
ax[0].plot(r,rho,linewidth=2.)
ax[1].plot(r,P_init,linewidth=2.)
ax[1].plot(r,P_ss,linewidth=2.)
ax[1].plot(r,P,linewidth=2.)
ax[2].plot(r,vp_init,linewidth=2.)
ax[2].plot(r,vp_ss,linewidth=2.)
ax[2].plot(r,vp,linewidth=2.)
ax[3].plot(r,vr_init,linewidth=2.)
ax[3].plot(r,vr,linewidth=2.)
ax[3].plot(r,vr_ss,linewidth=2.)
ax[0].set_xlim(1,20)
ax[0].set_ylim(0.0,.99)
ax[1].set_xlim(1,20)
ax[1].set_ylim(0.00001,.0019)
ax[2].set_xlim(1,20)
ax[2].set_ylim(-.009,0.0019)
ax[3].set_xlim(1,20)
ax[3].set_ylim(-.0029,0.0004)
ax[3].set_xlabel(r'$r/a$',size='x-large')
ax[0].set_ylabel(r'$\rho$',labelpad=20.,size='x-large',rotation='horizontal')
ax[1].set_ylabel(r'$P$',labelpad=20.,size='x-large',rotation='horizontal')
ax[2].set_ylabel(r'$\delta v_{\phi}$',labelpad=20.,size='x-large',rotation='horizontal')
ax[3].set_ylabel(r'$v_r$',labelpad=20.,size='x-large',rotation='horizontal')
ax[0].get_xaxis().set_visible(False)
ax[1].get_xaxis().set_visible(False)
ax[2].get_xaxis().set_visible(False)
y_formatter = matplotlib.ticker.FormatStrFormatter('%2.0e')
ax[1].yaxis.set_major_formatter(y_formatter)
ax[2].yaxis.set_major_formatter(y_formatter)
ax[3].yaxis.set_major_formatter(y_formatter)
plt.tight_layout()
plt.subplots_adjust(hspace = 0.0)
plt.savefig("SS1D.pdf")
