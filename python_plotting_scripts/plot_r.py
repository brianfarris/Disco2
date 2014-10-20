import math
from numpy import *
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import sys
import numpy as np
import readChkpt as rc

equat = False
if(len(sys.argv) < 2):
    print("\nGive me a checkpoint (.h5) file.\n")
    sys.exit()

filename = sys.argv[1]

dat = rc.readChkpt(filename)
t = dat[0]
r = dat[1]
rho = dat[4]
P = dat[5]
vr = dat[6]
vp = dat[7]

GAM = 1.3333333333

i = r.argmax()
R = np.linspace(r.min(), r.max(), 1000)
RHO = rho[i] * np.power(R/r[i], -0.5)
PPP = P[i] * np.power(R/r[i], -1.5)
URR = vr[i] * np.power(R/r[i], -0.5)
UPP = vp[i] * np.power(R/r[i], -1.5)

print("Plotting t = {0:g}".format(t))

xsca = 'log'
ysca = 'log'

#Plot.
plt.figure(figsize=(10,7))
plt.subplot(321)
plt.plot(r, rho, 'k+')
plt.plot(R, RHO, 'r')
plt.gca().set_xscale(xsca)
plt.gca().set_yscale(ysca)
plt.xlabel(u"$r$")
plt.ylabel(u"$\\rho_0$")
plt.subplot(322)
plt.plot(r, P, 'k+')
plt.plot(R, PPP, 'r')
plt.gca().set_xscale(xsca)
plt.gca().set_yscale(ysca)
plt.xlabel(u"$r$")
plt.ylabel(u"$P$")
plt.subplot(323)
if(vr<0.0).all() and ysca=='log':
    plt.plot(r, -vr, 'k+')
    plt.plot(R, -URR, 'r')
else:
    plt.plot(r, vr, 'k+')
    plt.plot(R, URR, 'r')
plt.gca().set_xscale(xsca)
plt.gca().set_yscale(ysca)
plt.xlabel(u"$r$")
plt.ylabel(u"$v^r$")
plt.subplot(324)
plt.plot(r, vp, 'k+')
plt.plot(R, UPP, 'r')
plt.gca().set_xscale(xsca)
plt.gca().set_yscale(ysca)
plt.xlabel(u"$r$")
plt.ylabel(u"$v^\phi$")
plt.subplot(325)
plt.plot(r, np.sqrt((vr*vr+r*r*vp*vp)/(GAM*P/rho)), 'k+')
plt.gca().set_xscale(xsca)
plt.gca().set_yscale(ysca)
plt.xlabel(u"$r$")
plt.ylabel(u"$M$")

plt.tight_layout()

plt.savefig("plot_r.png")

plt.show()

