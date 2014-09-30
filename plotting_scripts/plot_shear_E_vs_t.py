from matplotlib import pyplot as plt
import numpy as np
from scipy import integrate

t0 = 0.25

t = np.genfromtxt("DiagScalar.dat",usecols=0)
kinetic = np.genfromtxt("DiagScalar.dat",usecols=10)
internal = np.genfromtxt("DiagScalar.dat",usecols=11)
radiation_rate = np.genfromtxt("DiagScalar.dat",usecols=7)
radiated = np.insert(integrate.cumtrapz(radiation_rate,t),0,0.0)
total = np.genfromtxt("DiagScalar.dat",usecols=9)+radiated

fig,ax = plt.subplots(figsize=(8,8))
ax.plot(t/t0,kinetic,label='kinetic energy',linewidth=3.)
ax.plot(t/t0,internal,label='internal energy',linewidth=3.)
ax.plot(t/t0,radiated,label='radiated energy',linewidth=3.)
ax.plot(t/t0,total,label='total energy',linewidth=3.)


# set labels and titles
ax.set_xlabel(r'$t/t_0$',size='x-large')
ax.set_ylabel(r'$E$',rotation='horizontal',size='x-large')

# create legend
legend = ax.legend(loc='upper right', shadow=True)

ax.set_ylim(0,0.18)

plt.savefig('shear_E_vs_t.pdf')
