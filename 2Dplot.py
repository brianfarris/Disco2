import numpy 
from numpy import *
import scipy 
from scipy import *
import pylab 
from pylab import*
#import matplotlib.pyplot as plt
import h5py
from h5py import *
import triangle 
import triangle.plot

#from triangle import *



# h5 file under key 'EQUAT'
# 19 indices
#0 - r
#1 - phi
#2 - density all phi at each r
#3 - vr
#4 - vphi
#5 - Pressure
#...

###GLOBAL PARAMS
###PHYSICAL PARAMS
#initial binary separation
abin0 = 1.0
# viscous paramter
alpha=0.000
#Mach number in iso disk
Mach = 10.
#Binary mass ratio
q = 0.0000026125 #0.00001
#Disk mass draction of secondary
MdoMs = 1.00
# Sig0
Sig0 = MdoMs/(1.+1./q)/(pi*abin0*abin0)
#When to turn the live binary on
tmig_on = 0.0


###PLOTTING PARAMS
#plot log density?
logdens = False
#plot velocity vectors over density?
V_vectors = False
# color map
cmapd = 'cm.gray_r'
cmapdLog = 'cm.gist_heat'

####Read in h5 file(s)
fname="DiagEquat_"
#fname = "../LiveBinary/DiagEquat_circKepDecay_q0p001_MdoMs0p01_r128_phiNcap64_"
#fname = "../LiveBinary/DiagEquat_circKepDecay_q0p001_MdoMs0p01_r128_phiNcap128_"



#tsimA= 519.0039
#tsimA = 73.0041
#tsimA = 219.0097
tsimA = 75.0026

fnameA = fname + '%.4f' %tsimA + '.h5'
datA  = h5py.File(fnameA, 'r')


###Make the 2D grid of points
datT=transpose(datA['EQUAT'])
rr  = datT[0]
phi = datT[1]

#rph = transpose([rr,phi])
#rpdict = {'rphi':rph}

xx=multiply(rr, cos(phi))
yy=multiply(rr, sin(phi))

skip = 2
Nq=len(xx)/skip
xq=[]
yq=[]

for i in range(Nq):
  xq.append(xx[i+skip])
  yq.append(yy[i+skip])   

xy = [xx,yy]


### Make a triangulation of points for tricontour
xyt = transpose(xy)
A = dict(vertices=xyt)
B = triangle.triangulate(A)
triangles = B['triangles']

 #### Get the density and the velocity components 
 #### the fluid variable arrays are arranged so that
 #### they are 1D give each azimuthal value for a 
 ####given r and going outwards in r, thus the array below is Nphi*Nr long
dens = datT[2]
velr = datT[3]
velp = datT[4]

velx = velr * cos(phi) - velp * sin(phi)
vely = velp * cos(phi) + velr * sin(phi)



###POSTIONS OF BH's###
####################################                                            
### FILE NAME ###                                                               
file1 = "BinaryParams.dat"

#64 az res cap
#file1 = "../LiveBinary/BinaryParams_circKepDecay_q0p001_MdoMs0p01_r128_phiNcap64.dat"
#128 az res cap
#file1 = "../LiveBinary/BinaryParams_circKepDecay_q0p001_MdoMs0p01_r128_phiNcap128.dat"

####################################                                            
### READ IN DATA ###                                                            
####################################                                            
t       = loadtxt(file1, usecols=[0])
r_BH0   = loadtxt(file1, usecols=[1])
r_BH1   = loadtxt(file1, usecols=[2])
a_bin   = loadtxt(file1, usecols=[3])
phi_BH0 = loadtxt(file1, usecols=[4])
phi_BH1 = loadtxt(file1, usecols=[5])
ecc     = loadtxt(file1, usecols=[6])
Enrgy   = loadtxt(file1, usecols=[7])
Om      = loadtxt(file1, usecols=[11])
vrp     = loadtxt(file1, usecols=[10])
Tr      = -loadtxt(file1, usecols=[14])
dSTr    = loadtxt(file1, usecols=[15])



##abin = r_BH0 + r_BH1
norb = t/(2.*pi*a_bin**(3./2.))

### find the separation at the timestep tsimA when we are plotting ###
indexA = where(t <= tsimA )
iA=len(indexA[0])-1



tratio = tsimA/(t[iA])
tratiop1 = tsimA/t[iA+1]



if ( fabs(1.-tratio) > fabs(1.-tratiop1) ):
### polar coordinates of each BH
#phi0 = Omega_bin * tsimA
#phi0 = .5*(phi_BH0[iA]+ phi_BH0[iA+1])
	phi0 = phi_BH0[iA+1]
	phi1 = phi_BH1[iA+1]
#phi1 = phi0 + pi
	rbh0 = r_BH0[iA+1] 
	rbh1 = r_BH1[iA+1]
	### binary angular frequency at the current separation
	Omega_bin = sqrt(1./(a_bin[iA+1]*a_bin[iA+1]*a_bin[iA+1]))
	print "At index +1 ->  %.6f" %tratiop1
else:
	phi0 = phi_BH0[iA]
	phi1 = phi_BH1[iA]
#phi1 = phi0 + pi
	rbh0 = r_BH0[iA] 
	rbh1 = r_BH1[iA]
	### binary angular frequency at the current separation
	Omega_bin = sqrt(1./(a_bin[iA]*a_bin[iA]*a_bin[iA]))
	print "At index ->  %.6f" %tratio

	


### Cartesian coordinates of each BH
xbh0 = rbh0*cos(phi0)
ybh0 = rbh0*sin(phi0)
xbh1 = rbh1*cos(phi1)
ybh1 = rbh1*sin(phi1)

### DEVIATION FROM CM
#CMdev = 2.*cos(phi_BH1 - phi_BH0)




### e=0 barycentric a
##Om_bin = 1./(a_bin)**(3./2.)
akep= - 1./(1.+q)*1./(1.+1./q)*1./(2.*Enrgy)
## This will tell us if the code is updating a like we want but it wont 
## tell us if E is being updated properly by the caluclated forces 
## - we need another diagnostic for this - analytic migration rates
#Type I
rs = a_bin/(1.+q)
rs0 = abin0/(1.+q)
Om_disk_s = 1./(rs)**(3./2.)
csnd = sqrt(1./rs)/Mach
Mp = 1.0/(1.+q)
Ms = 1.0/(1.+1./q)
c1 = 1.56 #Numerical param see Ward 1997 #for Mplanet/Msun -> 0
vrI = -c1*q*(rs*Om_disk_s)*Sig0*rs*rs/Mp * (rs*Om_disk_s/csnd)**3

 
#Type II (Disk Dom)
c2= 1.5 #for Mplanet/Msun -> infinity
vrII = -c2*alpha*(rs*Om_disk_s)*((rs*Om_disk_s)/csnd)**(-2)
#Del0 = 0.1 #~Gap width
#Htdot = 8./27. * (Ms/Mp)**2 * Sig0 * rs**4 * Om_bin**2 * (rs/Del0)**3
#trII_LP86 = Ms*Om_bin*rs*rs/Htdot






#-------------------------------------------------------------------#
### Simulation Torques ###
### Int[ dSig * r * f_phi dr
#-------------------------------------------------------------------#
f1DVec = "DiagVector_%.4f.dat" %tsimA	

#64 az res cap
#f1DVec = "../LiveBinary/DiagVector_circKepDecay_q0p001_MdoMs0p01_r128_phiNcap64_%.4f.dat" %tsimA
#128 az res cap
#f1DVec = "../LiveBinary/DiagVector_circKepDecay_q0p001_MdoMs0p01_r128_phiNcap128_%.4f.dat" %tsimA


f1DSCal = "DiagScalar.dat" 

#64 az res cap	
#f1DSCal = "../LiveBinary/DiagScalar_circKepDecay_q0p001_MdoMs0p01_r128_phiNcap64.dat"
#128 az res cap
#f1DSCal = "../LiveBinary/DiagScalar_circKepDecay_q0p001_MdoMs0p01_r128_phiNcap128.dat"


#SigBkg = 1.0
#SigAzAv = loadtxt(f1DVec, usecols=[1]) ###Is this an azavg only? is it time avg?
#Dsig = SigAzAv - SigBkg 

# torque denisty from sim output
DTdr = loadtxt(f1DVec, usecols=[12])
dSDTdr = loadtxt(f1DVec, usecols=[20])

#time and integrated Torque from sim output
tScal =  loadtxt(f1DSCal, usecols=[0])  ## t in DiagScalar.dat
Tr_sp =  -loadtxt(f1DSCal, usecols=[12])  ## r integrated torque using Sig  - defined from bin params now
dSTr_sp =  loadtxt(f1DSCal, usecols=[20])  ## r integrated torque using deltaSig_phi - defined from bin params now
Tr_BHL = -loadtxt(f1DSCal, usecols=[18])  ## eqn 4 of BHL 2013
Mdsk_rs = Sig0*loadtxt(f1DSCal, usecols=[19])  ## eqn 4 of BHL 2013





### find the separation at the time tScal which is at the time where DigScalar dumps ###
Match = list(set(t) & set(tScal))
iS = []
for i in range(len(Match)):
	indS = where(t <= Match[i])  
	iS.append(len(indS[0]) - 1)

iS = sort(iS)

norb_Scal = tScal/(2*pi*(a_bin[iS])**(3./2.))


#Type II (Secondary Dom)
nus=alpha*csnd*csnd/Om_disk_s
trII_SD = rs[iS]*rs[iS]/nus[iS] /(6.*pi*(1.+q)) * 1./Mdsk_rs #Use simulation value of disk mass out to rs


## unpack the 1D radial array
r1D = []
r1D.append(rr[0]) 
for i in range(len(rr) -1):
	if (rr[i+1] - rr[i] != 0.0): 
		r1D.append(rr[i+1])	

#integrate in r
##DT=0.0
##for i in range(len(r1D)-1):
##	DT += (DTdr[i+1] -DTdr[i])*(r1D[i+1] - r1D[i]) ### AT ONLY ONE TIMESTEP THOUGH!
		

		
#NEED TO have rs, Om at same time steps as Tr
###vs_simTrq = 2.*Tr/(Ms*rs*Om_disk_s)
#vs_simTrq = Sig0*2.*(dSTr)/(Ms*a_bin[iS]/(1.+q) * Om_bin[iS] )
#vs_simTrq = Sig0*2.*(Tr)/( Ms*a_bin[iS]/(1.+q) * Om_bin[iS] )

###vs_simTrq = (-Sig0*Tr/Ms - a_bin[iS]/(1.+q)*a_bin[iS]/(1.+q)*  (Om_bin[iS]-Om_bin[iS-1])/(t[iS]-t[iS-1]))/(2.*a_bin[iS]/(1.+q) * Om_bin[iS])

## Assuming that Omega is Keplerian
for i in range(len(tScal)):
	if (tScal[i]<tmig_on):
		dSTr_sp[i] = 0.0
		Tr_sp[i] = 0.0
		Tr_BHL[i] = 0.0
		
for i in range(len(t)):
	if (t[i]<tmig_on):
		dSTr[i] = 0.0
		Tr[i] = 0.0

 
vs_simTrq = (2.*Sig0*Tr)/(Ms*a_bin/(1.+q) * Om) #/ (10.*10. - 0.05*0.05)*2.
vs_simTrq_sparse = (2.*Sig0*Tr_sp/((10.*10. - 0.05*0.05)*2.))/(Ms*a_bin[iS]/(1.+q) * Om[iS])
vs_BHL = (2.*Sig0*1.*Tr_BHL)/(Ms*a_bin[iS]/(1.+q) * Om[iS])
chkTr = Sig0*dSTr/Ms
#chkOm = a_bin[iS]/(1.+q)*a_bin[iS]/(1.+q)*(Om[iS]-Om[iS-1])/(t[iS]-t[iS-1])
#-------------------------------------------------------------------#



### Distance as a function of time
aTypeI = abin0 + t*vrI
#aTypeII_DD = rs0 +  t*vrII
#aTypeII_LP = abin0*(1. - t/trII_LP86)
#aTypeII_SD = rs0*(1. - t/trII_SD)

rsTypeII_DD  = zeros(len(t))
rsTypeII_DD[0] = rs0
rsTypeII_DD[1] = rsTypeII_DD[0] + t[0]*vrII[0]

for i in range(1, len(t)-1):
	rsTypeII_DD[i+1] = rsTypeII_DD[i] + (t[i]-t[i-1])*vrII[i]


rsTypeII_SD  = zeros(len(t))
rsTypeII_SD[0] = rs0
rsTypeII_SD[1] = rsTypeII_SD[0]*(1. - t[0]/trII_SD[0])

for i in range(1, len(tScal)-1):
	rsTypeII_SD[i+1] = rsTypeII_SD[i]*(1. - (t[i]-t[i-1])/trII_SD[i] )


rs_SimTrq  = zeros(len(t))
rs_SimTrq[0] = rs0
rs_SimTrq[1] = rs_SimTrq[0] + t[0]*vs_simTrq[0]

rs_SimTrq_sp  = zeros(len(tScal))
rs_SimTrq_sp[0] = rs0
rs_SimTrq_sp[1] = rs_SimTrq_sp[0] + tScal[0]*vs_simTrq_sparse[0]

rs_BHL  = zeros(len(tScal))
rs_BHL[0] = rs0
rs_BHL[1] = rs_SimTrq[0] + tScal[0]*vs_BHL[0]

for i in range(1, len(t)-1):
	rs_SimTrq[i+1] = rs_SimTrq[i] + (t[i]-t[i-1])*vs_simTrq[i]
	
for i in range(1, len(tScal)-1):
	rs_BHL[i+1] = rs_BHL[i] + (tScal[i]-tScal[i-1])*vs_BHL[i]
	rs_SimTrq_sp[i+1] = rs_SimTrq_sp[i] + (tScal[i]-tScal[i-1])*vs_simTrq_sparse[i]
 
#rs_SimTrq = rs0 +  tScal*vs_simTrq
#a_SimTrq = rs0 + tScal*vs_simTrq

#rsTypeII_SD = rs*(1. -  t/trII_SD)
###-------------------------------------###



###PLOTTING###
figure(figsize=[13,9])
subplot(221)
#CT = tripcolor(xx,yy,triangles,dens)
if (logdens==True):
    CT = tricontourf(xx,yy,triangles,log10(dens),200, cmap=cm.gist_heat, colorbar=True)
else:
    CT = tricontourf(xx,yy,triangles,dens,200, cmap=cm.gray_r, colorbar=True) 
CB = colorbar(CT)
if (V_vectors == True):
  quiver(xx,yy,velx,vely, color='gray')#, scale=20.)

### Plot the positions of the BH's
scatter(xbh0, ybh0, color="black")
scatter(xbh1, ybh1, color="black")
xlabel('x [r/a]')
ylabel('y [r/a]')

#xlim(-20,20)
#ylim(-20,20)

#xlim(-5,5)
#ylim(-5,5)

subplot(222)
##plot(norb, (r_BH0+r_BH1), color='black', linestyle='--')
#plot(norb, aTypeI, color='blue')
#plot(norb, rsTypeII_DD, color='green')
#plot(norb, rsTypeII_SD, color='red', linestyle='--')
#plot(norb, aTypeII_LP, color='blue')
#plot(norb,r_BH0, color='green')
#plot(norb, rs_BHL, color='green')

p1=plot(norb,r_BH1, color='black')
#p1=scatter(norb,r_BH1, marker='x', s=10, color='black')
#p2=scatter(norb_Scal, rs_SimTrq, marker='x', s=20, color='red')
#p2=plot(norb, rs_SimTrq, color='red')
plot(norb_Scal, rs_SimTrq_sp, color='green', linestyle=':')
xlim(tmig_on/(2.*pi), max(norb))

#ylim(min(r_BH1), max(r_BH1))
#figlegend((p1,p2),("$r^s_{sim}$","$r^s_{Trq}$"), (0.817,0.807))

### rs instead of a
#plot(norb, r_BH1, color='black')
#plot(norb, rsTypeII_SD, color='red', linestyle='--')

#ylim(.9996*a_bin[len(t)-1], 1.0004)
##axhline(y=abin0, color='black', linestyle=":")
ylabel('$r^s$')
xlabel('$N_{orb}$')
#xlim(0,100)

subplot(223)
#plot(norb,Om, color='black')
#plot(norb,Om_bin, color='green')
#plot(r1D, DTdr)
#plot(norb_Scal, chkD)
#plot(r1D, -dSDTdr)
plot(r1D, DTdr)
##plot(norb,vrp+abin0, color='black')
##plot(norb,phi_BH0, color='black')
#xlabel('$N_{orb}$')
xlabel('$r/a_0$')
ylabel('$T_r$')

#plot(norb, ecc)
#ylabel('e')
#xlim(abin0-5/Mach,abin0+5./Mach)

subplot(224)
#plot(norb, Enrgy/Enrgy[0], color='red')
#plot(norb, vs_simTrq, color='black')
plot(norb_Scal, vs_simTrq_sparse)
#plot(norb, chkOm, color='red')
#ylabel('E/$E_0$')
xlabel('$N_{orb}$')
#xlim(0,100)


show()
#savefig("../LiveBinary/LiveBin_Rk2_r10_128_PhiCap64_AllFon_tsim%.4f.png" %tsimA)
#savefig("../LiveBinary/LiveBin_circKepDecay_q0p001_MdoMs0p01_r128_phiNcap64_tsim%.4f.png" %tsimA)




