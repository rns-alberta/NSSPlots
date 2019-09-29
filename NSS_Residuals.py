#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from sklearn.preprocessing import PolynomialFeatures
from sklearn import linear_model
import numpy.polynomial.polynomial as poly

fig = plt.figure(figsize = (9,7))
ax1 = fig.add_subplot(111, projection='3d')

filename = 'NS_data_eosBBB1.txt'
data = np.loadtxt(str(filename), unpack=True)

fil = filename.replace('NS_data_eos', '')
name = fil.replace('.txt', '')

# Columns of data from the previous file 
Ec = data[0,:]*10**15       # central energy density
M = data[1,:]               # Total mass
M0 = data[2,:]              # Baryonic mass
Mstat = data[3,:]           # Mass when the NS is not rotating 
Mmax = data[4,:]            # Maximum mass for a given EOS
R = data[5,:]               # Radius 
Rratio = data[6,:]          # Ratio rp/re
Rstat = data[7,:]           # Radius when the NS is not rotating
Omega = data[8,:]*2*np.pi   # Angular velocity (rad/s), since data[8,:] is rotational frequency (in Hz)
Klim = data[9,:]            # Kepler limit (Hz)
J = data[10,:]              # Angular momentum
T = data[11,:]              # Rotational kinetic energy
W = data[12,:]              # Gravitational binding energy
Rmax = data[13,:]           # Maximum radius for a given EOS
Qmax = data[14,:]           # Ratio (MaxMass / MaxRadius) of the non rotating star

# Some constants
G = 6.67e-8          # in cm^3 g^-1 s^-2
c = 2.99792458e10    # in cm s^-1
Msun = 1             # in solar masses

# Converting quantities to CGS units
m = M * 1.9884e33   # Mass
r = R * 1.0e5       # Radius
rstat = Rstat * 1.0e5   # Radius of the non rotating NS
mstat = Mstat * 1.9884e33   # Mass of the non rotating NS
delta = Omega * ((rstat**3)/(G*mstat))**(0.5)   # Normalized omega
m0 = M0 * 1.9884e33     # Baryonic mass
normJ = (c*J)/(G*m0**2)     # Normalized angular momentum

print('Choose the residual plor to generate (Choose only integer numbers)')
print('1. For (M0-M)/M0')
print('2. For (M-M*)/M*')
print('3. For (R-R*)/R*')

# Choose the ssurface plot to make
method = input()
#method = 1

fil = filename.replace('NS_data_eos', '')
name = fil.replace('.txt', '')

if method == 1:
    x = delta**2
    y = G*mstat/(rstat*c**2)  
    z = (M0-M)/M0

    resid = np.zeros(len(x))

    for i in range(len(x)):
        Xc = x[i]
        Yc = y[i]
        term1 =  -1.51950180e-02 -2.10584774e-01*Xc + 8.92558965e-01*Yc
        term2 = 1.70935907*Xc**2 + 1.38563605*Xc*Yc -3.72546017*Yc**2
        term3 = -5.45957168e+00*Xc**3 -6.44297945e+00*Xc**2*Yc -3.94675715e+00*Xc*Yc**2 + 1.59905891e+01*Yc**3  
        term4 = 5.90617974e+00*Xc**4 + 9.34839417e+00*Xc**3*Yc + 3.94334721e+00*Xc**2*Yc**2 + 3.31392052e+00*Xc*Yc**3 -2.20445973e+01 *Yc**4
        Zc =  term1 + term2 + term3 + term4

        resid[i] =  100*(z[i] - Zc)
    
    ax1.scatter(x, y, resid, marker='o', label=str(name))
    ax1.set_xlabel(r'$\Omega^2 (R_*^3 / GM_*)$', fontsize='15')
    ax1.set_ylabel(r'$GM_*/R_*c^2$', fontsize='15')
    ax1.set_zlabel(r'Residuals in $(M_0-M)/M_0$', fontsize='15')
    ax1.view_init(azim=220)
    #ax1.view_init(elev=0)
    plt.legend()
    plt.savefig('M0-frac-residuals.png', format = 'png', transparent=False)
    plt.show()


elif method == 2:
    x = delta**2
    y = G*mstat/(rstat*c**2)
    z = (M-Mstat)/Mstat

    resid = np.zeros(len(x))

    for i in range(len(x)):
        Xc = x[i]
        Yc = y[i]
        term1 =  -0.00368173 + 0.10930922*Xc + 0.033793*Yc
        term2 = -0.96591924*Xc**2 -0.27353663*Xc*Yc -0.15198814*Yc**2
        term3 = 2.90307002*Xc**3 +  2.73955588*Xc**2*Yc -0.25523957*Xc*Yc**2 + 0.49056401*Yc**3  
        term4 = -2.97802269*Xc**4 -3.59644437*Xc**3*Yc + 0.08914017*Xc**2*Yc**2 + 2.69291532*Xc*Yc**3 -0.72895654*Yc**4
        Zc =  term1 + term2 + term3 + term4

        resid[i] =  100*(z[i] - Zc)

    ax1.scatter(x, y, resid, marker='o', label=str(name))
    ax1.set_xlabel(r'$\Omega^2 (R_*^3 / GM_*)$', fontsize='15')
    ax1.set_ylabel(r'$GM_*/R_*c^2$', fontsize='15')
    ax1.set_zlabel(r'Residuals in $(M_*-M)/M_*$', fontsize='15')
    ax1.view_init(azim=220)
    #ax1.view_init(elev=0)
    plt.legend()
    plt.savefig('Mstat-frac-residuals.png', format = 'png', transparent=False)
    plt.show()


elif method == 3:
    x = delta**2
    y = (M/R)/Qmax
    z = (R-Rstat)/Rstat

    resid = np.zeros(len(x))

    for i in range(len(x)):
        Xc = x[i]
        Yc = y[i]
        term1 =  -1.0790784  + 10.21820411*Xc +7.23015196*Yc
        term2 = -44.71532998*Xc**2 -39.82600955*Xc*Yc -16.94487249*Yc**2
        term3 = 121.99681083*Xc**3 + 73.35172903*Xc**2*Yc + 57.45075819*Xc*Yc**2 + 16.52530708 *Yc**3  
        term4 = -118.46931729*Xc**4 -63.91197667*Xc**3*Yc -40.44452759*Xc**2*Yc**2 -27.00592591*Xc*Yc**3 -5.72330642*Yc**4
        Zc =  term1 + term2 + term3 + term4
        
        resid[i] =  100*(z[i] - Zc)

    ax1.scatter(x, y, resid, marker='o', label=str(name))
    ax1.set_xlabel(r'$\Omega^2 (R_*^3 / GM_*)$', fontsize='15')
    ax1.set_ylabel(r'$(M/R)/(M_M/R_M)$', fontsize='15')
    ax1.set_zlabel(r'Residuals in $(R-R_*)/R_*$', fontsize='15')
    ax1.view_init(azim=220)
    #ax1.view_init(elev=0)
    plt.legend()
    plt.savefig('Rstat-frac-residuals.png', format = 'png', transparent=False)
    plt.show()


else:
    print('Enter an integer between 1 and 3.')
