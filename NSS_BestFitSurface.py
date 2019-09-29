#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
#from matplotlib import cm
#from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from sklearn.preprocessing import PolynomialFeatures
from sklearn import linear_model

# Size of the figure
fig = plt.figure(figsize = (9,7))

filename = 'NS_data_eosBBB1.txt'  #NS_data_eosBBB1
data = np.loadtxt(str(filename), unpack=True)

# Columns of data from the previous file 
Ec = data[0,:]*10**15       # central energy density
M = data[1,:]               # Total mass
M0 = data[2,:]              # Baryonic mass
Mstat = data[3,:]           # Mass when the NS is not rotating 
Mmax = data[4,:]            # Maximum mass for a given EOS
R = data[5,:]               # Radius 
Rratio = data[6,:]          # Ratio rp/re
Rstat = data[7,:]           # Radius when the NS is not rotating
Omega = data[8,:]*2*np.pi   # Angular velocity (rad/s), data[8,:] is spin frequency (in Hz)
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


def compute_surface(x, y, z):
    # Converting each column of data in each axis into arrays
    a1 = np.array(x)
    a2 = np.array(y)
    a3 = np.array(z)
    # Coverting the arrays into coordinate points in 3D
    datapoints = np.c_[a1,a2,a3]
    # Assigning the set of the first two coordinates to an array, X, and  
    # last coordinate to an array, Y
    X = datapoints[:,0:2]
    Y = datapoints[:,-1]
    # Degree of the polynomial fitting equation of the surface
    deg_of_poly = 4
    poly = PolynomialFeatures(degree=deg_of_poly)
    X_ = poly.fit_transform(X)
    # Fitting the linear model
    clf = linear_model.LinearRegression(fit_intercept=False)
    clf.fit(X_, Y)
    # Computing the surface
    N = 30
    Lengthx = 0.4
    Lengthy = 1.0
    predict_x0, predict_x1 = np.meshgrid(np.linspace(0, Lengthx, N), np.linspace(0, Lengthy, N))
    predict_x = np.concatenate((predict_x0.reshape(-1, 1), predict_x1.reshape(-1, 1)), axis=1)
    predict_x_ = poly.fit_transform(predict_x)
    predict_y = clf.predict(predict_x_)

    #print(poly.powers_)
    print(clf.coef_)
    print(poly.get_feature_names())
    rsquared = clf.score(X_,Y, sample_weight=None)
    print('R^2= '+repr(rsquared))
    
    return  predict_x0, predict_x, predict_y, predict_x1

print('Choose the best surface to generate (Choose only integer numbers)')
print('1. For (M0-M)/M0')
print('2. For (M-M*)/M*')
print('3. For (R-R*)/R*')
print('4. For M/Mmax')

# Choose the ssurface plot to make
method = input()
#method = 3

if filename[8] == 'e':
    fil = filename.replace('NS_data_eos', '')
    name = fil.replace('.txt', '')
else:
    name = 'All EOS'

if method == 1:
    x = delta**2
    y = (M/R)/Qmax  #y = G*mstat/(rstat*c**2)
    z = (M0-M)/M0
    predict_x0, predict_x, predict_y, predict_x1 = compute_surface(x, y, z)
    # Plotting surface and points
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(predict_x0, predict_x1, predict_y.reshape(predict_x0.shape), rstride=1, cstride=1, alpha=0.5)
    ax.scatter(x, y, z, c='b', marker='o', label=str(name))
    ax.set_xlabel(r'$\Omega^2 (R_*^3 / GM_*)$', fontsize='15')
    ax.set_ylabel(r'$(M/R)/(M_{max}/R_{max})$', fontsize='15')
    ax.set_zlabel(r'$(M_0-M)/M_0$', fontsize='15')
    plt.legend()
    #Rotation angle
    # ax.view_init(azim=140)
    #Limits in axes
    #ax.axes.set_xlim3d(left=0, right=100) 
    #ax.axes.set_ylim3d(bottom=0, top=2) 
    #ax.axes.set_zlim3d(bottom=0, top=100)
    plt.savefig('M0-frac-bestfitsurface.png', format = 'png', transparent=False)
    plt.show()

elif method == 2:
    x = delta**2             
    y = (M/R)/Qmax
    z = (M-Mstat)/Mstat
    predict_x0, predict_x, predict_y, predict_x1 = compute_surface(x, y, z)
    # Plotting surface and points
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(predict_x0, predict_x1, predict_y.reshape(predict_x0.shape), rstride=1, cstride=1, alpha=0.5)
    ax.scatter(x, y, z, c='b', marker='o', label=str(name))
    ax.set_xlabel(r'$\Omega^2 (R_*^3 / GM_*)$', fontsize='15')
    ax.set_ylabel(r'$(M/R)/(M_{max}/R_{max})$', fontsize='15')
    ax.set_zlabel(r'$(M-M_*)/M_*$', fontsize='15')
    plt.legend()
    #Rotation angle
    # ax.view_init(azim=140)
    #Limits in axes
    #ax.axes.set_xlim3d(left=0, right=100) 
    #ax.axes.set_ylim3d(bottom=0, top=2) 
    #ax.axes.set_zlim3d(bottom=0, top=100)
    plt.savefig('Mstat-frac-bestfitsurface.png', format = 'png', transparent=False)
    plt.show()

elif method == 3:
    x = delta**2
    y = (M/R)/Qmax 
    z = (R-Rstat)/Rstat
    predict_x0, predict_x, predict_y, predict_x1 = compute_surface(x, y, z)
    # Plotting surface and points
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(predict_x0, predict_x1, predict_y.reshape(predict_x0.shape), rstride=1, cstride=1, alpha=0.5)
    ax.scatter(x, y, z, c='b', marker='o', label=str(name))
    ax.set_xlabel(r'$\Omega^2 (R_*^3 / GM_*)$', fontsize='15')
    ax.set_ylabel(r'$(M/R)/(M_{max}/R_{max})$', fontsize='15')
    ax.set_zlabel(r'$(R-R_*)/R_*$', fontsize='15')
    plt.legend()
    #Rotation angle
    # ax.view_init(azim=140)
    #Limits in axes
    #ax.axes.set_xlim3d(left=0, right=100) 
    #ax.axes.set_ylim3d(bottom=0, top=2) 
    #ax.axes.set_zlim3d(bottom=0, top=100)
    plt.savefig('Rstat-frac-bestfitsurface.png', format = 'png', transparent=False)
    plt.show()

elif method == 4:
    x = delta**2
    y = (M/R)/Qmax
    z = M/Mmax
    predict_x0, predict_x, predict_y, predict_x1 = compute_surface(x, y, z)
    # Plotting surface and points
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(predict_x0, predict_x1, predict_y.reshape(predict_x0.shape), rstride=1, cstride=1, alpha=0.5)
    ax.scatter(x, y, z, c='b', marker='o', label=str(name))
    ax.set_xlabel(r'$\Omega^2 (R_*^3 / GM_*)$', fontsize='15')
    ax.set_ylabel(r'$(M/R)/(M_{max}/R_{max})$', fontsize='15')
    ax.set_zlabel(r'$M/M_{max}$', fontsize='15')
    plt.legend()
    #Rotation angle
    # ax.view_init(azim=140)
    #Limits in axes
    #ax.axes.set_xlim3d(left=0, right=100) 
    #ax.axes.set_ylim3d(bottom=0, top=2) 
    #ax.axes.set_zlim3d(bottom=0, top=100)
    plt.savefig('M-Mmax-bestfitsurface.png', format = 'png', transparent=False)
    plt.show()

else:
    print('Enter a number from 1 to 4')