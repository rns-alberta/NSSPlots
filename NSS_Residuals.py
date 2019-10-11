
# -*- coding: utf-8 -*-

from builtins import input
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from sklearn.preprocessing import PolynomialFeatures
from sklearn import linear_model
#import numpy.polynomial.polynomial as poly

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

##################################################

data_all = np.loadtxt('NS_data_allEOS.txt', unpack=True)

# Columns of data from the previous file 
Ec_all = data_all[0,:]*10**15       # central energy density
M_all= data_all[1,:]               # Total mass
M0_all = data_all[2,:]              # Baryonic mass
Mstat_all = data_all[3,:]           # Mass when the NS is not rotating 
Mmax_all = data_all[4,:]            # Maximum mass for a given EOS
R_all = data_all[5,:]               # Radius 
Rratio_all = data_all[6,:]          # Ratio rp/re
Rstat_all = data_all[7,:]           # Radius when the NS is not rotating
Omega_all = data_all[8,:]*2*np.pi   # Angular velocity (rad/s), data[8,:] is spin frequency (in Hz)
Klim_all = data_all[9,:]            # Kepler limit (Hz)
J_all = data_all[10,:]              # Angular momentum
T_all = data_all[11,:]              # Rotational kinetic energy
W_all = data_all[12,:]              # Gravitational binding energy
Rmax_all = data_all[13,:]           # Maximum radius for a given EOS
Qmax_all = data_all[14,:]           # Ratio (MaxMass / MaxRadius) of the non rotating star

# Converting quantities to CGS units
m_all = M_all * 1.9884e33   # Mass
r_all = R_all * 1.0e5       # Radius
rstat_all = Rstat_all * 1.0e5   # Radius of the non rotating NS
mstat_all = Mstat_all * 1.9884e33   # Mass of the non rotating NS
delta_all = Omega_all * ((rstat_all**3)/(G*mstat_all))**(0.5)   # Normalized omega
m0_all = M0_all * 1.9884e33     # Baryonic mass
normJ_all = (c*J_all)/(G*m0_all**2)     # Normalized angular momentum

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
    coefs = clf.coef_
    vari = poly.get_feature_names()
    
    return  predict_x0, predict_x, predict_y, predict_x1, coefs, vari



print('Choose the residual plot to generate (Choose only integer numbers)')
print('1. For (M0-M)/M0')
print('2. For (M-M*)/M*')
print('3. For (R-R*)/R*')

# Choose the surface plot to make
method = eval(input())
#method = 3

fil = filename.replace('NS_data_eos', '')
name = fil.replace('.txt', '')


if method == 1:
    x = delta**2
    y = G*mstat/(rstat*c**2)  
    z = (M0-M)/M0

    resid2 = np.zeros(len(x))

    predict_x0, predict_x, predict_y, predict_x1, coefs, vari = compute_surface(x, y, z)
      
    vari = [v.replace('^', '**') for v in vari]   
    vari = [v.replace(' ', '*') for v in vari]      

    a = 0
    for i in range(len(x)):
        globals()[vari[1]] = x[i]      
        globals()[vari[2]] = y[i]  
        
        for j in range(0, len(vari)):
            a = a + ( eval(vari[j]) * coefs[j] )

        resid2[i] = (z[i] - a)
        if j == len(vari)-1:
            a=0
        
    ax1.scatter(x, y, resid2, marker='o', label=str(name))
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

    resid2 = np.zeros(len(x))

    predict_x0, predict_x, predict_y, predict_x1, coefs, vari = compute_surface(x, y, z)
      
    vari = [v.replace('^', '**') for v in vari]   
    vari = [v.replace(' ', '*') for v in vari]      

    a = 0
    for i in range(len(x)):
        globals()[vari[1]] = x[i]      
        globals()[vari[2]] = y[i]  
        
        for j in range(0, len(vari)):
            a = a + ( eval(vari[j]) * coefs[j] )

        resid2[i] = (z[i] - a)
        if j == len(vari)-1:
            a=0
        
    ax1.scatter(x, y, resid2, marker='o', label=str(name))
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

    resid2 = np.zeros(len(x))

    predict_x0, predict_x, predict_y, predict_x1, coefs, vari = compute_surface(x, y, z)
      
    vari = [v.replace('^', '**') for v in vari]   
    vari = [v.replace(' ', '*') for v in vari]      

    a = 0
    for i in range(len(x)):
        globals()[vari[1]] = x[i]      
        globals()[vari[2]] = y[i]  
        
        for j in range(0, len(vari)):
            a = a + ( eval(vari[j]) * coefs[j] )

        resid2[i] = (z[i] - a)
        if j == len(vari)-1:
            a=0
    
    ax1.scatter(x, y, resid2, marker='o', label=str(name))
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
