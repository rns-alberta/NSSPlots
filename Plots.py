#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Code to compare various plots given files with different quantities 
of rotating neutron stars, considering ten EOS

"""

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure(figsize = (9,7))

# Files to read the data from
dataAB1 = np.loadtxt('NS_data_eosABPR1.txt', unpack=True)
dataAB2 = np.loadtxt('NS_data_eosABPR2.txt', unpack=True)
dataAB3 = np.loadtxt('NS_data_eosABPR3.txt', unpack=True)
dataAP = np.loadtxt('NS_data_eosAPR.txt', unpack=True)
dataB1 = np.loadtxt('NS_data_eosBBB1.txt', unpack=True)
dataB2 = np.loadtxt('NS_data_eosBBB2.txt', unpack=True)
dataHL1 = np.loadtxt('NS_data_eosHLPS1.txt', unpack=True)
dataHL2 = np.loadtxt('NS_data_eosHLPS2.txt', unpack=True)
dataHL3 = np.loadtxt('NS_data_eosHLPS3.txt', unpack=True)
dataL = np.loadtxt('NS_data_eosL.txt', unpack=True)

# Some constants
G = 6.67e-8
c = 2.99792458e10
Msun = 1

# Extracting columns from the data files of each EOS
# EOS ABPR1
ec_AB1 = dataAB1[0,:]*10**15   # central energy density (in g/cm^3)
M_AB1 = dataAB1[1,:]           # Total mass (in Msun)
M0_AB1 = dataAB1[2,:]          # Baryonic mass (in Msun)
Mstat_AB1 = dataAB1[3,:]       # Mass when the NS is not rotating (in Msun)
Mmax_AB1 = dataAB1[4,:]        # Maximum mass for a given EOS (in Msun)
R_AB1 = dataAB1[5,:]           # Radius (in km)
Rratio_AB1 = dataAB1[6,:]      # Ratio rp/re
Rstat_AB1 = dataAB1[7,:]       # Radius when the NS is not rotating (in km)
freq_AB1 = dataAB1[8,:]        # Rotational frequency (in Hz)
kfreq_AB1 = dataAB1[9,:]       # Kepler limit (in Hz)
J_AB1 = dataAB1[10,:]          # Angular momentum (in cm^2 g/s)
T_AB1 = dataAB1[11,:]          # Rotational kinetic energy (in g)
W_AB1 = dataAB1[12,:]          # Gravitational binding energy (in g)
Rmax_AB1 = dataAB1[13,:]       # Maximum radius of the non rotating star (in km)
Qmax_AB1 = dataAB1[14,:]       # MaxMass / MaxRadius of the non rotating star

# EOS ABPR2
ec_AB2 = dataAB2[0,:]*10**15   # central energy density (in g/cm^3)
M_AB2 = dataAB2[1,:]           # Total mass (in Msun)
M0_AB2 = dataAB2[2,:]          # Baryonic mass (in Msun)
Mstat_AB2 = dataAB2[3,:]       # Mass when the NS is not rotating (in Msun)
Mmax_AB2 = dataAB2[4,:]        # Maximum mass for a given EOS (in Msun)
R_AB2 = dataAB2[5,:]           # Radius (in km)
Rratio_AB2 = dataAB2[6,:]      # Ratio rp/re
Rstat_AB2 = dataAB2[7,:]       # Radius when the NS is not rotating (in km)
freq_AB2 = dataAB2[8,:]        # Rotational frequency (in Hz)
kfreq_AB2 = dataAB2[9,:]       # Kepler limit (in Hz)
J_AB2 = dataAB2[10,:]          # Angular momentum (in cm^2 g/s)
T_AB2 = dataAB2[11,:]          # Rotational kinetic energy (in g)
W_AB2 = dataAB2[12,:]          # Gravitational binding energy (in g)
Rmax_AB2 = dataAB2[13,:]       # Maximum radius of the non rotating star (in km)
Qmax_AB2 = dataAB2[14,:]       # MaxMass / MaxRadius of the non rotating star

# EOS ABPR3
ec_AB3 = dataAB3[0,:]*10**15   # central energy density (in g/cm^3)
M_AB3 = dataAB3[1,:]           # Total mass (in Msun)
M0_AB3 = dataAB3[2,:]          # Baryonic mass (in Msun)
Mstat_AB3 = dataAB3[3,:]       # Mass when the NS is not rotating (in Msun)
Mmax_AB3 = dataAB3[4,:]        # Maximum mass for a given EOS (in Msun)
R_AB3 = dataAB3[5,:]           # Radius (in km)
Rratio_AB3 = dataAB3[6,:]      # Ratio rp/re
Rstat_AB3 = dataAB3[7,:]       # Radius when the NS is not rotating (in km)
freq_AB3 = dataAB3[8,:]        # Rotational frequency (in Hz)
kfreq_AB3 = dataAB3[9,:]       # Kepler limit (in Hz)
J_AB3 = dataAB3[10,:]          # Angular momentum (in cm^2 g/s)
T_AB3 = dataAB3[11,:]          # Rotational kinetic energy (in g)
W_AB3 = dataAB3[12,:]          # Gravitational binding energy (in g)
Rmax_AB3 = dataAB3[13,:]       # Maximum radius of the non rotating star (in km)
Qmax_AB3 = dataAB3[14,:]       # MaxMass / MaxRadius of the non rotating star

# EOS APR
ec_AP = dataAP[0,:]*10**15   # central energy density (in g/cm^3)
M_AP = dataAP[1,:]           # Total mass (in Msun)
M0_AP = dataAP[2,:]          # Baryonic mass (in Msun)
Mstat_AP = dataAP[3,:]       # Mass when the NS is not rotating  (in Msun)
Mmax_AP = dataAP[4,:]        # Maximum mass for a given EOS (in Msun)
R_AP = dataAP[5,:]           # Radius (in km)
Rratio_AP = dataAP[6,:]      # Ratio rp/re
Rstat_AP = dataAP[7,:]       # Radius when the NS is not rotating (in km)
freq_AP = dataAP[8,:]        # Rotational frequency (in Hz)
kfreq_AP = dataAP[9,:]       # Kepler limit (in Hz)
J_AP = dataAP[10,:]          # Angular momentum (in cm^2 g/s)
T_AP = dataAP[11,:]          # Rotational kinetic energy (in g)
W_AP = dataAP[12,:]          # Gravitational binding energy (in g)
Rmax_AP = dataAP[13,:]       # Maximum radius of the non rotating star (in km)
Qmax_AP = dataAP[14,:]       # MaxMass / MaxRadius of the non rotating star

# EOS BBB1
ec_B1 = dataB1[0,:]*10**15   # central energy density (in g/cm^3)
M_B1 = dataB1[1,:]           # Total mass (in Msun)
M0_B1 = dataB1[2,:]          # Baryonic mass (in Msun)
Mstat_B1 = dataB1[3,:]       # Mass when the NS is not rotating  (in Msun)
Mmax_B1 = dataB1[4,:]        # Maximum mass for a given EOS (in Msun)
R_B1 = dataB1[5,:]           # Radius (in km)
Rratio_B1 = dataB1[6,:]      # Ratio rp/re
Rstat_B1 = dataB1[7,:]       # Radius when the NS is not rotating (in km)
freq_B1 = dataB1[8,:]        # Rotational frequency (in Hz)
kfreq_B1 = dataB1[9,:]       # Kepler limit (in Hz)
J_B1 = dataB1[10,:]          # Angular momentum (in cm^2 g/s)
T_B1 = dataB1[11,:]          # Rotational kinetic energy  (in g)
W_B1 = dataB1[12,:]          # Gravitational binding energy (in g)
Rmax_B1 = dataB1[13,:]       # Maximum radius of the non rotating star (in km)
Qmax_B1 = dataB1[14,:]       # MaxMass / MaxRadius of the non rotating star

# EOS BBB2
ec_B2 = dataB2[0,:]*10**15   # central energy density (in g/cm^3)
M_B2 = dataB2[1,:]           # Total mass (in Msun)
M0_B2 = dataB2[2,:]          # Baryonic mass (in Msun)
Mstat_B2 = dataB2[3,:]       # Mass when the NS is not rotating (in Msun)
Mmax_B2 = dataB2[4,:]        # Maximum mass for a given EOS (in Msun)
R_B2 = dataB2[5,:]           # Radius (in km)
Rratio_B2 = dataB2[6,:]      # Ratio rp/re
Rstat_B2 = dataB2[7,:]       # Radius when the NS is not rotating (in km)
freq_B2 = dataB2[8,:]        # Rotational frequency (in Hz)
kfreq_B2 = dataB2[9,:]       # Kepler limit (in Hz)
J_B2 = dataB2[10,:]          # Angular momentum (in cm^2 g/s)
T_B2 = dataB2[11,:]          # Rotational kinetic energy (in g)
W_B2 = dataB2[12,:]          # Gravitational binding energy (in g)
Rmax_B2 = dataB2[13,:]       # Maximum radius of the non rotating star (in km)
Qmax_B2 = dataB2[14,:]       # MaxMass / MaxRadius of the non rotating star

# EOS HLPS1
ec_HL1 = dataHL1[0,:]*10**15   # central energy density (in g/cm^3)
M_HL1 = dataHL1[1,:]           # Total mass (in Msun)
M0_HL1 = dataHL1[2,:]          # Baryonic mass (in Msun)
Mstat_HL1 = dataHL1[3,:]       # Mass when the NS is not rotating (in Msun)
Mmax_HL1 = dataHL1[4,:]        # Maximum mass for a given EOS (in Msun)
R_HL1 = dataHL1[5,:]           # Radius (in km)
Rratio_HL1 = dataHL1[6,:]      # Ratio rp/re
Rstat_HL1 = dataHL1[7,:]       # Radius when the NS is not rotating (in km)
freq_HL1 = dataHL1[8,:]        # Rotational frequency (in Hz)
kfreq_HL1 = dataHL1[9,:]       # Kepler limit (in Hz)
J_HL1 = dataHL1[10,:]          # Angular momentum (in cm^2 g/s)
T_HL1 = dataHL1[11,:]          # Rotational kinetic energy (in g)
W_HL1 = dataHL1[12,:]          # Gravitational binding energy (in g)
Rmax_HL1 = dataHL1[13,:]       # Maximum radius of the non rotating star (in km)
Qmax_HL1 = dataHL1[14,:]       # MaxMass / MaxRadius of the non rotating star

# EOS HLPS2
ec_HL2 = dataHL2[0,:]*10**15   # central energy density (in g/cm^3)
M_HL2 = dataHL2[1,:]           # Total mass (in Msun)
M0_HL2 = dataHL2[2,:]          # Baryonic mass (in Msun)
Mstat_HL2 = dataHL2[3,:]       # Mass when the NS is not rotating (in Msun)
Mmax_HL2 = dataHL2[4,:]        # Maximum mass for a given EOS (in Msun)
R_HL2 = dataHL2[5,:]           # Radius (in km)
Rratio_HL2 = dataHL2[6,:]      # Ratio rp/re
Rstat_HL2 = dataHL2[7,:]       # Radius when the NS is not rotating (in km)
freq_HL2 = dataHL2[8,:]        # Rotational frequency (in Hz)
kfreq_HL2 = dataHL2[9,:]       # Kepler limit (in Hz)
J_HL2 = dataHL2[10,:]          # Angular momentum (in cm^2 g/s)
T_HL2 = dataHL2[11,:]          # Rotational kinetic energy (in g)
W_HL2 = dataHL2[12,:]          # Gravitational binding energy (in g)
Rmax_HL2 = dataHL2[13,:]       # Maximum radius of the non rotating star (in km)
Qmax_HL2 = dataHL2[14,:]       # MaxMass / MaxRadius of the non rotating star

# EOS HLPS3
ec_HL3 = dataHL3[0,:]*10**15   # central energy density (in g/cm^3)
M_HL3 = dataHL3[1,:]           # Total mass (in Msun)
M0_HL3 = dataHL3[2,:]          # Baryonic mass (in Msun)
Mstat_HL3 = dataHL3[3,:]       # Mass when the NS is not rotating (in Msun)
Mmax_HL3 = dataHL3[4,:]        # Maximum mass for a given EOS (in Msun)
R_HL3 = dataHL3[5,:]           # Radius (in km)
Rratio_HL3 = dataHL3[6,:]      # Ratio rp/re
Rstat_HL3 = dataHL3[7,:]       # Radius when the NS is not rotating (in km)
freq_HL3 = dataHL3[8,:]        # Rotational frequency (in Hz)
kfreq_HL3 = dataHL3[9,:]       # Kepler limit (in Hz)
J_HL3 = dataHL3[10,:]          # Angular momentum (in cm^2 g/s)
T_HL3 = dataHL3[11,:]          # Rotational kinetic energy (in g)
W_HL3 = dataHL3[12,:]          # Gravitational binding energy (in g)
Rmax_HL3 = dataHL3[13,:]       # Maximum radius of the non rotating star (in km)
Qmax_HL3 = dataHL3[14,:]       # MaxMass / MaxRadius of the non rotating star

# EOS L
ec_L = dataL[0,:]*10**15   # central energy density (in g/cm^3)
M_L = dataL[1,:]           # Total mass (in Msun)
M0_L = dataL[2,:]          # Baryonic mass (in Msun)
Mstat_L = dataL[3,:]       # Mass when the NS is not rotating (in Msun)
Mmax_L = dataL[4,:]        # MAximum mass for a given EOS (in Msun)
R_L = dataL[5,:]           # Radius (in km)
Rratio_L = dataL[6,:]      # Ratio rp/re
Rstat_L = dataL[7,:]       # Radius when the NS is not rotating (in km)
freq_L = dataL[8,:]        # Rotational frequency (in Hz)
kfreq_L = dataL[9,:]       # Kepler limit (in Hz)
J_L = dataL[10,:]          # Angular momentum (in cm^2 g/s)
T_L = dataL[11,:]          # Rotational kinetic energy (in g)
W_L = dataL[12,:]          # Gravitational binding energy (in g)
Rmax_L = dataL[13,:]       # Maximum radius of the non rotating star (in km)
Qmax_L = dataL[14,:]       # MaxMass / MaxRadius of the non rotating star

# Converting parameters to CGS units
# EOS ABPR1
m_AB1 = M_AB1 * 1.9884e33   # Mass (in g)
r_AB1 = R_AB1 * 1.0e5       # Radius (in cm)
mmax_AB1 = Mmax_AB1 * 1.9884e33   # Maximum mass for a given EOS (in g)
rstat_AB1 = Rstat_AB1 * 1.0e5   # Radius of the non rotating NS (in cm)
mstat_AB1 = Mstat_AB1 * 1.9884e33   # Mass of the non rotating NS (in g)
delta_AB1 = freq_AB1 * (2*np.pi) * ((rstat_AB1**3)/(G*mstat_AB1))**(0.5)   # Normalized angular frequency
m0_AB1 = M0_AB1 * 1.9884e33     # Baryonic mass
normJ_AB1 = (c*J_AB1)/(G*m0_AB1**2)     # Normalized angular momentum
a_AB1 = (c*J_AB1)/(G*m_AB1**2)         # spin parameter
# EOS ABPR2
m_AB2 = M_AB2 * 1.9884e33   # Mass 
r_AB2 = R_AB2 * 1.0e5       # Radius 
mmax_AB2 = Mmax_AB2 * 1.9884e33   # Maximum mass for a given EOS (in g)
rstat_AB2 = Rstat_AB2 * 1.0e5   # Radius of the non rotating NS
mstat_AB2 = Mstat_AB2 * 1.9884e33   # Mass of the non rotating NS
delta_AB2 = freq_AB2 * (2*np.pi) * ((rstat_AB2**3)/(G*mstat_AB2))**(0.5)   # Normalized angular frequency
m0_AB2 = M0_AB2 * 1.9884e33     # Baryonic mass
normJ_AB2 = (c*J_AB2)/(G*m0_AB2**2)     # Normalized angular momentum
a_AB2 = (c*J_AB2)/(G*m_AB2**2)         # spin parameter
# EOS ABPR3
m_AB3 = M_AB3 * 1.9884e33   # Mass
r_AB3 = R_AB3 * 1.0e5       # Radius
mmax_AB3 = Mmax_AB3 * 1.9884e33   # Maximum mass for a given EOS (in g)
rstat_AB3 = Rstat_AB3 * 1.0e5   # Radius of the non rotating NS
mstat_AB3 = Mstat_AB3 * 1.9884e33   # Mass of the non rotating NS
delta_AB3 = freq_AB3 * (2*np.pi) * ((rstat_AB3**3)/(G*mstat_AB3))**(0.5)   # Normalized angular frequency
m0_AB3 = M0_AB3 * 1.9884e33     # Baryonic mass
normJ_AB3 = (c*J_AB3)/(G*m0_AB3**2)     # Normalized angular momentum
a_AB3 = (c*J_AB3)/(G*m_AB3**2)         # spin parameter
# EOS APR
m_AP = M_AP * 1.9884e33   # Mass
r_AP = R_AP * 1.0e5       # Radius
mmax_AP = Mmax_AP * 1.9884e33   # Maximum mass for a given EOS
rstat_AP = Rstat_AP * 1.0e5   # Radius of the non rotating NS
mstat_AP = Mstat_AP * 1.9884e33   # Mass of the non rotating NS
delta_AP = freq_AP * (2*np.pi) * ((rstat_AP**3)/(G*mstat_AP))**(0.5)   # Normalized angular frequency
m0_AP = M0_AP * 1.9884e33     # Baryonic mass
normJ_AP = (c*J_AP)/(G*m0_AP**2)     # Normalized angular momentum
a_AP = (c*J_AP)/(G*m_AP**2)         # spin parameter
# EOS BBB1
m_B1 = M_B1 * 1.9884e33   # Mass
r_B1 = R_B1 * 1.0e5       # Radius
mmax_B1 = Mmax_B1 * 1.9884e33   # Maximum mass for a given EOS
rstat_B1 = Rstat_B1 * 1.0e5   # Radius of the non rotating NS
mstat_B1 = Mstat_B1 * 1.9884e33   # Mass of the non rotating NS
delta_B1 = freq_B1 * (2*np.pi) * ((rstat_B1**3)/(G*mstat_B1))**(0.5)   # Normalized angular frequency
m0_B1 = M0_B1 * 1.9884e33     # Baryonic mass
normJ_B1 = (c*J_B1)/(G*m0_B1**2)     # Normalized angular momentum
a_B1 = (c*J_B1)/(G*m_B1**2)         # spin parameter
# EOS BBB2
m_B2 = M_B2 * 1.9884e33   # Mass
r_B2 = R_B2 * 1.0e5       # Radius
mmax_B2 = Mmax_B2 * 1.9884e33   # Maximum mass for a given EOS
rstat_B2 = Rstat_B2 * 1.0e5   # Radius of the non rotating NS
mstat_B2 = Mstat_B2 * 1.9884e33   # Mass of the non rotating NS
delta_B2 = freq_B2 * (2*np.pi) * ((rstat_B2**3)/(G*mstat_B2))**(0.5)   # Normalized angular frequency
m0_B2 = M0_B2 * 1.9884e33     # Baryonic mass
normJ_B2 = (c*J_B2)/(G*m0_B2**2)     # Normalized angular momentum
a_B2 = (c*J_B2)/(G*m_B2**2)         # spin parameter
# EOS HLPS1
m_HL1 = M_HL1 * 1.9884e33   # Mass
r_HL1 = R_HL1 * 1.0e5       # Radius
mmax_HL1 = Mmax_HL1 * 1.9884e33   # Maximum mass for a given EOS
rstat_HL1 = Rstat_HL1 * 1.0e5   # Radius of the non rotating NS
mstat_HL1 = Mstat_HL1 * 1.9884e33   # Mass of the non rotating NS
delta_HL1 = freq_HL1 * (2*np.pi) * ((rstat_HL1**3)/(G*mstat_HL1))**(0.5)   # Normalized angular frequency
m0_HL1 = M0_HL1 * 1.9884e33     # Baryonic mass
normJ_HL1 = (c*J_HL1)/(G*m0_HL1**2)     # Normalized angular momentum
a_HL1 = (c*J_HL1)/(G*m_HL1**2)         # spin parameter
# EOS HLPS2
m_HL2 = M_HL2 * 1.9884e33   # Mass
r_HL2 = R_HL2 * 1.0e5       # Radius
mmax_HL2 = Mmax_HL2 * 1.9884e33   # Maximum mass for a given EOS
rstat_HL2 = Rstat_HL2 * 1.0e5   # Radius of the non rotating NS
mstat_HL2 = Mstat_HL2 * 1.9884e33   # Mass of the non rotating NS
delta_HL2 = freq_HL2 * (2*np.pi) * ((rstat_HL2**3)/(G*mstat_HL2))**(0.5)   # Normalized angular frequency
m0_HL2 = M0_HL2 * 1.9884e33     # Baryonic mass
normJ_HL2 = (c*J_HL2)/(G*m0_HL2**2)     # Normalized angular momentum
a_HL2 = (c*J_HL2)/(G*m_HL2**2)         # spin parameter
# EOS HLPS3
m_HL3 = M_HL3 * 1.9884e33   # Mass
r_HL3 = R_HL3 * 1.0e5       # Radius
mmax_HL3 = Mmax_HL3 * 1.9884e33   # Maximum mass for a given EOS
rstat_HL3 = Rstat_HL3 * 1.0e5   # Radius of the non rotating NS
mstat_HL3 = Mstat_HL3 * 1.9884e33   # Mass of the non rotating NS
delta_HL3 = freq_HL3 * (2*np.pi) * ((rstat_HL3**3)/(G*mstat_HL3))**(0.5)   # Normalized angular frequency
m0_HL3 = M0_HL3 * 1.9884e33     # Baryonic mass
normJ_HL3 = (c*J_HL3)/(G*m0_HL3**2)     # Normalized angular momentum
a_HL3 = (c*J_HL3)/(G*m_HL3**2)         # spin parameter
# EOS L
m_L = M_L * 1.9884e33   # Mass
r_L = R_L * 1.0e5       # Radius
mmax_L = Mmax_L * 1.9884e33   # Maximum mass for a given EOS
rstat_L = Rstat_L * 1.0e5   # Radius of the non rotating NS
mstat_L = Mstat_L * 1.9884e33   # Mass of the non rotating NS
delta_L = freq_L * (2*np.pi) * ((rstat_L**3)/(G*mstat_L))**(0.5)   # Normalized angular frequency
m0_L = M0_L * 1.9884e33     # Baryonic mass
normJ_L = (c*J_L)/(G*m0_L**2)     # Normalized angular momentum
a_L = (c*J_L)/(G*m_L**2)         # spin parameter

# Menu
print('Choose the plot to generate (Choose only integer numbers)')
print('1. Sequences of Mass vs energy density')
print('2. Mass vs Radius')
print('3. Frequency vs Mass')
print('4. (M-M*)/M* vs frequency')
print('5. (M-M*)/M* vs normalized omega')
print('6. (M-M_0)/M_0 vs frequency')
print('7. (M-M_0)/M_0 vs normalized omega')
print('8. M/Mmax vs normalized omega')
print('9. (R-R*)/R* vs frequency')
print('10. (R-R*)/R* vs normalized omega')
print('11. Compactness vs frequency')
print('12. Compactness vs normalized omega')
print('13. Frequency vs Normalized angular momentum')
print('14. Frequency vs T/W')
print('15. 3D plot of frequency, normalized omega and normalized compactness ')
print('16. 3D plot of spin parameter, normalized omega and M/Mmax')
print('---------------------------------------------------------------')
print('17. 3D plot of M/Mmax, normalized omega and normal compactness')
print('18. 3D plot of (R-R*)/R*, normalized omega and normal compactness')
print('19. 3D plot of (M-M*)/M*, normalized omega and normal compactness')
print('20. 3D plot of (M-M_0)/M_0, normalized omega and normal compactness')
print('---------------------------------------------------------------')
print('21. 3D plot of M/Mmax, normalized omega and normalized compactness')
print('22. 3D plot of (R-R*)/R*, normalized omega and normalized compactness')
print('23. 3D plot of (M-M*)/M*, normalized omega and normalized compactness')
print('24. 3D plot of (M-M_0)/M_0, normalized omega and normalized compactness')
print('---------------------------------------------------------------')
print('25. 3D plot of M/Mmax, normalized omega and max normalized compactness')
print('26. 3D plot of (R-R*)/R*, normalized omega and max normalized compactness')
print('27. 3D plot of (M-M*)/M*, normalized omega and max normalized compactness')
print('28. 3D plot of (M-M_0)/M_0, normalized omega and max normalized compactness')
print('---------------------------------------------------------------')

# Choose the number of the desired plot 
method = input()
#method = 17

if method == 1:
    plt.scatter(ec_AB1, M_AB1, s=3, label = 'EOS ABPR1')
    plt.scatter(ec_AB2, M_AB2, s=3, label = 'EOS ABPR2')
    plt.scatter(ec_AB3, M_AB3, s=3, label = 'EOS ABPR3')
    plt.scatter(ec_AP, M_AP, s=3, label = 'EOS APR')
    plt.scatter(ec_B1, M_B1, s=3, label = 'EOS BBB1')
    plt.scatter(ec_B2, M_B2, s=3, label = 'EOS BBB2')
    plt.scatter(ec_HL1, M_HL1, s=3, label = 'EOS HLPS1')
    plt.scatter(ec_HL2, M_HL2, s=3, label = 'EOS HLPS2')
    plt.scatter(ec_HL3, M_HL3, s=3, label = 'EOS HLPS3')
    plt.scatter(ec_L, M_L, s=3, label = 'EOS L')
    plt.xlabel(r'$\varepsilon_c$ [g cm$^{-3}$]', fontsize='15')
    plt.ylabel(r'$M$ [M$_\odot$]', fontsize='15')
    plt.legend()
    plt.show()

elif method == 2:
    plt.scatter(R_AB1, M_AB1, s=3, label = 'EOS ABPR1')
    plt.scatter(R_AB2, M_AB2, s=3, label = 'EOS ABPR2')
    plt.scatter(R_AB3, M_AB3, s=3, label = 'EOS ABPR3')
    plt.scatter(R_AP, M_AP, s=3, label = 'EOS APR')
    plt.scatter(R_B1, M_B1, s=3, label = 'EOS BBB1')
    plt.scatter(R_B2, M_B2, s=3, label = 'EOS BBB2')
    plt.scatter(R_HL1, M_HL1, s=3, label = 'EOS HLPS1')
    plt.scatter(R_HL2, M_HL2, s=3, label = 'EOS HLPS2')
    plt.scatter(R_HL3, M_HL3, s=3, label = 'EOS HLPS3')
    plt.scatter(R_L, M_L, s=3, label = 'EOS L')
    plt.xlabel(r'$R$ [km]', fontsize='15')
    plt.ylabel(r'$M$ [M$_\odot$]', fontsize='15')
    plt.legend()
    plt.show()
    
elif method == 3:
    plt.scatter(M_AB1, freq_AB1, s=3, label = 'EOS ABPR1')
    plt.scatter(M_AB2, freq_AB2, s=3, label = 'EOS ABPR2')
    plt.scatter(M_AB3, freq_AB3, s=3, label = 'EOS ABPR3')
    plt.scatter(M_AP, freq_AP, s=3, label = 'EOS APR')
    plt.scatter(M_B1, freq_B1, s=3, label = 'EOS BBB1')
    plt.scatter(M_B2, freq_B2, s=3, label = 'EOS BBB2')
    plt.scatter(M_HL1, freq_HL1, s=3, label = 'EOS HLPS1')
    plt.scatter(M_HL2, freq_HL2, s=3, label = 'EOS HLPS2')
    plt.scatter(M_HL3, freq_HL3, s=3, label = 'EOS HLPS3')
    plt.scatter(M_L, freq_L, s=3, label = 'EOS L')
    plt.xlabel(r'$M$ [M$_\odot$]', fontsize='15')
    plt.ylabel(r'$\nu$ [Hz]', fontsize='15')
    plt.legend()
    plt.show()
    
elif method == 4:
    plt.scatter(freq_AB1, (M_AB1-Mstat_AB1)/Mstat_AB1, s=3, label = 'EOS ABPR1')
    plt.scatter(freq_AB2, (M_AB2-Mstat_AB2)/Mstat_AB2, s=3, label = 'EOS ABPR2')
    plt.scatter(freq_AB3, (M_AB3-Mstat_AB3)/Mstat_AB3, s=3, label = 'EOS ABPR3')
    plt.scatter(freq_AP, (M_AP-Mstat_AP)/Mstat_AP, s=3, label = 'EOS APR')
    plt.scatter(freq_B1, (M_B1-Mstat_B1)/Mstat_B1, s=3, label = 'EOS BBB1')
    plt.scatter(freq_B2, (M_B2-Mstat_B2)/Mstat_B2, s=3, label = 'EOS BBB2')
    plt.scatter(freq_HL1, (M_HL1-Mstat_HL1)/Mstat_HL1, s=3, label = 'EOS HLPS1')
    plt.scatter(freq_HL2, (M_HL2-Mstat_HL2)/Mstat_HL2, s=3, label = 'EOS HLPS2')
    plt.scatter(freq_HL3, (M_HL3-Mstat_HL3)/Mstat_HL3, s=3, label = 'EOS HLPS3')
    plt.scatter(freq_L, (M_L-Mstat_L)/Mstat_L, s=3, label = 'EOS L')
    plt.xlabel(r'$\nu$ [Hz]', fontsize='15')
    plt.ylabel(r'$(M-M_*)/M_*$', fontsize='15')
    plt.legend()
    plt.show()
    
elif method == 5:
    plt.scatter(delta_AB1, (M_AB1-Mstat_AB1)/Mstat_AB1, s=3, label = 'EOS ABPR1')
    plt.scatter(delta_AB2, (M_AB2-Mstat_AB2)/Mstat_AB2, s=3, label = 'EOS ABPR2')
    plt.scatter(delta_AB3, (M_AB3-Mstat_AB3)/Mstat_AB3, s=3, label = 'EOS ABPR3')
    plt.scatter(delta_AP, (M_AP-Mstat_AP)/Mstat_AP, s=3, label = 'EOS APR')
    plt.scatter(delta_B1, (M_B1-Mstat_B1)/Mstat_B1, s=3, label = 'EOS BBB1')
    plt.scatter(delta_B2, (M_B2-Mstat_B2)/Mstat_B2, s=3, label = 'EOS BBB2')
    plt.scatter(delta_HL1, (M_HL1-Mstat_HL1)/Mstat_HL1, s=3, label = 'EOS HLPS1')
    plt.scatter(delta_HL2, (M_HL2-Mstat_HL2)/Mstat_HL2, s=3, label = 'EOS HLPS2')
    plt.scatter(delta_HL3, (M_HL3-Mstat_HL3)/Mstat_HL3, s=3, label = 'EOS HLPS3')
    plt.scatter(delta_L, (M_L-Mstat_L)/Mstat_L, s=3, label = 'EOS L')
    plt.xlabel(r'$\Omega (R_*^3 / GM_*)^{1/2}$', fontsize='15')
    plt.ylabel(r'$(M-M_*)/M_*$', fontsize='15')
    plt.legend()
    plt.show()
    
elif method == 6:
    plt.scatter(freq_AB1, (M0_AB1-M_AB1)/M0_AB1, s=3, label = 'EOS ABPR1')
    plt.scatter(freq_AB2, (M0_AB2-M_AB2)/M0_AB2, s=3, label = 'EOS ABPR2')
    plt.scatter(freq_AB3, (M0_AB3-M_AB3)/M0_AB3, s=3, label = 'EOS ABPR3')
    plt.scatter(freq_AP, (M0_AP-M_AP)/M0_AP, s=3, label = 'EOS APR')
    plt.scatter(freq_B1, (M0_B1-M_B1)/M0_B1, s=3, label = 'EOS BBB1')
    plt.scatter(freq_B2, (M0_B2-M_B2)/M0_B2, s=3, label = 'EOS BBB2')
    plt.scatter(freq_HL1, (M0_HL1-M_HL1)/M0_HL1, s=3, label = 'EOS HLPS1')
    plt.scatter(freq_HL2, (M0_HL2-M_HL2)/M0_HL2, s=3, label = 'EOS HLPS2')
    plt.scatter(freq_HL3, (M0_HL3-M_HL3)/M0_HL3, s=3, label = 'EOS HLPS3')
    plt.scatter(freq_L, (M0_L-M_L)/M0_L, s=3, label = 'EOS L')
    plt.xlabel(r'$\nu$ [Hz]', fontsize='15')
    plt.ylabel(r'$(M_0-M)/M_0$', fontsize='15')
    plt.legend()
    plt.show()
    
elif method == 7:
    plt.scatter(delta_AB1, (M0_AB1-M_AB1)/M0_AB1, s=3, label = 'EOS ABPR1')
    plt.scatter(delta_AB2, (M0_AB2-M_AB2)/M0_AB2, s=3, label = 'EOS ABPR2')
    plt.scatter(delta_AB3, (M0_AB3-M_AB3)/M0_AB3, s=3, label = 'EOS ABPR3')
    plt.scatter(delta_AP, (M0_AP-M_AP)/M0_AP, s=3, label = 'EOS APR')
    plt.scatter(delta_B1, (M0_B1-M_B1)/M0_B1, s=3, label = 'EOS BBB1')
    plt.scatter(delta_B2, (M0_B2-M_B2)/M0_B2, s=3, label = 'EOS BBB2')
    plt.scatter(delta_HL1, (M0_HL1-M_HL1)/M0_HL1, s=3, label = 'EOS HLPS1')
    plt.scatter(delta_HL2, (M0_HL2-M_HL2)/M0_HL2, s=3, label = 'EOS HLPS2')
    plt.scatter(delta_HL3, (M0_HL3-M_HL3)/M0_HL3, s=3, label = 'EOS HLPS3')
    plt.scatter(delta_L, (M0_L-M_L)/M0_L, s=3, label = 'EOS L')
    plt.xlabel(r'$\Omega (R_*^3 / GM_*)^{1/2}$', fontsize='15')
    plt.ylabel(r'$(M_0-M)/M_0$', fontsize='15')
    plt.legend()
    plt.show()
  
elif method == 8:
    plt.scatter(delta_AB1, M_AB1/Mmax_AB1, s=3, label = 'EOS ABPR1')
    plt.scatter(delta_AB2, M_AB2/Mmax_AB2, s=3, label = 'EOS ABPR2')
    plt.scatter(delta_AB3, M_AB3/Mmax_AB3, s=3, label = 'EOS ABPR3')
    plt.scatter(delta_AP, M_AP/Mmax_AP, s=3, label = 'EOS APR')
    plt.scatter(delta_B1, M_B1/Mmax_B1, s=3, label = 'EOS BBB1')
    plt.scatter(delta_B2, M_B2/Mmax_B2, s=3, label = 'EOS BBB2')
    plt.scatter(delta_HL1, M_HL1/Mmax_HL1, s=3, label = 'EOS HLPS1')
    plt.scatter(delta_HL2, M_HL2/Mmax_HL2, s=3, label = 'EOS HLPS2')
    plt.scatter(delta_HL3, M_HL3/Mmax_HL3, s=3, label = 'EOS HLPS3')
    plt.scatter(delta_L, M_L/Mmax_L, s=3, label = 'EOS L')
    plt.xlabel(r'$\Omega (R_*^3 / GM_*)^{1/2}$', fontsize='15')
    plt.ylabel(r'$M/M_{max}$', fontsize='15')
    plt.legend()
    plt.show()
    
elif method == 9:
    plt.scatter(freq_AB1, (R_AB1-Rstat_AB1)/Rstat_AB1, s=3, label = 'EOS ABPR1')
    plt.scatter(freq_AB2, (R_AB2-Rstat_AB2)/Rstat_AB2, s=3, label = 'EOS ABPR2')
    plt.scatter(freq_AB3, (R_AB3-Rstat_AB3)/Rstat_AB3, s=3, label = 'EOS ABPR3')
    plt.scatter(freq_AP, (R_AP-Rstat_AP)/Rstat_AP, s=3, label = 'EOS APR')
    plt.scatter(freq_B1, (R_B1-Rstat_B1)/Rstat_B1, s=3, label = 'EOS BBB1')
    plt.scatter(freq_B2, (R_B2-Rstat_B2)/Rstat_B2, s=3, label = 'EOS BBB2')
    plt.scatter(freq_HL1, (R_HL1-Rstat_HL1)/Rstat_HL1, s=3, label = 'EOS HLPS1')
    plt.scatter(freq_HL2, (R_HL2-Rstat_HL2)/Rstat_HL2, s=3, label = 'EOS HLPS2')
    plt.scatter(freq_HL3, (R_HL3-Rstat_HL3)/Rstat_HL3, s=3, label = 'EOS HLPS3')
    plt.scatter(freq_L, (R_L-Rstat_L)/Rstat_L, s=3, label = 'EOS L')
    plt.xlabel(r'$\nu$ [Hz]', fontsize='15')
    plt.ylabel(r'$(R-R_*)/R_*$', fontsize='15')
    plt.legend()
    plt.show()

elif method == 10:
    plt.scatter(delta_AB1, (R_AB1-Rstat_AB1)/Rstat_AB1, s=3, label = 'EOS ABPR1')
    plt.scatter(delta_AB2, (R_AB2-Rstat_AB2)/Rstat_AB2, s=3, label = 'EOS ABPR2')
    plt.scatter(delta_AB3, (R_AB3-Rstat_AB3)/Rstat_AB3, s=3, label = 'EOS ABPR3')
    plt.scatter(delta_AP, (R_AP-Rstat_AP)/Rstat_AP, s=3, label = 'EOS APR')
    plt.scatter(delta_B1, (R_B1-Rstat_B1)/Rstat_B1, s=3, label = 'EOS BBB1')
    plt.scatter(delta_B2, (R_B2-Rstat_B2)/Rstat_B2, s=3, label = 'EOS BBB2')
    plt.scatter(delta_HL1, (R_HL1-Rstat_HL1)/Rstat_HL1, s=3, label = 'EOS HLPS1')
    plt.scatter(delta_HL2, (R_HL2-Rstat_HL2)/Rstat_HL2, s=3, label = 'EOS HLPS2')
    plt.scatter(delta_HL3, (R_HL3-Rstat_HL3)/Rstat_HL3, s=3, label = 'EOS HLPS3')
    plt.scatter(delta_L, (R_L-Rstat_L)/Rstat_L, s=3, label = 'EOS L')
    plt.xlabel(r'$\Omega (R_*^3 / GM_*)^{1/2}$', fontsize='15')
    plt.ylabel(r'$(R-R_*)/R_*$', fontsize='15')
    plt.legend()
    plt.show()
    
elif method == 11:
    plt.scatter(freq_AB1, (M_AB1/Msun)*(1.47/R_AB1), s=3, label = 'EOS ABPR1')
    plt.scatter(freq_AB2, (M_AB2/Msun)*(1.47/R_AB2), s=3, label = 'EOS ABPR2')
    plt.scatter(freq_AB3, (M_AB3/Msun)*(1.47/R_AB3), s=3, label = 'EOS ABPR3')
    plt.scatter(freq_AP, (M_AP/Msun)*(1.47/R_AP), s=3, label = 'EOS APR')
    plt.scatter(freq_B1, (M_B1/Msun)*(1.47/R_B1), s=3, label = 'EOS BBB1')
    plt.scatter(freq_B2, (M_B2/Msun)*(1.47/R_B2), s=3, label = 'EOS BBB2')
    plt.scatter(freq_HL1, (M_HL1/Msun)*(1.47/R_HL1), s=3, label = 'EOS HLPS1')
    plt.scatter(freq_HL2, (M_HL2/Msun)*(1.47/R_HL2), s=3, label = 'EOS HLPS2')
    plt.scatter(freq_HL3, (M_HL3/Msun)*(1.47/R_HL3), s=3, label = 'EOS HLPS3')
    plt.scatter(freq_L, (M_L/Msun)*(1.47/R_L), s=3, label = 'EOS L')
    plt.xlabel(r'$\nu$ [Hz]', fontsize='15')
    plt.ylabel(r'$\zeta$', fontsize='15')
    plt.legend()
    plt.show()
    
elif method == 12:
    plt.scatter(delta_AB1, (M_AB1/Msun)*(1.47/R_AB1), s=3, label = 'EOS ABPR1')
    plt.scatter(delta_AB2, (M_AB2/Msun)*(1.47/R_AB2), s=3, label = 'EOS ABPR2')
    plt.scatter(delta_AB3, (M_AB3/Msun)*(1.47/R_AB3), s=3, label = 'EOS ABPR3')
    plt.scatter(delta_AP, (M_AP/Msun)*(1.47/R_AP), s=3, label = 'EOS APR')
    plt.scatter(delta_B1, (M_B1/Msun)*(1.47/R_B1), s=3, label = 'EOS BBB1')
    plt.scatter(delta_B2, (M_B2/Msun)*(1.47/R_B2), s=3, label = 'EOS BBB2')
    plt.scatter(delta_HL1, (M_HL1/Msun)*(1.47/R_HL1), s=3, label = 'EOS HLPS1')
    plt.scatter(delta_HL2, (M_HL2/Msun)*(1.47/R_HL2), s=3, label = 'EOS HLPS2')
    plt.scatter(delta_HL3, (M_HL3/Msun)*(1.47/R_HL3), s=3, label = 'EOS HLPS3')
    plt.scatter(delta_L, (M_L/Msun)*(1.47/R_L), s=3, label = 'EOS L')
    plt.xlabel(r'$\Omega (R_*^3 / GM_*)^{1/2}$', fontsize='15')
    plt.ylabel(r'$\zeta$', fontsize='15')
    plt.legend()
    plt.show()
    
elif method == 13:
    plt.scatter(normJ_AB1, freq_AB1, s=3, label = 'EOS ABPR1')
    plt.scatter(normJ_AB2, freq_AB2, s=3, label = 'EOS ABPR2')
    plt.scatter(normJ_AB3, freq_AB3, s=3, label = 'EOS ABPR3')
    plt.scatter(normJ_AP, freq_AP, s=3, label = 'EOS APR')
    plt.scatter(normJ_B1, freq_B1, s=3, label = 'EOS BBB1')
    plt.scatter(normJ_B2, freq_B2, s=3, label = 'EOS BBB2')
    plt.scatter(normJ_HL1, freq_HL1, s=3, label = 'EOS HLPS1')
    plt.scatter(normJ_HL2, freq_HL2, s=3, label = 'EOS HLPS2')
    plt.scatter(normJ_HL3, freq_HL3, s=3, label = 'EOS HLPS3')
    plt.scatter(normJ_L, freq_L, s=3, label = 'EOS L')
    plt.xlabel(r'$cJ/GM_0^2$', fontsize='15')
    plt.ylabel(r'$\nu$ [Hz]', fontsize='15')
    plt.legend()
    plt.show()
    
elif method == 14:
    plt.scatter(T_AB1/W_AB1, freq_AB1, s=3, label = 'EOS ABPR1')
    plt.scatter(T_AB2/W_AB2, freq_AB2, s=3, label = 'EOS ABPR2')
    plt.scatter(T_AB3/W_AB3, freq_AB3, s=3, label = 'EOS ABPR3')
    plt.scatter(T_AP/W_AP, freq_AP, s=3, label = 'EOS APR')
    plt.scatter(T_B1/W_B1, freq_B1, s=3, label = 'EOS BBB1')
    plt.scatter(T_B2/W_B2, freq_B2, s=3, label = 'EOS BBB2')
    plt.scatter(T_HL1/W_HL1, freq_HL1, s=3, label = 'EOS HLPS1')
    plt.scatter(T_HL2/W_HL2, freq_HL2, s=3, label = 'EOS HLPS2')
    plt.scatter(T_HL3/W_HL3, freq_HL3, s=3, label = 'EOS HLPS3')
    plt.scatter(T_L/W_L, freq_L, s=3, label = 'EOS L')
    plt.xlabel(r'$T/W$', fontsize='15')
    plt.ylabel(r'$\nu$ [Hz]', fontsize='15')
    plt.legend()
    plt.show()

elif method == 15:
    ax = fig.gca(projection='3d')
    ax.scatter(delta_AB1**2, (M_AB1/R_AB1)/Qmax_AB1, freq_AB1, label = 'EOS ABPR1')
    ax.scatter(delta_AB2**2, (M_AB2/R_AB2)/Qmax_AB2, freq_AB2, label = 'EOS ABPR2')
    ax.scatter(delta_AB3**2, (M_AB3/R_AB3)/Qmax_AB3, freq_AB3, label = 'EOS ABPR3')
    ax.scatter(delta_AP**2, (M_AP/R_AP)/Qmax_AP, freq_AP, label = 'EOS APR')
    ax.scatter(delta_B1**2, (M_B1/R_B1)/Qmax_B1, freq_B1, label = 'EOS BBB1')
    ax.scatter(delta_B2**2, (M_B2/R_B2)/Qmax_B2, freq_B2, label = 'EOS BBB2')
    ax.scatter(delta_HL1**2, (M_HL1/R_HL1)/Qmax_HL1, freq_HL1, label = 'EOS HLPS1')
    ax.scatter(delta_HL2**2, (M_HL2/R_HL2)/Qmax_HL2, freq_HL2, label = 'EOS HLPS2')
    ax.scatter(delta_HL3**2, (M_HL3/R_HL3)/Qmax_HL3, freq_HL3, label = 'EOS HLPS3')
    ax.scatter(delta_L**2, (M_L/R_L)/Qmax_L, freq_L, label = 'EOS L')
    ax.set_xlabel(r'$\Omega^2 (R_*^3 / GM_*)$', fontsize='15')
    ax.set_ylabel(r'$(M/R)/(M_{max}/R_{max})$', fontsize='15')
    ax.set_zlabel(r'$\nu$ [Hz]', fontsize='15')
    plt.gca().invert_xaxis()
    ax.view_init(azim=110)
    plt.legend()
    plt.show()

elif method == 16:
    ax = fig.gca(projection='3d')
    ax.scatter(delta_AB1**2, M_AB1/Mmax_AB1, a_AB1, label = 'EOS ABPR1')
    ax.scatter(delta_AB2**2, M_AB2/Mmax_AB2, a_AB2, label = 'EOS ABPR2')
    ax.scatter(delta_AB3**2, M_AB3/Mmax_AB3, a_AB3, label = 'EOS ABPR3')
    ax.scatter(delta_AP**2, M_AP/Mmax_AP, a_AP, label = 'EOS APR')
    ax.scatter(delta_B1**2, M_B1/Mmax_B1, a_B1, label = 'EOS BBB1')
    ax.scatter(delta_B2**2, M_B2/Mmax_B2, a_B2, label = 'EOS BBB2')
    ax.scatter(delta_HL1**2, M_HL1/Mmax_HL1, a_HL1, label = 'EOS HLPS1')
    ax.scatter(delta_HL2**2, M_HL2/Mmax_HL2, a_HL2, label = 'EOS HLPS2')
    ax.scatter(delta_HL3**2, M_HL3/Mmax_HL3, a_HL3, label = 'EOS HLPS3')
    ax.scatter(delta_L**2, M_L/Mmax_L, a_L, label = 'EOS L')
    ax.set_xlabel(r'$\Omega (R_*^3 / GM_*)^{1/2}$', fontsize='15')
    ax.set_ylabel(r'$M/M_{max}$', fontsize='15')
    ax.set_zlabel(r'$a$', fontsize='15')
    #ax.view_init(azim=90)
    plt.gca().invert_xaxis()
    plt.legend()
    plt.show()
    
elif method == 17:
    ax = fig.gca(projection='3d')
    ax.scatter(delta_AB1**2, G*mstat_AB1/(rstat_AB1*c**2), M_AB1/Mmax_AB1, label = 'EOS ABPR1')
    ax.scatter(delta_AB2**2, G*mstat_AB2/(rstat_AB2*c**2), M_AB2/Mmax_AB2, label = 'EOS ABPR2')
    ax.scatter(delta_AB3**2, G*mstat_AB3/(rstat_AB3*c**2), M_AB3/Mmax_AB3, label = 'EOS ABPR3')
    ax.scatter(delta_AP**2, G*mstat_AP/(rstat_AP*c**2), M_AP/Mmax_AP, label = 'EOS APR')
    ax.scatter(delta_B1**2, G*mstat_B1/(rstat_B1*c**2), M_B1/Mmax_B1, label = 'EOS BBB1')
    ax.scatter(delta_B2**2, G*mstat_B2/(rstat_B2*c**2), M_B2/Mmax_B2, label = 'EOS BBB2')
    ax.scatter(delta_HL1**2, G*mstat_HL1/(rstat_HL1*c**2), M_HL1/Mmax_HL1, label = 'EOS HLPS1')
    ax.scatter(delta_HL2**2, G*mstat_HL2/(rstat_HL2*c**2), M_HL2/Mmax_HL2, label = 'EOS HLPS2')
    ax.scatter(delta_HL3**2, G*mstat_HL3/(rstat_HL3*c**2), M_HL3/Mmax_HL3, label = 'EOS HLPS3')
    ax.scatter(delta_L**2, G*mstat_L/(rstat_L*c**2), M_L/Mmax_L, label = 'EOS L')
    ax.set_xlabel(r'$\Omega^2 (R_*^3 / GM_*)$', fontsize='15')
    ax.set_ylabel(r'$\zeta$', fontsize='15')
    ax.set_zlabel(r'$M/M_{max}$', fontsize='15')
    plt.legend()
    plt.show()
    
elif method == 18:
    ax = fig.gca(projection='3d')
    ax.scatter(delta_AB1**2, G*mstat_AB1/(rstat_AB1*c**2), (R_AB1-Rstat_AB1)/Rstat_AB1, label = 'EOS ABPR1')
    ax.scatter(delta_AB2**2, G*mstat_AB2/(rstat_AB2*c**2), (R_AB2-Rstat_AB2)/Rstat_AB2, label = 'EOS ABPR2')
    ax.scatter(delta_AB3**2, G*mstat_AB3/(rstat_AB3*c**2), (R_AB3-Rstat_AB3)/Rstat_AB3, label = 'EOS ABPR3')
    ax.scatter(delta_AP**2, G*mstat_AP/(rstat_AP*c**2), (R_AP-Rstat_AP)/Rstat_AP, label = 'EOS APR')
    ax.scatter(delta_B1**2, G*mstat_B1/(rstat_B1*c**2), (R_B1-Rstat_B1)/Rstat_B1, label = 'EOS BBB1')
    ax.scatter(delta_B2**2, G*mstat_B2/(rstat_B2*c**2), (R_B2-Rstat_B2)/Rstat_B2, label = 'EOS BBB2')
    ax.scatter(delta_HL1**2, G*mstat_HL1/(rstat_HL1*c**2), (R_HL1-Rstat_HL1)/Rstat_HL1, label = 'EOS HLPS1')
    ax.scatter(delta_HL2**2, G*mstat_HL2/(rstat_HL2*c**2), (R_HL2-Rstat_HL2)/Rstat_HL2, label = 'EOS HLPS2')
    ax.scatter(delta_HL3**2, G*mstat_HL3/(rstat_HL3*c**2), (R_HL3-Rstat_HL3)/Rstat_HL3, label = 'EOS HLPS3')
    ax.scatter(delta_L**2, G*mstat_L/(rstat_L*c**2), (R_L-Rstat_L)/Rstat_L, label = 'EOS L')
    ax.set_xlabel(r'$\Omega^2 (R_*^3 / GM_*)$', fontsize='15')
    ax.set_ylabel(r'$\zeta$', fontsize='15')
    ax.set_zlabel(r'$(R-R_*)/R_*$', fontsize='15')
    plt.gca().invert_yaxis()
    plt.legend()
    plt.show()
    
elif method == 19:
    ax = fig.gca(projection='3d')
    ax.scatter(delta_AB1**2, G*mstat_AB1/(rstat_AB1*c**2),(M_AB1-Mstat_AB1)/Mstat_AB1, label = 'EOS ABPR1')
    ax.scatter(delta_AB2**2, G*mstat_AB2/(rstat_AB2*c**2),(M_AB2-Mstat_AB2)/Mstat_AB2, label = 'EOS ABPR2')
    ax.scatter(delta_AB3**2, G*mstat_AB3/(rstat_AB3*c**2),(M_AB3-Mstat_AB3)/Mstat_AB3, label = 'EOS ABPR3')
    ax.scatter(delta_AP**2, G*mstat_AP/(rstat_AP*c**2),(M_AP-Mstat_AP)/Mstat_AP, label = 'EOS APR')
    ax.scatter(delta_B1**2, G*mstat_B1/(rstat_B1*c**2),(M_B1-Mstat_B1)/Mstat_B1, label = 'EOS BBB1')
    ax.scatter(delta_B2**2, G*mstat_B2/(rstat_B2*c**2), (M_B2-Mstat_B2)/Mstat_B2, label = 'EOS BBB2')
    ax.scatter(delta_HL1**2, G*mstat_HL1/(rstat_HL1*c**2),(M_HL1-Mstat_HL1)/Mstat_HL1, label = 'EOS HLPS1')
    ax.scatter(delta_HL2**2, G*mstat_HL2/(rstat_HL2*c**2),(M_HL2-Mstat_HL2)/Mstat_HL2, label = 'EOS HLPS2')
    ax.scatter(delta_HL3**2, G*mstat_HL3/(rstat_HL3*c**2),(M_HL3-Mstat_HL3)/Mstat_HL3, label = 'EOS HLPS3')
    ax.scatter(delta_L**2, G*mstat_L/(rstat_L*c**2),(M_L-Mstat_L)/Mstat_L, label = 'EOS L')
    ax.set_xlabel(r'$\Omega^2 (R_*^3 / GM_*)$', fontsize='15')
    ax.set_ylabel(r'$\zeta$', fontsize='15')
    ax.set_zlabel(r'$(M-M_*)/M_*$', fontsize='15')
    ax.view_init(azim=145)
    plt.legend()
    plt.show()
    
elif method == 20:
    ax = fig.gca(projection='3d')
    ax.scatter(delta_AB1**2, G*mstat_AB1/(rstat_AB1*c**2), (M0_AB1-M_AB1)/M0_AB1, label = 'EOS ABPR1')
    ax.scatter(delta_AB2**2, G*mstat_AB2/(rstat_AB2*c**2), (M0_AB2-M_AB2)/M0_AB2, label = 'EOS ABPR2')
    ax.scatter(delta_AB3**2, G*mstat_AB3/(rstat_AB3*c**2), (M0_AB3-M_AB3)/M0_AB3, label = 'EOS ABPR3')
    ax.scatter(delta_AP**2, G*mstat_AP/(rstat_AP*c**2), (M0_AP-M_AP)/M0_AP, label = 'EOS APR')
    ax.scatter(delta_B1**2, G*mstat_B1/(rstat_B1*c**2), (M0_B1-M_B1)/M0_B1, label = 'EOS BBB1')
    ax.scatter(delta_B2**2, G*mstat_B2/(rstat_B2*c**2), (M0_B2-M_B2)/M0_B2, label = 'EOS BBB2')
    ax.scatter(delta_HL1**2, G*mstat_HL1/(rstat_HL1*c**2), (M0_HL1-M_HL1)/M0_HL1, label = 'EOS HLPS1')
    ax.scatter(delta_HL2**2, G*mstat_HL2/(rstat_HL2*c**2), (M0_HL2-M_HL2)/M0_HL2, label = 'EOS HLPS2')
    ax.scatter(delta_HL3**2, G*mstat_HL3/(rstat_HL3*c**2), (M0_HL3-M_HL3)/M0_HL3, label = 'EOS HLPS3')
    ax.scatter(delta_L**2, G*mstat_L/(rstat_L*c**2), (M0_L-M_L)/M0_L, label = 'EOS L')
    ax.set_xlabel(r'$\Omega^2 (R_*^3 / GM_*)$', fontsize='15')
    ax.set_ylabel(r'$\zeta$', fontsize='15')
    ax.set_zlabel(r'$(M_0-M)/M_0$', fontsize='15')
    plt.gca().invert_yaxis()
    plt.legend()
    plt.show()
    
#############################
 
elif method == 21:
    ax = fig.gca(projection='3d')
    ax.scatter(delta_AB1**2, (m_AB1*rstat_AB1)/(mstat_AB1*r_AB1), M_AB1/Mmax_AB1, label = 'EOS ABPR1')
    ax.scatter(delta_AB2**2, (m_AB2*rstat_AB2)/(mstat_AB2*r_AB2), M_AB2/Mmax_AB2, label = 'EOS ABPR2')
    ax.scatter(delta_AB3**2, (m_AB3*rstat_AB3)/(mstat_AB3*r_AB3), M_AB3/Mmax_AB3, label = 'EOS ABPR3')
    ax.scatter(delta_AP**2, (m_AP*rstat_AP)/(mstat_AP*r_AP), M_AP/Mmax_AP, label = 'EOS APR')
    ax.scatter(delta_B1**2, (m_B1*rstat_B1)/(mstat_B1*r_B1), M_B1/Mmax_B1, label = 'EOS BBB1')
    ax.scatter(delta_B2**2, (m_B2*rstat_B2)/(mstat_B2*r_B2), M_B2/Mmax_B2, label = 'EOS BBB2')
    ax.scatter(delta_HL1**2, (m_HL1*rstat_HL1)/(mstat_HL1*r_HL1), M_HL1/Mmax_HL1, label = 'EOS HLPS1')
    ax.scatter(delta_HL2**2, (m_HL2*rstat_HL2)/(mstat_HL2*r_HL2), M_HL2/Mmax_HL2, label = 'EOS HLPS2')
    ax.scatter(delta_HL3**2, (m_HL3*rstat_HL3)/(mstat_HL3*r_HL3), M_HL3/Mmax_HL3, label = 'EOS HLPS3')
    ax.scatter(delta_L**2, (m_L*rstat_L)/(mstat_L*r_L), M_L/Mmax_L, label = 'EOS L')
    ax.set_xlabel(r'$\Omega^2 (R_*^3 / GM_*)$', fontsize='15')
    ax.set_ylabel(r'$MR*/RM*$', fontsize='15')
    ax.set_zlabel(r'$M/M_{max}$', fontsize='15')
    plt.legend()
    plt.show()
    
elif method == 22:
    ax = fig.gca(projection='3d')
    ax.scatter(delta_AB1**2, (m_AB1*rstat_AB1)/(mstat_AB1*r_AB1), (R_AB1-Rstat_AB1)/Rstat_AB1, label = 'EOS ABPR1')
    ax.scatter(delta_AB2**2, (m_AB2*rstat_AB2)/(mstat_AB2*r_AB2), (R_AB2-Rstat_AB2)/Rstat_AB2, label = 'EOS ABPR2')
    ax.scatter(delta_AB3**2, (m_AB3*rstat_AB3)/(mstat_AB3*r_AB3), (R_AB3-Rstat_AB3)/Rstat_AB3, label = 'EOS ABPR3')
    ax.scatter(delta_AP**2, (m_AP*rstat_AP)/(mstat_AP*r_AP), (R_AP-Rstat_AP)/Rstat_AP, label = 'EOS APR')
    ax.scatter(delta_B1**2, (m_B1*rstat_B1)/(mstat_B1*r_B1), (R_B1-Rstat_B1)/Rstat_B1, label = 'EOS BBB1')
    ax.scatter(delta_B2**2, (m_B2*rstat_B2)/(mstat_B2*r_B2), (R_B2-Rstat_B2)/Rstat_B2, label = 'EOS BBB2')
    ax.scatter(delta_HL1**2, (m_HL1*rstat_HL1)/(mstat_HL1*r_HL1), (R_HL1-Rstat_HL1)/Rstat_HL1, label = 'EOS HLPS1')
    ax.scatter(delta_HL2**2, (m_HL2*rstat_HL2)/(mstat_HL2*r_HL2), (R_HL2-Rstat_HL2)/Rstat_HL2, label = 'EOS HLPS2')
    ax.scatter(delta_HL3**2, (m_HL3*rstat_HL3)/(mstat_HL3*r_HL3), (R_HL3-Rstat_HL3)/Rstat_HL3, label = 'EOS HLPS3')
    ax.scatter(delta_L**2, (m_L*rstat_L)/(mstat_L*r_L), (R_L-Rstat_L)/Rstat_L, label = 'EOS L')
    ax.set_xlabel(r'$\Omega^2 (R_*^3 / GM_*)$', fontsize='15')
    ax.set_ylabel(r'$MR_*/RM_*$', fontsize='15')
    ax.set_zlabel(r'$(R-R_*)/R_*$', fontsize='15')
    plt.gca().invert_yaxis()
    plt.legend()
    plt.show()
    
elif method == 23:
    ax = fig.gca(projection='3d')
    ax.scatter(delta_AB1**2, (m_AB1*rstat_AB1)/(mstat_AB1*r_AB1),(M_AB1-Mstat_AB1)/Mstat_AB1, label = 'EOS ABPR1')
    ax.scatter(delta_AB2**2, (m_AB2*rstat_AB2)/(mstat_AB2*r_AB2),(M_AB2-Mstat_AB2)/Mstat_AB2, label = 'EOS ABPR2')
    ax.scatter(delta_AB3**2, (m_AB3*rstat_AB3)/(mstat_AB3*r_AB3),(M_AB3-Mstat_AB3)/Mstat_AB3, label = 'EOS ABPR3')
    ax.scatter(delta_AP**2, (m_AP*rstat_AP)/(mstat_AP*r_AP),(M_AP-Mstat_AP)/Mstat_AP, label = 'EOS APR')
    ax.scatter(delta_B1**2, (m_B1*rstat_B1)/(mstat_B1*r_B1),(M_B1-Mstat_B1)/Mstat_B1, label = 'EOS BBB1')
    ax.scatter(delta_B2**2, (m_B2*rstat_B2)/(mstat_B2*r_B2), (M_B2-Mstat_B2)/Mstat_B2, label = 'EOS BBB2')
    ax.scatter(delta_HL1**2, (m_HL1*rstat_HL1)/(mstat_HL1*r_HL1),(M_HL1-Mstat_HL1)/Mstat_HL1, label = 'EOS HLPS1')
    ax.scatter(delta_HL2**2, (m_HL2*rstat_HL2)/(mstat_HL2*r_HL2),(M_HL2-Mstat_HL2)/Mstat_HL2, label = 'EOS HLPS2')
    ax.scatter(delta_HL3**2, (m_HL3*rstat_HL3)/(mstat_HL3*r_HL3),(M_HL3-Mstat_HL3)/Mstat_HL3, label = 'EOS HLPS3')
    ax.scatter(delta_L**2, (m_L*rstat_L)/(mstat_L*r_L),(M_L-Mstat_L)/Mstat_L, label = 'EOS L')
    ax.set_xlabel(r'$\Omega^2 (R_*^3 / GM_*)$', fontsize='15')
    ax.set_ylabel(r'$MR_*/RM_*$', fontsize='15')
    ax.set_zlabel(r'$(M-M_*)/M_*$', fontsize='15')
    ax.view_init(azim=145)
    plt.legend()
    plt.show()
    
elif method == 24:
    ax = fig.gca(projection='3d')
    ax.scatter(delta_AB1**2, (m_AB1*rstat_AB1)/(mstat_AB1*r_AB1), (M0_AB1-M_AB1)/M0_AB1, label = 'EOS ABPR1')
    ax.scatter(delta_AB2**2, (m_AB2*rstat_AB2)/(mstat_AB2*r_AB2), (M0_AB2-M_AB2)/M0_AB2, label = 'EOS ABPR2')
    ax.scatter(delta_AB3**2, (m_AB3*rstat_AB3)/(mstat_AB3*r_AB3), (M0_AB3-M_AB3)/M0_AB3, label = 'EOS ABPR3')
    ax.scatter(delta_AP**2, (m_AP*rstat_AP)/(mstat_AP*r_AP), (M0_AP-M_AP)/M0_AP, label = 'EOS APR')
    ax.scatter(delta_B1**2, (m_B1*rstat_B1)/(mstat_B1*r_B1), (M0_B1-M_B1)/M0_B1, label = 'EOS BBB1')
    ax.scatter(delta_B2**2, (m_B2*rstat_B2)/(mstat_B2*r_B2), (M0_B2-M_B2)/M0_B2, label = 'EOS BBB2')
    ax.scatter(delta_HL1**2, (m_HL1*rstat_HL1)/(mstat_HL1*r_HL1), (M0_HL1-M_HL1)/M0_HL1, label = 'EOS HLPS1')
    ax.scatter(delta_HL2**2, (m_HL2*rstat_HL2)/(mstat_HL2*r_HL2), (M0_HL2-M_HL2)/M0_HL2, label = 'EOS HLPS2')
    ax.scatter(delta_HL3**2, (m_HL3*rstat_HL3)/(mstat_HL3*r_HL3), (M0_HL3-M_HL3)/M0_HL3, label = 'EOS HLPS3')
    ax.scatter(delta_L**2, (m_L*rstat_L)/(mstat_L*r_L), (M0_L-M_L)/M0_L, label = 'EOS L')
    ax.set_xlabel(r'$\Omega^2 (R_*^3 / GM_*)$', fontsize='15')
    ax.set_ylabel(r'$MR_*/RM_*$', fontsize='15')
    ax.set_zlabel(r'$(M_0-M)/M_0$', fontsize='15')
    plt.gca().invert_yaxis()
    plt.legend()
    plt.show()
    
    
#############################
 
elif method == 25:
    ax = fig.gca(projection='3d')
    ax.scatter(delta_AB1**2, (M_AB1/R_AB1)/Qmax_AB1, M_AB1/Mmax_AB1, label = 'EOS ABPR1')
    ax.scatter(delta_AB2**2, (M_AB2/R_AB2)/Qmax_AB2, M_AB2/Mmax_AB2, label = 'EOS ABPR2')
    ax.scatter(delta_AB3**2, (M_AB3/R_AB3)/Qmax_AB3, M_AB3/Mmax_AB3, label = 'EOS ABPR3')
    ax.scatter(delta_AP**2, (M_AP/R_AP)/Qmax_AP, M_AP/Mmax_AP, label = 'EOS APR')
    ax.scatter(delta_B1**2, (M_B1/R_B1)/Qmax_B1, M_B1/Mmax_B1, label = 'EOS BBB1')
    ax.scatter(delta_B2**2, (M_B2/R_B2)/Qmax_B2, M_B2/Mmax_B2, label = 'EOS BBB2')
    ax.scatter(delta_HL1**2, (M_HL1/R_HL1)/Qmax_HL1, M_HL1/Mmax_HL1, label = 'EOS HLPS1')
    ax.scatter(delta_HL2**2, (M_HL2/R_HL2)/Qmax_HL2, M_HL2/Mmax_HL2, label = 'EOS HLPS2')
    ax.scatter(delta_HL3**2, (M_HL3/R_HL3)/Qmax_HL3, M_HL3/Mmax_HL3, label = 'EOS HLPS3')
    ax.scatter(delta_L**2, (M_L/R_L)/Qmax_L, M_L/Mmax_L, label = 'EOS L')
    ax.set_xlabel(r'$\Omega^2 (R_*^3 / GM_*)$', fontsize='15')
    ax.set_ylabel(r'$(M/R)/(M_{max}/R_{max})$', fontsize='15')
    ax.set_zlabel(r'$M/M_{max}$', fontsize='15')
    plt.legend()
    plt.show()
    
elif method == 26:
    ax = fig.gca(projection='3d')
    ax.scatter(delta_AB1**2, (M_AB1/R_AB1)/Qmax_AB1, (R_AB1-Rstat_AB1)/Rstat_AB1, label = 'EOS ABPR1')
    ax.scatter(delta_AB2**2, (M_AB2/R_AB2)/Qmax_AB2, (R_AB2-Rstat_AB2)/Rstat_AB2, label = 'EOS ABPR2')
    ax.scatter(delta_AB3**2, (M_AB3/R_AB3)/Qmax_AB3, (R_AB3-Rstat_AB3)/Rstat_AB3, label = 'EOS ABPR3')
    ax.scatter(delta_AP**2, (M_AP/R_AP)/Qmax_AP, (R_AP-Rstat_AP)/Rstat_AP, label = 'EOS APR')
    ax.scatter(delta_B1**2, (M_B1/R_B1)/Qmax_B1, (R_B1-Rstat_B1)/Rstat_B1, label = 'EOS BBB1')
    ax.scatter(delta_B2**2, (M_B2/R_B2)/Qmax_B2, (R_B2-Rstat_B2)/Rstat_B2, label = 'EOS BBB2')
    ax.scatter(delta_HL1**2, (M_HL1/R_HL1)/Qmax_HL1, (R_HL1-Rstat_HL1)/Rstat_HL1, label = 'EOS HLPS1')
    ax.scatter(delta_HL2**2, (M_HL2/R_HL2)/Qmax_HL2, (R_HL2-Rstat_HL2)/Rstat_HL2, label = 'EOS HLPS2')
    ax.scatter(delta_HL3**2, (M_HL3/R_HL3)/Qmax_HL3, (R_HL3-Rstat_HL3)/Rstat_HL3, label = 'EOS HLPS3')
    ax.scatter(delta_L**2, (M_L/R_L)/Qmax_L, (R_L-Rstat_L)/Rstat_L, label = 'EOS L')
    ax.set_xlabel(r'$\Omega^2 (R_*^3 / GM_*)$', fontsize='15')
    ax.set_ylabel(r'$(M/R)/(M_{max}/R_{max})$', fontsize='15')
    ax.set_zlabel(r'$(R-R_*)/R_*$', fontsize='15')
    plt.gca().invert_xaxis()
    ax.view_init(azim=145)
    plt.legend()
    plt.show()
    
elif method == 27:
    ax = fig.gca(projection='3d')
    ax.scatter(delta_AB1**2, (M_AB1/R_AB1)/Qmax_AB1,(M_AB1-Mstat_AB1)/Mstat_AB1, label = 'EOS ABPR1')
    ax.scatter(delta_AB2**2, (M_AB2/R_AB2)/Qmax_AB2,(M_AB2-Mstat_AB2)/Mstat_AB2, label = 'EOS ABPR2')
    ax.scatter(delta_AB3**2, (M_AB3/R_AB3)/Qmax_AB3,(M_AB3-Mstat_AB3)/Mstat_AB3, label = 'EOS ABPR3')
    ax.scatter(delta_AP**2, (M_AP/R_AP)/Qmax_AP,(M_AP-Mstat_AP)/Mstat_AP, label = 'EOS APR')
    ax.scatter(delta_B1**2, (M_B1/R_B1)/Qmax_B1,(M_B1-Mstat_B1)/Mstat_B1, label = 'EOS BBB1')
    ax.scatter(delta_B2**2, (M_B2/R_B2)/Qmax_B2, (M_B2-Mstat_B2)/Mstat_B2, label = 'EOS BBB2')
    ax.scatter(delta_HL1**2, (M_HL1/R_HL1)/Qmax_HL1,(M_HL1-Mstat_HL1)/Mstat_HL1, label = 'EOS HLPS1')
    ax.scatter(delta_HL2**2, (M_HL2/R_HL2)/Qmax_HL2,(M_HL2-Mstat_HL2)/Mstat_HL2, label = 'EOS HLPS2')
    ax.scatter(delta_HL3**2, (M_HL3/R_HL3)/Qmax_HL3,(M_HL3-Mstat_HL3)/Mstat_HL3, label = 'EOS HLPS3')
    ax.scatter(delta_L**2, (M_L/R_L)/Qmax_L,(M_L-Mstat_L)/Mstat_L, label = 'EOS L')
    ax.set_xlabel(r'$\Omega^2 (R_*^3 / GM_*)$', fontsize='15')
    ax.set_ylabel(r'$(M/R)/(M_{max}/R_{max})$', fontsize='15')
    ax.set_zlabel(r'$(M-M_*)/M_*$', fontsize='15')
    ax.view_init(azim=145)
    plt.legend()
    plt.show()
    
elif method == 28:
    ax = fig.gca(projection='3d')
    ax.scatter(delta_AB1**2, (M_AB1/R_AB1)/Qmax_AB1, (M0_AB1-M_AB1)/M0_AB1, label = 'EOS ABPR1')
    ax.scatter(delta_AB2**2, (M_AB2/R_AB2)/Qmax_AB2, (M0_AB2-M_AB2)/M0_AB2, label = 'EOS ABPR2')
    ax.scatter(delta_AB3**2, (M_AB3/R_AB3)/Qmax_AB3, (M0_AB3-M_AB3)/M0_AB3, label = 'EOS ABPR3')
    ax.scatter(delta_AP**2, (M_AP/R_AP)/Qmax_AP, (M0_AP-M_AP)/M0_AP, label = 'EOS APR')
    ax.scatter(delta_B1**2, (M_B1/R_B1)/Qmax_B1, (M0_B1-M_B1)/M0_B1, label = 'EOS BBB1')
    ax.scatter(delta_B2**2, (M_B2/R_B2)/Qmax_B2, (M0_B2-M_B2)/M0_B2, label = 'EOS BBB2')
    ax.scatter(delta_HL1**2, (M_HL1/R_HL1)/Qmax_HL1, (M0_HL1-M_HL1)/M0_HL1, label = 'EOS HLPS1')
    ax.scatter(delta_HL2**2, (M_HL2/R_HL2)/Qmax_HL2, (M0_HL2-M_HL2)/M0_HL2, label = 'EOS HLPS2')
    ax.scatter(delta_HL3**2, (M_HL3/R_HL3)/Qmax_HL3, (M0_HL3-M_HL3)/M0_HL3, label = 'EOS HLPS3')
    ax.scatter(delta_L**2, (M_L/R_L)/Qmax_L, (M0_L-M_L)/M0_L, label = 'EOS L')
    ax.set_xlabel(r'$\Omega^2 (R_*^3 / GM_*)$', fontsize='15')
    ax.set_ylabel(r'$(M/R)/(M_{max}/R_{max})$', fontsize='15')
    ax.set_zlabel(r'$(M_0-M)/M_0$', fontsize='15')
    plt.gca().invert_yaxis()
    ax.view_init(azim=145)
    plt.legend()
    plt.show()
    
else:
    print('Enter a number from 1 to 28')
