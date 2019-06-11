# NS_data

Output files from the sequences.c code that computes a series of constant 
baryonic mass tracks of neutron stars with changing spin velocity. Each column in the file represents:

1. Central energy density
2. Total mass
3. Baryonic mass
4. Mass when the NS is not rotating 
5. Maximum mass of a non rotating NS for a given EOS (obtained from maxmass.c)
6. Radius 
7. Ratio r_p/r_e
8. Radius when the NS is not rotating
9. Spin frequency
10. Kepler limit
11. Angular momentum
12. Rotational kinetic energy
13. Gravitational binding energy
14. Maximum radius of the non rotating star for a given EOS (obtained from maxmass.c)
15. Ratio of the maximum mass and the maximum radius of the non rotating star

The python (version 2.7) code Plots.py generates various plots combining the data from the ten files. 
To run it from the Terminal the following command line is used "python Plots.py" and by typing a number from 1 to 18 
and "enter" we get the chosen plot.

The code "Plot_oneEOS.py" can be run in the same way, but this only generates the plot using one equation of state, 
instead of the 10 in the code "Plots.py"
