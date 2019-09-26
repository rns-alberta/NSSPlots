# NSSPlots

This code makes use of the output files from the NSSS code that computes a series of constant 
baryonic mass tracks of neutron stars with changing spin velocity.

There are four Python codes: NSS_Plots, NSS_oneEOS, NSS_BestFitSurface and NSS_Residuals. The first one generates 
various plots combining the data obtained from ten different equations of state (EOS); the second one generates the 
same plots as the previous one, but using just one EOS, instead of the ten; the third one creates the best fit surface, 
and its equating for the sequences using one EOS, which are produced by NSSS; lastly, the fourth code computes 
the residuals between each set of sequences and the best fit surface of all the ten EOS sequences. 

To run it from the Terminal the following command line is used 
python NSS_Plots.py
type an integer number and "enter" to obtain the chosen plot from a list on the terminal.
