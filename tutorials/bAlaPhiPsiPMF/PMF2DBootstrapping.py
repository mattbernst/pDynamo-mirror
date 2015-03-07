"""2-D PMFs - histogram the generated data."""

from Definitions import *

# . Get the system.
system = Unpickle ( os.path.join ( outPath, "bAla.pkl" ) )
system.Summary ( )

# . Get the list of trajectory file names.
fileNames = glob.glob ( os.path.join ( outPath, "bAla_phi_*_psi_*.trj" ) )
fileNames.sort ( )

# . Calculate the PMF.
state = WHAM_Bootstrapping ( fileNames                           ,
                             bins                 = [ 100, 100 ] ,
                             logFrequency         =          1   ,
                             maximumIterations    =       1000   ,
                             resamples            =        100   ,
                             rmsGradientTolerance =     1.0e-3   ,
                             temperature          =     300.0    )

# . Write the PMF to a file.
histogram = state["Histogram"]
pmf       = state["PMF"      ]
sE        = state["Standard Error"]
lower     = state["Lower 95% Confidence Interval"]
upper     = state["Upper 95% Confidence Interval"]
histogram.ToTextFileWithData ( os.path.join ( outPath, "phiPsiPMFWithErrors.dat" ), [ pmf, sE, lower, upper ], format = "{:20.3f} {:20.3f} {:20.3f} {:20.3f} {:20.3f} {:20.3f}\n" )
