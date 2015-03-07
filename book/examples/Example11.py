"""Example 11."""

from Definitions import *

# . Define the molecule and its energy model.
molecule = MOLFile_ToSystem ( os.path.join ( molPath, "cyclohexane_chair.mol" ) )
molecule.DefineMMModel ( MMModelOPLS ( "bookSmallExamples" ) )
molecule.DefineNBModel ( NBModelFull ( ) )
molecule.Summary ( )

# . Determine the starting energy.
eStart = molecule.Energy ( )

# . Optimization.
BakerSaddleOptimize_SystemGeometry ( molecule,
                                     logFrequency         =      1 ,
                                     maximumIterations    =    250 ,
                                     rmsGradientTolerance = 1.0e-3 )

# . Determine the final energy.
eStop = molecule.Energy ( )

# . Print the energy change.
logFile.Paragraph ( "Energy change after search = {:20.4f}\n".format ( eStop - eStart ) )

# . Save the coordinates.
molecule.label = "Cyclohexane saddle conformation - OPLS-AA optimized."
XYZFile_FromSystem ( os.path.join ( scratchPath, "cyclohexane_saddle.xyz" ), molecule )
