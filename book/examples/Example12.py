"""Example 12."""

from Definitions import *

# . Define the molecule and its energy model.
molecule = MOLFile_ToSystem ( os.path.join ( molPath, "cyclohexane_chair.mol" ) )
molecule.coordinates3 = XYZFile_ToCoordinates3 ( os.path.join ( scratchPath, "cyclohexane_saddle.xyz" ) )
molecule.DefineMMModel ( MMModelOPLS ( "bookSmallExamples" ) )
molecule.DefineNBModel ( NBModelFull ( ) )
molecule.Summary ( )

# . Calculate an energy.
molecule.Energy ( )

# . Create an output trajectory.
trajectory = SystemGeometryTrajectory ( os.path.join ( scratchPath, "cyclohexane_sdpath.trj" ), molecule, mode = "w" )

# . Optimization.
SteepestDescentPath_SystemGeometry ( molecule,                       \
                                     functionStep      = 2.0,        \
                                     logFrequency      =  10,        \
                                     maximumIterations = 400,        \
                                     pathStep          = 0.025,      \
                                     saveFrequency     =  10,        \
                                     trajectory        = trajectory, \
                                     useMassWeighting  = True        )
