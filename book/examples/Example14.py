"""Example 14."""

from Definitions import *

# . Define the molecule and its energy model.
molecule = MOLFile_ToSystem ( os.path.join ( molPath, "cyclohexane_chair.mol" ) )
molecule.coordinates3 = XYZFile_ToCoordinates3 ( os.path.join ( xyzPath, "cyclohexane_chair.xyz" ) )
molecule.DefineMMModel ( MMModelOPLS ( "bookSmallExamples" ) )
molecule.DefineNBModel ( NBModelFull ( ) )
molecule.Summary ( )

# . Calculate the normal modes.
NormalModes_SystemGeometry ( molecule, modify = "project" )

# . Create an output trajectory.
trajectory = SystemGeometryTrajectory ( os.path.join ( scratchPath, "cyclohexane_chair_mode7.trj" ), molecule, mode = "w" )

# . Generate a trajectory for one of the modes.
NormalModesTrajectory_SystemGeometry ( molecule,           \
                                       trajectory,         \
                                       mode   =  7,        \
                                       cycles = 10,        \
                                       frames = 21,        \
                                       temperature = 600.0 )
