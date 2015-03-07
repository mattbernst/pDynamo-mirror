"""Example 13."""

from Definitions import *

# . Define the molecule and its energy model.
molecule = MOLFile_ToSystem ( os.path.join ( molPath, "cyclohexane_chair.mol" ) )
molecule.DefineMMModel ( MMModelOPLS ( "bookSmallExamples" ) )
molecule.DefineNBModel ( NBModelFull ( ) )
molecule.Summary ( )

# . Assign the reactant and product coordinates.
reactants = XYZFile_ToCoordinates3 ( os.path.join ( xyzPath, "cyclohexane_chair.xyz"     ) )
products  = XYZFile_ToCoordinates3 ( os.path.join ( xyzPath, "cyclohexane_twistboat.xyz" ) )

# . Get an initial path.
trajectoryPath = os.path.join ( scratchPath, "cyclohexane_sawpath.trj" )
GrowingStringInitialPath ( molecule, 11, reactants, products, trajectoryPath )

# . Optimization.
trajectory = SystemGeometryTrajectory ( trajectoryPath, molecule, mode = "a+" )
ChainOfStatesOptimizePath_SystemGeometry ( molecule, trajectory ,
                                           logFrequency         =   1 ,
                                           maximumIterations    = 100 ,
                                           rmsGradientTolerance = 0.1 )
