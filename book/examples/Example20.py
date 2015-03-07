"""Example 20."""

from Definitions import *

# . Define the box side (in Angstroms).
BOXSIDE = 18.641

# . Define the MM and NB models.
mmModel = MMModelOPLS ( "bookSmallExamples" )
nbModel = NBModelABFS ( )

# . Generate the solvent.
solvent = MOLFile_ToSystem ( os.path.join ( molPath, "water216_cubicBox.mol" ) )
solvent.DefineSymmetry ( crystalClass = CrystalClassCubic ( ), a = BOXSIDE )
solvent.DefineMMModel  ( mmModel )
solvent.DefineNBModel  ( nbModel )
solvent.Summary ( )

# . Save the system for later use.
Pickle ( os.path.join ( scratchPath, "water216_cubicBox.pkl" ), solvent )

# . Define a normal deviate generator in a given state.
normalDeviateGenerator = NormalDeviateGenerator.WithRandomNumberGenerator ( RandomNumberGenerator.WithSeed ( 491831 ) )

# . Equilibration.
VelocityVerletDynamics_SystemGeometry ( solvent                                ,
                                        logFrequency              =        500 ,
                                        normalDeviateGenerator    = normalDeviateGenerator ,
                                        steps                     =       5000 ,
                                        timeStep                  =      0.001 ,
                                        temperatureScaleFrequency =        100 ,
                                        temperatureScaleOption    = "constant" ,
                                        temperatureStart          =      300.0 )

# . Data-collection.
trajectory = SystemGeometryTrajectory ( os.path.join ( scratchPath, "water216_cubicBox.trj" ), solvent, mode = "w" )
VelocityVerletDynamics_SystemGeometry ( solvent              ,
                                        logFrequency =   500 ,
                                        steps        = 10000 ,
                                        timeStep     = 0.001 ,
                                        trajectories = [ ( trajectory, 50 ) ] )
