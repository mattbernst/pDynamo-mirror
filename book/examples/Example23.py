"""Example 23."""

from Definitions import *

# . Define some parameters.
DINCREMENT    =   1.0
DMINIMUM      =   1.5
DNAME         = "dOH"
FORCECONSTANT =  20.0
NWINDOWS      =     5

# . Define the atom indices.
OXYGEN   =  5
HYDROGEN = 17

# . Define the MM and NB models.
mmModel = MMModelOPLS ( "bookSmallExamples" )
nbModel = NBModelFull ( )

# . Generate the molecule.
molecule = MOLFile_ToSystem ( os.path.join ( molPath, "bala_c7eq.mol" ) )
molecule.DefineMMModel ( mmModel )
molecule.DefineNBModel ( nbModel )
molecule.Summary ( )

# . Read in the starting coordinates.
molecule.coordinates3 = XYZFile_ToCoordinates3 ( os.path.join ( xyzPath, "bala_1pt5.xyz" ) )

# . Define a constraint container and assign it to the system.
constraints = SoftConstraintContainer ( )
molecule.DefineSoftConstraints ( constraints )

# . Save the molecule definition.
Pickle ( os.path.join ( scratchPath, "bala_example23.pkl" ), molecule )

# . Define a random number generator.
randomNumberGenerator  = RandomNumberGenerator ( )
normalDeviateGenerator = NormalDeviateGenerator.WithRandomNumberGenerator ( randomNumberGenerator )

# . Loop over the values of the distance.
for i in range ( NWINDOWS ):

    # . Reset the random number generator.
    randomNumberGenerator.SetSeed ( 291731 + i )

    # . Calculate the new constraint distance.
    distance = DINCREMENT * float ( i ) + DMINIMUM

    # . Define a new constraint.
    scModel    = SoftConstraintEnergyModelHarmonic ( distance, FORCECONSTANT )
    constraint = SoftConstraintDistance ( OXYGEN, HYDROGEN, scModel )
    constraints[DNAME] = constraint

    # . Equilibration.
    LeapFrogDynamics_SystemGeometry ( molecule                        ,
                                      logFrequency           =   1000 ,
                                      normalDeviateGenerator = normalDeviateGenerator ,
                                      steps                  =  50000 ,
                                      temperature            =  300.0 ,
                                      temperatureControl     =   True ,
                                      temperatureCoupling    =    0.1 ,
                                      timeStep               =  0.001 )

    # . Data-collection.
    trajectory = SystemSoftConstraintTrajectory ( os.path.join ( scratchPath, "bala_window{:d}.trj".format ( i ) ), molecule, mode = "w" )
    LeapFrogDynamics_SystemGeometry ( molecule                     ,
                                      logFrequency        =   1000 ,
                                      steps               = 100000 ,
                                      temperature         =  300.0 ,
                                      temperatureControl  =   True ,
                                      temperatureCoupling =    0.1 ,
                                      timeStep            =  0.001 ,
                                      trajectories        = [ ( trajectory, 1 ) ] )
