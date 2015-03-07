"""Example 25."""

from Definitions import *

# . Read in the system.
solution = Unpickle ( os.path.join ( pklPath, "ch4_water215_cubicBox_mc.pkl" ) )
solution.Summary ( )

# . Define a random number generator.
randomNumberGenerator = RandomNumberGenerator.WithSeed ( 899311 )

# . Do a Monte Carlo simulation.
trajectory = SystemGeometryTrajectory ( os.path.join ( scratchPath, "ch4_water215_cubicBox_mc.trj" ), solution, mode = "w" )
MonteCarlo_SystemGeometry ( solution                                        ,
                            blocks                =                     20  ,
                            moves                 =                 100000  ,
                            randomNumberGenerator = randomNumberGenerator   ,
                            trajectories          = [ ( trajectory, 100 ) ] )
