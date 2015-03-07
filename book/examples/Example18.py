"""Example 18."""

from Definitions import *

# . Define various parameters.
CLUSTERSIZE     = 3.0
ENERGYTOLERANCE = 0.01
FORCECONSTANT   = 5.0
NTRIALS         = 100

# . Set up the system.
cluster = MOLFile_ToSystem ( os.path.join ( molPath, "argon13.mol" ) )
cluster.DefineMMModel ( MMModelOPLS ( "lennardJones" ) )
cluster.DefineNBModel ( NBModelFull ( ) )
cluster.Summary ( )

# . Define tether constraints for each atom.
origin            = Vector3.Null ( )
tetherEnergyModel = SoftConstraintEnergyModelHarmonicRange ( 0.0, 0.5 * CLUSTERSIZE, FORCECONSTANT )
tethers           = SoftConstraintContainer ( )
for i in range ( len ( cluster.atoms ) ):
    tethers["{:d}".format ( i )] = SoftConstraintTether ( i, origin, tetherEnergyModel )

# . Define a random number generator.
randomNumberGenerator  = RandomNumberGenerator ( )
normalDeviateGenerator = NormalDeviateGenerator.WithRandomNumberGenerator ( randomNumberGenerator )

# . Initialize lists to keep energies.
pe0 = []
pe1 = []
pe2 = []
pe3 = []

# . Loop over the trials.
for i in range ( NTRIALS ):

    # . Reset the cluster coordinates and cluster constraints.
    cluster.coordinates3 = MOLFile_ToCoordinates3 ( os.path.join ( molPath, "argon13.mol" ) )
    cluster.DefineSoftConstraints ( tethers )

    # . Initialize some variables for the trial.
    randomNumberGenerator.SetSeed ( 957612 + i )
    tStart = 300.0 * ( randomNumberGenerator.NextReal ( ) + 1.0 )

    # . Do a short dynamics to generate a new structure.
    VelocityVerletDynamics_SystemGeometry ( cluster                                ,
                                            log                       =       None ,
                                            normalDeviateGenerator    = normalDeviateGenerator ,
                                            steps                     =      10000 ,
                                            timeStep                  =      0.001 ,
                                            temperatureScaleFrequency =        100 ,
                                            temperatureScaleOption    = "constant" ,
                                            temperatureStart          =     tStart )

    # . Save the starting coordinates and energy.
    temporary3 = Clone ( cluster.coordinates3 )
    cluster.DefineSoftConstraints ( None )
    pe0.append ( cluster.Energy ( log = None ) )

    # . Minimization.
    cluster.DefineSoftConstraints ( tethers )
    ConjugateGradientMinimize_SystemGeometry ( cluster                       ,
                                               log                  =  None  ,
                                               maximumIterations    = 10000  ,
                                               rmsGradientTolerance = 1.0e-4 )
    cluster.DefineSoftConstraints ( None )
    ConjugateGradientMinimize_SystemGeometry ( cluster                       ,
                                               log                  =  None  ,
                                               maximumIterations    = 10000  ,
                                               rmsGradientTolerance = 1.0e-4 )
    pe1.append ( cluster.Energy ( log = None ) )

    # . Simulated annealing from the starting coordinates.
    cluster.coordinates3 = temporary3
    cluster.DefineSoftConstraints ( tethers )
    VelocityVerletDynamics_SystemGeometry ( cluster                                   ,
                                            log                       =          None ,
                                            steps                     =         40000 ,
                                            timeStep                  =         0.001 ,
                                            temperatureScaleFrequency =            10 ,
                                            temperatureScaleOption    = "exponential" ,
                                            temperatureStart          = tStart        ,
                                            temperatureStop           = tStart * math.exp ( - 10.0 ) )
    cluster.DefineSoftConstraints ( None )
    pe2.append ( cluster.Energy ( log = None ) )

    # . Minimization of the annealed structure.
    ConjugateGradientMinimize_SystemGeometry ( cluster                       ,
                                               log                  =  None  ,
                                               maximumIterations    = 10000  ,
                                               rmsGradientTolerance = 1.0e-4 )
    pe3.append ( cluster.Energy ( log = None ) )

# . Prepare the energies for statistics.
stpe1 = Statistics ( pe1 )
stpe2 = Statistics ( pe2 )
stpe3 = Statistics ( pe3 )

# . Output the results.
table = logFile.GetTable ( columns = [ 10, 20, 20, 20, 20 ] )
table.Start   ( )
table.Title   ( "Optimization Results" )
table.Heading ( "Attempt"          )
table.Heading ( "Initial Energy"   )
table.Heading ( "Minimized Energy" )
table.Heading ( "Annealed Energy"  )
table.Heading ( "Final Energy"     )
for i in range ( NTRIALS ):
    table.Entry ( "{:d}"    .format (     i  ) )
    table.Entry ( "{:20.3f}".format ( pe0[i] ) )
    table.Entry ( "{:20.3f}".format ( pe1[i] ) )
    table.Entry ( "{:20.3f}".format ( pe2[i] ) )
    table.Entry ( "{:20.3f}".format ( pe3[i] ) )
table.Entry ( "Minimum Energies:", alignment = "l", columnSpan = 2 )
table.Entry ( "{:20.3f}".format ( stpe1.minimum ) )
table.Entry ( "{:20.3f}".format ( stpe2.minimum ) )
table.Entry ( "{:20.3f}".format ( stpe3.minimum ) )
table.Entry ( "Frequencies:",      alignment = "l", columnSpan = 2 )
table.Entry ( "{:d}".format ( stpe1.Count ( stpe1.minimum, tolerance = ENERGYTOLERANCE ) ) )
table.Entry ( "{:d}".format ( stpe2.Count ( stpe2.minimum, tolerance = ENERGYTOLERANCE ) ) )
table.Entry ( "{:d}".format ( stpe3.Count ( stpe3.minimum, tolerance = ENERGYTOLERANCE ) ) )
table.Stop ( )
