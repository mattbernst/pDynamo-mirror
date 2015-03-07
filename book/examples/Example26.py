"""Example 26."""

from Definitions import *

# . Set some parameters.
NAME        = "example26.trj"
NLAMBDAS    =    21
SOLUTE      =     0
TEMPERATURE = 300.0
DLAMBDA     =   1.0 / float ( NLAMBDAS - 1 )
RT          = CONSTANT_MOLAR_GAS * TEMPERATURE / 1000.0

# . Read in the system.
solution = Unpickle ( os.path.join ( pklPath, "water216_cubicBox_mc.pkl" ) )
solution.Summary ( )

# . Define a random number generator.
randomNumberGenerator = RandomNumberGenerator ( )

# . Initialize the dictionary that will hold the free energy values.
dG = {}

# . Perform simulations at different coupling constants.
for i in range ( NLAMBDAS - 1, -1, -1 ):

    # . Reset the random number generator.
    randomNumberGenerator.SetSeed ( 622199 + i )

    # . Get the value of the coupling parameter.
    LAMBDA = float ( i ) * DLAMBDA

    # . Scale the solute's charge parameters.
    MonteCarlo_ScaleIsolateInteractionParameters ( solution, SOLUTE, chargeScale = LAMBDA )

    # . Equilibration.
    MonteCarlo_SystemGeometry ( solution                                      ,
                                blocks                = 10                    ,
                                moves                 = 100000                ,
                                randomNumberGenerator = randomNumberGenerator ,
                                temperature           = TEMPERATURE           )

    # . Data-collection.
    mcData = SystemGeometryTrajectory ( os.path.join ( scratchPath, NAME ), solution, mode = "w" )
    MonteCarlo_SystemGeometry ( solution                                      ,
                                blocks                = 20                    ,
                                moves                 = 100000                ,
                                randomNumberGenerator = randomNumberGenerator ,
                                temperature           = TEMPERATURE           ,
                                trajectories          = [ ( mcData, 100 ) ]   )

    # . Define a trajectory object for reading.
    mcData = SystemGeometryTrajectory ( os.path.join ( scratchPath, NAME ), solution, mode = "r" )

    # . Initialize the accumulators.
    gB = gF = 0.0

    # . Loop over the frames in the trajectory.
    while mcData.RestoreOwnerData ( ):

        # . Get the interaction energy at i.
        MonteCarlo_ScaleIsolateInteractionParameters ( solution, SOLUTE, chargeScale = LAMBDA )
        eI = MonteCarlo_IsolateInteractionEnergy ( solution, SOLUTE )

        # . Calculate the energy at i-1.
        if i > 0:
            MonteCarlo_ScaleIsolateInteractionParameters ( solution, SOLUTE, chargeScale = LAMBDA - DLAMBDA )
            eJ = MonteCarlo_IsolateInteractionEnergy ( solution, SOLUTE )
            gB += math.exp ( - ( eJ - eI ) / RT )

        # . Calculate the energy at i+1.
        if i < ( NLAMBDAS - 1 ):
            MonteCarlo_ScaleIsolateInteractionParameters ( solution, SOLUTE, chargeScale = LAMBDA + DLAMBDA )
            eJ = MonteCarlo_IsolateInteractionEnergy ( solution, SOLUTE )
            gF += math.exp ( - ( eJ - eI ) / RT )

    # . Scale and save the values.
    if len ( mcData ) == 0:
        gB = 1.0
        gF = 1.0
    else:
        gB /= float ( len ( mcData ) )
        gF /= float ( len ( mcData ) )
    if i > 0:                dG[(i,i-1)] = - RT * math.log ( gB )
    if i < ( NLAMBDAS - 1 ): dG[(i,i+1)] = - RT * math.log ( gF )

# . Output the results.
table = logFile.GetTable ( columns = [ 12, 12, 16, 16, 16 ] )
table.Start   ( )
table.Title   ( "Calculated Free Energies" )
table.Heading ( "Lambda I"     )
table.Heading ( "Lambda J"     )
table.Heading ( "dG (I->J)"    )
table.Heading ( "dG (I<-J)"    )
table.Heading ( "dG (average)" )
dGijTot = dGjiTot = 0.0
for j in range ( NLAMBDAS - 2, -1, -1 ):
    i = j + 1
    dGij = dG[(i,j)]
    dGji = dG[(j,i)]
    dGijTot += dGij
    dGjiTot += dGji
    table.Entry ( "{:12.4f}".format ( float ( i ) * DLAMBDA ) )
    table.Entry ( "{:12.4f}".format ( float ( j ) * DLAMBDA ) )
    table.Entry ( "{:16.4e}".format ( dGij                  ) )
    table.Entry ( "{:16.4e}".format ( dGji                  ) )
    table.Entry ( "{:16.4e}".format ( 0.5 * ( dGij - dGji ) ) )
table.Entry ( "Total:}", alignment = "l", columnSpan = 2 )
table.Entry ( "{:16.3f}".format ( dGijTot ) )
table.Entry ( "{:16.3f}".format ( dGjiTot ) )
table.Entry ( "{:16.3f}".format ( 0.5 * ( dGijTot - dGjiTot ) ) )
table.Stop ( )
