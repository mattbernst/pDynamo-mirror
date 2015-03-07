"""Example 15."""

from Definitions import *

# . Methods.
def FreeEnergies ( fileName, temperatures, symmetryNumber = 1 ):
    """Calculate the potential energy for a system and its
       Gibbs free energies at several temperatures."""

    # . Define the molecule and its energy model.
    molecule = MOLFile_ToSystem ( os.path.join ( molPath, fileName + ".mol" ) )
    molecule.coordinates3 = XYZFile_ToCoordinates3 ( os.path.join ( xyzPath, fileName + ".xyz" ) )
    molecule.DefineMMModel ( MMModelOPLS ( "bookSmallExamples" ) )
    molecule.DefineNBModel ( NBModelFull ( ) )
    molecule.Summary ( )

    # . Calculate the energy and normal modes.
    e = molecule.Energy ( )
    NormalModes_SystemGeometry ( molecule, modify = "project" )

    # . Loop over the temperatures.
    g = []
    for T in temperatures:
        tdics = ThermodynamicsRRHO_SystemGeometry ( molecule,                        \
                                                    pressure       = 1.0,            \
                                                    symmetryNumber = symmetryNumber, \
                                                    temperature    = T               )
        g.append ( tdics["Gibbs Free Energy"] )

    # . Return the energies.
    return ( e, g )

# . Create a sequence of temperatures.
temperatures = [ 100.0 * i for i in range ( 1, 11 ) ]

# . Get the energies for the boat and chair structures.
( eB, gBoat  ) = FreeEnergies ( "cyclohexane_twistboat", temperatures, symmetryNumber = 4 )
( eC, gChair ) = FreeEnergies ( "cyclohexane_chair",     temperatures, symmetryNumber = 6 )
deltaE = ( eB - eC )

# . Output the equilibrium constants.
table = logFile.GetTable ( columns = [ 25, 25 ] )
table.Start   ( )
table.Title   ( "Equilibrium Constants (Chair -> Twist Boat)" )
table.Heading ( "Temperature" )
table.Heading ( "Log K"       )
for ( T, gC, gB ) in zip ( temperatures, gChair, gBoat ):
    RT     = ( CONSTANT_MOLAR_GAS * T ) / 1000.0
    log10K = math.log10 ( math.e ) * ( - ( gB - gC + deltaE ) / RT )
    table.Entry ( "{:.4f}".format ( T      ) )
    table.Entry ( "{:.6g}".format ( log10K ) )
table.Stop ( )
