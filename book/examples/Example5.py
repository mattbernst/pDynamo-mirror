"""Example 5."""

from Definitions import *

# . Define the energy models.
energyModels = [ QCModelMNDO ( "am1"  ),
                 QCModelMNDO ( "mndo" ),
                 QCModelMNDO ( "pm3"  ) ]

# . Get the filename.
fileName = os.path.join ( xyzPath, "water.xyz" )

# . Loop over the energy models.
results = []
for model in energyModels:
    molecule = XYZFile_ToSystem ( fileName )
    molecule.DefineQCModel ( model )
    molecule.Summary ( )
    energy  = molecule.Energy ( )
    charges = molecule.AtomicCharges ( )
    dipole  = molecule.DipoleMoment  ( )
    results.append ( ( model.label, energy, charges, dipole.Norm2 ( ) ) )

# . Output the results.
table = logFile.GetTable ( columns = [ 10, 20, 20, 20, 20, 20 ] )
table.Start  ( )
table.Title  ( "Energy Model Results for Water" )
table.Heading ( "Model"  )
table.Heading ( "Energy" )
table.Heading ( "Charges", columnSpan = 3 )
table.Heading ( "Dipole" )
for ( label, energy, charges, dipole ) in results:
    table.Entry ( label, alignment = "l" )
    table.Entry ( "{:.1f}".format ( energy ) )
    for charge in charges: table.Entry ( "{:.3f}".format ( charge ) )
    table.Entry ( "{:.3f}".format ( dipole ) )
table.Stop ( )
