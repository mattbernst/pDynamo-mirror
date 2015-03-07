"""Example 27."""

from Definitions import *

# . Define some parameters.
MOLECULENAME =    "water"
NLINEAR      =          6
NMOLECULES   = NLINEAR**3

# . Define the MM and NB models.
mmModel = MMModelOPLS ( "bookSmallExamples" )
nbModel = NBModelMonteCarlo ( )

# . Define the solvent molecule.
molecule = MOLFile_ToSystem ( os.path.join ( molPath, MOLECULENAME + ".mol" ) )
molecule.Summary ( )

# . Build the cubic system.
solvent = BuildCubicSolventBox ( molecule, NMOLECULES )
solvent.label = "Water Box"
solvent.DefineMMModel ( mmModel )
solvent.DefineNBModel ( nbModel )
solvent.Summary ( )

# . Do Monte Carlo simulations to equilibrate the system.
MonteCarlo_SystemGeometry ( solvent           ,
                            blocks   =     10 ,
                            moves    = 100000 ,
                            pressure = 1000.0 )
MonteCarlo_SystemGeometry ( solvent           ,
                            blocks   =     10 ,
                            moves    = 100000 )

# . Calculate and print the final density.
mass    = sum ( solvent.atoms.GetItemAttributes ( "mass" ) )
volume  = solvent.symmetryParameters.volume
density = ( mass / volume ) * ( UNITS_MASS_AMU_TO_KG * 1.0e+30 )
logFile.Paragraph ( "Solvent density = {:.2f} kg m^-3.".format ( density ) )

# . Save the system.
solvent.configuration.Clear ( )
Pickle ( os.path.join ( scratchPath, "{:s}{:d}_cubicBox_mc.pkl".format ( MOLECULENAME, NMOLECULES ) ), solvent )
