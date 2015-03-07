"""Example 27.

Updated version using different box builder and Monte Carlo refinement.
"""

from Definitions import *

# . Define some parameters.
_DENSITY      = 1000.0 # . Density of water (kg m^-3).
_MOLECULENAME = "water"
_NMOLECULES   = 216

# . Define the solvent MM and NB models.
mmModel = MMModelOPLS ( "bookSmallExamples" )
nbModel = NBModelMonteCarlo ( )

# . Define the solvent molecule.
molecule = MOLFile_ToSystem ( os.path.join ( molPath, _MOLECULENAME + ".mol" ) )
molecule.Summary ( )

# . Get the dimensions of a cubic box with the required number of molecules.
a = SolventCubicBoxDimensions ( molecule, _NMOLECULES, _DENSITY )

# . Create a symmetryParameters instance with the correct dimensions.
symmetryParameters = SymmetryParameters ( )
symmetryParameters.SetCrystalParameters ( a, a, a, 90.0, 90.0, 90.0 )

# . Create the basic solvent box.
solvent = BuildSolventBox ( CrystalClassCubic ( ), symmetryParameters, molecule, _DENSITY )
solvent.label = "Water Box"
solvent.DefineMMModel ( mmModel )
solvent.DefineNBModel ( nbModel )
solvent.Summary ( )

# . Do Monte Carlo simulations to equilibrate the system.
MonteCarlo_SystemGeometry ( solvent,                  \
                            blocks          =     50, \
                            log             =   None, \
                            moves           =   1000, \
                            volumeFrequency =       0 )
MonteCarlo_SystemGeometry ( solvent,                  \
                            blocks          =     20, \
                            moves           = 100000  )

# . Calculate and print the final density.
logFile.Paragraph ( "Solvent density = {:.2f} kg m^-3.".format ( SystemDensity ( solvent ) ) )

# . Save the system.
solvent.configuration.Clear ( )
Pickle ( os.path.join ( scratchPath, "{:s}{:d}_cubicBox_mc.pkl".format ( _MOLECULENAME, _NMOLECULES ) ), solvent )
