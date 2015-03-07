"""Modification of Example 5 that includes some DFT QC Models."""

import os, os.path

from pBabel    import XYZFile_ToSystem
from pCore     import logFile
from pMolecule import DIISSCFConverger, QCModelDFT, QCModelMNDO

# . The DFT QC Model is QCModelDFT. Its constructor takes the following arguments:
#
#   accuracy         - the DFT accuracy. One of "Very Low", "Low", "Medium" (default), "High", "Very High".
#   converger        - the SCF converger (currently only DIIS).
#   densityBasis     - the basis set for Coulomb fitting. The default is "demon".
#   functional       - the DFT functional. Currently only "lda" (default) and "blyp".
#   orbitalBasis     - the orbital basis set. The default is "321g".
#   qcChargeModel    - the charge model to use in QC/MM interactions. Use the default for the moment which is "Lowdin".
#
# . Good combinations of orbital and density basis sets are (in increasing order of cost):
#
#   321g/demon; 631gs/ahlrichs; svp/weigend.
#
# . To perform QC/MM calculations, proceed in the usual way. Any of the NB models that can be used with
#   QCModelMNDO can also be used with QCModelDFT.
#

# . Define the SCF converger.
converger = DIISSCFConverger ( densityTolerance = 1.0e-8, maximumSCFCycles = 25 )

# . Define the energy models.
energyModels = [ QCModelMNDO ( "am1"  ),
                 QCModelDFT  ( converger = converger, densityBasis = "demon",    functional = "lda",  orbitalBasis = "321g"  ), \
                 QCModelDFT  ( converger = converger, densityBasis = "demon",    functional = "blyp", orbitalBasis = "321g"  ), \
                 QCModelDFT  ( converger = converger, densityBasis = "ahlrichs", functional = "lda",  orbitalBasis = "631gs" ), \
                 QCModelDFT  ( converger = converger, densityBasis = "ahlrichs", functional = "blyp", orbitalBasis = "631gs" ), \
                 QCModelDFT  ( converger = converger, densityBasis = "weigend",  functional = "lda",  orbitalBasis = "svp"   ), \
                 QCModelDFT  ( converger = converger, densityBasis = "weigend",  functional = "blyp", orbitalBasis = "svp"   )  ]

# . Get the fileName.
fileName = os.path.join ( os.getenv ( "PDYNAMO_PMOLECULE" ), "data", "xyz", "water.xyz" )

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
table = logFile.GetTable ( columns = [ 20, 20, 20, 20, 20, 20 ] )
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
