"""Example 19."""

from Definitions import *

# . Define various parameters.
BUFFER     = 4.0
CINCREMENT = 1.0
CSTART     = 0.0
NENERGIES  =  40

# . Define the energy models.
mmModel = MMModelOPLS ( "protein" )
nbFull  = NBModelFull ( )
nbABFS  = NBModelABFS ( )

# . Set up the system.
molecule = PDBFile_ToSystem ( os.path.join ( pdbPath, "crambin.pdb" ), useComponentLibrary = True )
molecule.DefineMMModel ( mmModel )
molecule.DefineNBModel ( nbFull  )
molecule.Summary ( )

# . Stop if not all coordinates are defined.
BuildHydrogenCoordinates3FromConnectivity ( molecule )
if molecule.coordinates3.numberUndefined > 0:
    logFile.Paragraph ( "The system has undefined coordinates." )
# . Continue.
else:
    # . Get the energy with a full model.
    eF = molecule.Energy ( log = None )

    # . Reset the NB model for the molecule.
    molecule.DefineNBModel ( nbABFS )

    # . Initialize the cutoff.
    cut = CSTART

    # . Output the energy difference for each cutoff.
    table = logFile.GetTable ( columns = [ 20, 20 ] )
    table.Start   ( )
    table.Title   ( "Cutoff/Full Energy Difference" )
    table.Heading ( "Inner Cutoff" )
    table.Heading ( "Difference"   )
    for i in range ( NENERGIES ):
        cut += CINCREMENT
        nbABFS.ClearPairwiseInteractions ( )
        nbABFS.SetOptions ( innerCutoff = cut, outerCutoff = cut + BUFFER, listCutoff = cut + BUFFER )
        eT   = molecule.Energy ( log = None )
        table.Entry ( "{:.1f}".format ( cut     ) )
        table.Entry ( "{:.4f}".format ( eT - eF ) )
    table.Stop ( )
