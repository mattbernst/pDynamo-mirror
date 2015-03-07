"""Setup the proteins."""

from Definitions import *

# . PDB files.
_PDBPaths = ( "1UAO", "2E4E" )

# . Loop over structures.
for pdbPath in _PDBPaths:

    # . Define the energy models.
    mmModel = MMModelOPLS ( "protein" )
    nbModel = NBModelABFS ( )

    # . Set up the systems.
    system = PDBFile_ToSystem ( os.path.join ( dataPath, pdbPath + ".pdb" ), modelNumber = 1, useComponentLibrary = True )
    BuildHydrogenCoordinates3FromConnectivity ( system )
    system.DefineMMModel ( mmModel )
    system.DefineNBModel ( nbModel )
    system.Summary ( )

    # . Get an initial energy.
    system.Energy ( doGradients = True )

    # . Save.
    mmModel.ClearModelBuildingData ( )
    system.configuration.Clear ( )
    Pickle ( os.path.join ( outPath, pdbPath + ".pkl" ), system )
