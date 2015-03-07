"""Read a PDB file and convert it to a PDB model."""

import os, os.path, sys
sys.path.append ( os.path.join ( os.path.dirname ( os.path.realpath ( __file__ ) ), "..", "data" ) )

from Definitions import dataPath, outPath
from pBabel      import PDBFile_ToPDBModel, PDBModel_ToModelFile
from pCore       import logFile

# . File names.
fileName = os.path.join ( dataPath, "1CDK.pdb" )

# . Read the PDB file and create a PDB model.
model = PDBFile_ToPDBModel ( fileName )
logFile.LineBreak ( )
model.Summary ( )

# . Write out the PDB model.
( head, tail ) = os.path.split ( fileName )
outName        = os.path.join ( outPath, tail[0:-4] + ".model" )
PDBModel_ToModelFile ( outName, model )
