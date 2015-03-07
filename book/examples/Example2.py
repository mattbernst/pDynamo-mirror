"""Example 2."""

from Definitions import *

# . Define the molecule.
moleculeName = "bala_c7eq"
smiles       = "CC(=O)NC(C)C(=O)NC"

# . Initialize a list to contain the molecules.
molecules = []

# . Read all molecules.
molecules.append ( MOLFile_ToSystem ( os.path.join ( molPath, moleculeName + ".mol" ) ) )
molecules.append ( PDBFile_ToSystem ( os.path.join ( pdbPath, moleculeName + ".pdb" ) ) )
molecules.append ( XYZFile_ToSystem ( os.path.join ( xyzPath, moleculeName + ".xyz" ) ) )

# . Generate a molecule from a SMILES string.
molecules.append ( SMILES_ToSystem ( smiles ) )

# . Print summaries of the molecules.
for molecule in molecules:
    molecule.Summary ( )
