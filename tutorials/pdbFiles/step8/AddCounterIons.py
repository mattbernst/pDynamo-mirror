"""Add counterions to the system."""

import os, os.path, sys
sys.path.append ( os.path.join ( os.path.dirname ( os.path.realpath ( __file__ ) ), "..", "data" ) )

from Definitions      import dataPath, outPath
from pBabel           import MOLFile_ToSystem, PDBFile_FromSystem
from pCore            import Pickle, Unpickle
from pMolecule        import MMModelOPLS
from pMoleculeScripts import AddCounterIons

# . Parameters.
# . Box sizes.
_XBox = 75.0
_YBox = 60.0
_ZBox = 60.0

# . Number and type of ions to add.
_NNegative   = 27
_NPositive   = 24
_NegativeIon = "chloride"
_PositiveIon = "potassium"

# . Reorient option.
_Reorient = False

# . Define the solvent MM and NB models.
mmModel = MMModelOPLS ( "bookSmallExamples" )

# . Get the system.
system = Unpickle ( os.path.join ( outPath, "step7.pkl" ) )
system.Summary ( )

# . Reorient the system if necessary (see the results of GetSolvationInformation.py).
masses = system.atoms.GetItemAttributes ( "mass" )
if _Reorient: system.coordinates3.ToPrincipalAxes ( weights = masses )

# . Get the positive and negative ions.
if _NNegative > 0:
    anion = MOLFile_ToSystem ( os.path.join ( dataPath, _NegativeIon + ".mol" ) )
    anion.DefineMMModel ( mmModel )
    anion.Summary ( )
if _NPositive > 0:
    cation = MOLFile_ToSystem ( os.path.join ( dataPath, _PositiveIon + ".mol" ) )
    cation.DefineMMModel ( mmModel )
    cation.Summary ( )

# . Add the counterions.
newSystem = AddCounterIons ( system, _NNegative, anion, _NPositive, cation, ( _XBox, _YBox, _ZBox ) )

# . Save the combined system.
newSystem.configuration.Clear ( )
Pickle ( os.path.join ( outPath, "step8_a.pkl" ), newSystem )

# . Print PDB file.
PDBFile_FromSystem ( os.path.join ( outPath, "step8_a.pdb" ), newSystem )
