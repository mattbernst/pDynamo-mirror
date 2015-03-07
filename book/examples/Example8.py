"""Example 8."""

from Definitions import *

# . Define the MM, NB and QC models.
mmModel = MMModelOPLS ( "bookSmallExamples" )
nbModel = NBModelFull ( )
qcModel = QCModelMNDO ( )

# . Define the molecule.
molecule = MOLFile_ToSystem ( os.path.join ( molPath, "waterDimer_cs.mol" ) )

# . Define the selection for the first molecule.
firstWater = Selection.FromIterable ( [ 0, 1, 2 ] )

# . Define the energy model.
molecule.DefineMMModel ( mmModel )
molecule.DefineQCModel ( qcModel, qcSelection = firstWater )
molecule.DefineNBModel ( nbModel )
molecule.Summary ( )

# . Calculate an energy.
molecule.Energy ( )
