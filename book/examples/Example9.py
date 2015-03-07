"""Example 9."""

from Definitions import *

# . Define the MM, NB and QC models.
mmModel = MMModelOPLS ( "bookSmallExamples" )
nbModel = NBModelFull ( )
qcModel = QCModelMNDO ( )

# . Define the molecule.
molecule = MOLFile_ToSystem ( os.path.join ( molPath, "bala_c7eq.mol" ) )

# . Define the selection for the first molecule.
methylGroup = Selection.FromIterable ( [ 10, 11, 12, 13 ] )

# . Define the energy model.
molecule.DefineMMModel ( mmModel )
molecule.DefineQCModel ( qcModel, qcSelection = methylGroup )
molecule.DefineNBModel ( nbModel )
molecule.Summary ( )

# . Calculate an energy.
molecule.Energy ( )
