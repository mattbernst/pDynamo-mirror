"""Definition of simple QC/MM systems for testing."""

import math, os.path

from pBabel    import MOLFile_ToSystem
from pCore     import logFile, LogFileActive, Selection
from pMolecule import ElectronicState, MMModelOPLS, NBModelFull

#===================================================================================================================================
# . Parameters for system definitions.
#===================================================================================================================================
# . Paths.
_dataPath = os.path.join ( os.getenv ( "PDYNAMO_PMOLECULE" ), "data", "mol" )

# . Molecule definitions.
_keywordLabels = ( "label", "dataPath", "fileName", "qcCharge", "multiplicity", "qcSelection" )
_moleculeData  = ( ( "Alanine Dipeptide"   , _dataPath, "bAla_c7eq"        , 0, 1, Selection.FromIterable ( range ( 10, 14 ) ) ), \
                   ( "Cyclohexane 6"       , _dataPath, "cyclohexane"      , 0, 1, Selection.FromIterable ( range (  6     ) ) ), \
                   ( "Cyclohexane 9"       , _dataPath, "cyclohexane"      , 0, 1, Selection.FromIterable ( range (  9     ) ) ), \
                   ( "Tyrosine Dipeptide 1", _dataPath, "tyrosineDipeptide", 0, 1, Selection.FromIterable ( range (  6, 14 ) + range ( 22, 29 ) ) ), \
                   ( "Tyrosine Dipeptide 2", _dataPath, "tyrosineDipeptide", 1, 2, Selection.FromIterable ( range (  6, 14 ) + range ( 22, 29 ) ) ), \
                   ( "Water Dimer 1"       , _dataPath, "waterDimer_cs"    , 0, 1, Selection.FromIterable ( range (  3     ) ) ), \
                   ( "Water Dimer 2"       , _dataPath, "waterDimer_cs"    , 0, 1, Selection.FromIterable ( range (  3,  6 ) ) )  )

#===================================================================================================================================
# . Class for a QC/MM test system.
#===================================================================================================================================
class QCMMTestSystem ( object ):
    """QC/MM test system."""

    def __init__ ( self, **kwargs ):
        """Constructor."""
        for ( attribute, value ) in kwargs.iteritems ( ):
            setattr ( self, attribute, value )
        self.mmModel  = MMModelOPLS ( "protein" )
        self.nbModel  = NBModelFull ( )

    def GetSystem ( self, doQCMM = True, log = logFile, nbModel = None, qcModel = None ):
        """Get the system with the energy model defined."""
        # . Basic setup.
        molecule       = MOLFile_ToSystem  ( os.path.join ( self.dataPath, self.fileName + ".mol" ) )
        molecule.label = self.label
        molecule.DefineMMModel ( self.mmModel )
        # . Set up the QC model.
        if qcModel is not None:
            molecule.electronicState = ElectronicState ( charge = self.qcCharge, multiplicity = self.multiplicity )
            if doQCMM: molecule.DefineQCModel ( qcModel, qcSelection = self.qcSelection )
            else:      molecule.DefineQCModel ( qcModel )
        # . Set up the NB model.
        if ( qcModel is None ) or doQCMM:
            if nbModel is None: molecule.DefineNBModel  ( self.nbModel )
            else:               molecule.DefineNBModel  (      nbModel )
        # . Summary.
        if LogFileActive ( log ):
            molecule.Summary ( log = log )
            log.Paragraph ( "\nFormula = " + molecule.atoms.FormulaString ( ) + "." )
        # . Finish up.
        return molecule

#===================================================================================================================================
# . Set up the test systems (as a dictionary).
#===================================================================================================================================
qcmmTestSystems = {}
for values in _moleculeData:
    kwargs = { key : value for ( key, value ) in zip ( _keywordLabels, values ) }
    qcmmTestSystems[kwargs["label"]] = QCMMTestSystem ( **kwargs )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
