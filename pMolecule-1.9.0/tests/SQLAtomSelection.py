"""Test SQL atom selection."""

import glob, os

from pBabel           import PDBFile_ToSystem
from pCore            import TestCase
from pMolecule        import AtomSelection, MMModelOPLS, NBModelABFS, SQLAtomSelector
from pMoleculeScripts import BuildHydrogenCoordinates3FromConnectivity

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class SQLAtomSelectionTest ( TestCase ):
    """Test SQL atom selection."""

    def runTest ( self ):
        """The test."""

        # . Paths.
        dataPath = os.path.join ( os.getenv ( "PDYNAMO_PMOLECULE" ), "data", "pdb" )
        log = self.GetLog ( )

        # . Models.
        mmModel = MMModelOPLS ( "protein" )
        nbModel = NBModelABFS ( )

        # . Get all files.
        pdbFiles = glob.glob ( os.path.join ( dataPath, "*.pdb" ) )
        pdbFiles.sort ( )

        # . Read all PDB files.
        numberFailed = 0
        for pdbFile in pdbFiles:

            if log is not None: log.Text ( "\nProcessing " + pdbFile + ":\n" )
            system = PDBFile_ToSystem ( pdbFile, log = log, useComponentLibrary = True )
            BuildHydrogenCoordinates3FromConnectivity ( system )
            try:

                # . Setup.
                if system.coordinates3.numberUndefined > 0: raise
                system.DefineMMModel ( mmModel, log = log )
                system.DefineNBModel ( nbModel )
                system.Summary ( log = log )
                system.Energy  ( log = log, doGradients = True )

                # . Selection.
                selector              = SQLAtomSelector ( system )

                # . Standard selections.
                aromatics             = selector.aromatics
                backbone              = selector.backbone
                boundaryAtoms         = selector.boundaryAtoms
                counterions           = selector.counterions
                heavyAtoms            = selector.heavyAtoms
                hydrogens             = selector.hydrogens
                mmAtoms               = selector.mmAtoms
                polymerAtoms1         = selector.linearPolymerAtoms
                protein               = selector.protein
                qcAtoms               = selector.qcAtoms
                ringAtoms             = selector.ringAtoms
                water                 = selector.water

                # . Where selections.
                nearOrigin1           = selector.Where ( "X*X + Y*Y + Z*Z < 25.0" )
                positive              = selector.Where ( "Charge > 0.0" )
                threonines1           = selector.Where ( "Path LIKE '%:THR.%:%'" )
                threonines2           = selector.Where ( "ResNam='THR'" )

                # . Atom selection methods.
                nearOrigin2           = nearOrigin1.Within ( 5.0 ).ByComponent ( )
                polymerAtoms2         = threonines1.ByLinearPolymer ( )
                neighbors             = threonines1.ByBondedNeighbor ( iterations = 3 )

                # . Atom selection operators.
                complementNearOrigin2 = ~ nearOrigin2
                null                  = ( nearOrigin2 & complementNearOrigin2 )
                total1                = ( nearOrigin2 | complementNearOrigin2 )
                total2                = ( nearOrigin2 ^ complementNearOrigin2 )

                # . Basic checks.
                n    = len ( system.atoms )
                isOK = ( len ( boundaryAtoms ) == 0 ) and \
                       ( len ( qcAtoms       ) == 0 ) and \
                       ( len ( mmAtoms       ) == n ) and \
                       ( len ( heavyAtoms    )  + len ( hydrogens     ) == n ) and \
                       ( len ( polymerAtoms1 ) == len ( polymerAtoms2 )      ) and \
                       ( len ( counterions   )  + len ( protein       ) + len ( water ) == n ) and \
                       ( len ( threonines1   ) == len ( threonines2   )      ) and \
                       ( len ( null          ) == 0 ) and \
                       ( len ( total1        ) == len ( total2        )      ) and \
                       ( len ( total1        ) == n )
                if not isOK: raise
 
            except Exception as e:
                numberFailed += 1
                if log is not None: log.Text ( "\nError occurred> " +  e.args[0] + "\n" )

        # . Summary of results.
        if log is not None:
            summary = log.GetSummary ( )
            summary.Start ( "SQL Atom Selection Tests" )
            summary.Entry ( "Successes", "{:d}".format ( len ( pdbFiles ) - numberFailed ) )
            summary.Entry ( "Failures" , "{:d}".format (                    numberFailed ) )
            summary.Stop  ( )

        # . Success/failure.
        self.assertTrue ( ( numberFailed == 0 ) )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":

    # . Run the test.
    test = SQLAtomSelectionTest ( )
    test.run ( )
