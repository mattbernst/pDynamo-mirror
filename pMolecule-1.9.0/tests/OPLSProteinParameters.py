"""Test OPLS protein parameters."""

import glob, os

from pBabel  import MOLFile_ToSystem, XYZFile_ToCoordinates3
from pCore   import TestCase
from pMolecule import MMModelOPLS, NBModelFull

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class OPLSProteinParameterTest ( TestCase ):
    """Test OPLS protein parameters."""

    def runTest ( self ):
        """The test."""

        # . Paths.
        dataPath = os.path.join ( os.getenv ( "PDYNAMO_ROOT" ), "molecularStructures", "aminoAcids", "mol" )
        log = self.GetLog ( )

        # . Models.
        mmModel = MMModelOPLS ( "protein" )
        nbModel = NBModelFull ( )

        # . Get all files.
        molFiles = glob.glob ( os.path.join ( dataPath, "*.mol" ) )
        molFiles.sort ( )

        # . Read all mol files.
        numberFailed = 0
        for molFile in molFiles:

            if log is not None: log.Text ( "\nProcessing " + molFile + ":\n" )
            molecule = MOLFile_ToSystem ( molFile )
            try:
                molecule.DefineMMModel ( mmModel, log = log )
                molecule.DefineNBModel ( nbModel )
                molecule.Summary ( log = log )
                molecule.Energy  ( log = log, doGradients = True )
            except Exception as e:
                numberFailed += 1
                if log is not None: log.Text ( "\nError occurred> " +  e.args[0] + "\n" )

        # . Summary of results.
        if log is not None:
            summary = log.GetSummary ( )
            summary.Start ( "OPLS Protein Parameter Tests" )
            summary.Entry ( "Successes", "{:d}".format ( len ( molFiles ) - numberFailed ) )
            summary.Entry ( "Failures" , "{:d}".format (                    numberFailed ) )
            summary.Stop  ( )

        # . Success/failure.
        self.assertTrue ( ( numberFailed == 0 ) )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":

    # . Run the test.
    test = OPLSProteinParameterTest ( )
    test.run ( )
