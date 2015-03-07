"""Test to read and write SMILES."""

# . The SMILES module needs work.
# . No attempt is made to make sure everything is OK.

import glob, os, os.path

from pBabel import SMILESConnectivityError, SMILESReaderError, SMILES_FromSystem, SMILES_ToSystem
from pCore  import logFile, LogFileActive, TestCase

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class SMILESReadWriteTest ( TestCase ):
    """A test case for reading and writing SMILES."""

    def runTest ( self ):
        """The test."""

        # . Initialization.
        isOK         = True
        log          = self.GetLog ( )
        numberErrors = 0

        # . Define the SMILES to process.
        smilesPath  = os.path.join ( os.getenv ( "PDYNAMO_PBABEL" ), "data", "smiles" )
        smilesFiles = glob.glob ( os.path.join ( smilesPath, "*.smi*" ) )
        smiles = []
        for smilesFile in smilesFiles:
            for line in open ( smilesFile, "r" ):
                words = line.split ( )
                if ( len ( words ) > 0 ) and ( not words[0].startswith ( "#" ) ): smiles.append ( words[0] )

        # . Loop over the examples.
        results = []
        if log is not None: log.Separator ( )
        for oldSmiles in smiles:
            try:
                system    = SMILES_ToSystem ( oldSmiles, log = log )
                system.Summary ( log = log )
                newSmiles = SMILES_FromSystem ( system, log = log )
                results.append ( ( oldSmiles, newSmiles ) )
            except ( SMILESConnectivityError, SMILESReaderError ) as e:
                if log is not None: log.Text ( "\nSMILES error> " + e.args[0] + "\n" )
            except Exception as e:
                numberErrors += 1
                if log is not None: log.Text ( "\nError occurred> " +  e.args[0] + "\n" )
            if log is not None: log.Separator ( )

        # . Write out the old and new SMILES.
        if log is not None:
            oldLength = 20 # . Minimum widths.
            newLength = 20
            for ( oldSmiles, newSmiles ) in results:
                oldLength = max ( oldLength, len ( oldSmiles ) )
                newLength = max ( newLength, len ( newSmiles ) )

            # . Define a table for the results.
            table = log.GetTable ( columns = [ oldLength + 2, newLength + 2 ] )
            table.Start ( )
            table.Title ( "SMILES Comparison" )
            table.Heading ( "Old" )
            table.Heading ( "New" )
            for ( oldSmiles, newSmiles ) in results:
                table.Entry ( oldSmiles, alignment = "left" )
                table.Entry ( newSmiles, alignment = "left" )
            table.Stop ( )

        # . Success/failure.
        self.assertTrue ( isOK and ( numberErrors == 0 ) )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":

    # . Run the test.
    test = SMILESReadWriteTest ( )
    test.runTest ( )
