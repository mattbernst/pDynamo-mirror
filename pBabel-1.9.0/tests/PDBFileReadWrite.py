"""Read and write PDB files."""

import glob, os, os.path

from pBabel           import PDBFile_FromSystem, PDBFile_ToSystem
from pCore            import logFile, LogFileActive, TestCase
from pMoleculeScripts import PrintComponentFrequency

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class PDBFileReadWriteTest ( TestCase ):
    """A test case for reading and writing PDB files."""

    def runTest ( self ):
        """The test."""

        # . Initialization.
        isOK         = True
        numberErrors = 0

        # . Output setup.
        sourcePath = os.path.join ( os.getenv ( "PDYNAMO_PBABEL" ), "data", "pdb" )
        outPath    = None
        if self.resultPath is not None:
            outPath = os.path.join ( self.resultPath, "pdb" )
        log = self.GetLog ( )

        # . Set up the output directory.
        if outPath is not None:
            if not os.path.exists ( outPath ): os.mkdir ( outPath )
            outFiles = glob.glob ( os.path.join ( outPath, "*.pdb" ) )
            for outFile in outFiles: os.remove ( outFile )

        # . File names.
        pdbFiles = glob.glob ( os.path.join ( sourcePath, "*.pdb" ) )
        pdbFiles.sort ( )

        # . Loop over the files.
        for pdbFile in pdbFiles:

            # . Process file name.
            ( head, fileName ) = os.path.split ( pdbFile )

            # . Header.
            if log is not None:
                log.Paragraph ( "Processing " + fileName + ":" )

            try:

                # . Reading.
                system = PDBFile_ToSystem ( pdbFile, log = log )
                system.Summary ( log = log )
                if log is not None:
                    PrintComponentFrequency ( system.sequence, log = log )
                    log.LineBreak ( )

                # . Writing.
                if outPath is not None:
                    PDBFile_FromSystem ( os.path.join ( outPath, fileName ), system )

            # . Error.
            except Exception as e:
                numberErrors += 1
                if log is not None: log.Text ( "\nError occurred> " +  e.args[0] + "\n" )

        # . Success/failure.
        self.assertTrue ( isOK and ( numberErrors == 0 ) )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":

    # . Run the test.
    test = PDBFileReadWriteTest ( )
    test.run ( )
