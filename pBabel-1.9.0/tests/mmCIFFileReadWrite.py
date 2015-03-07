"""Read and write mmCIF files."""

import glob, os, os.path

from pBabel import mmCIFFile_FromSystem, mmCIFFile_ToSystem
from pCore  import logFile, LogFileActive, TestCase

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class mmCIFFileReadWriteTest ( TestCase ):
    """A test case for reading and writing mmcif files."""

    def runTest ( self ):
        """The test."""

        # . Initialization.
        isOK         = True
        numberErrors = 0

        # . Output setup.
        sourcePath = os.path.join ( os.getenv ( "PDYNAMO_PBABEL" ), "data", "mmcif" )
        outPath    = None
        if self.resultPath is not None:
            outPath = os.path.join ( self.resultPath, "mmcif" )
        log = self.GetLog ( )

        # . Set up the output directory.
        if outPath is not None:
            if not os.path.exists ( outPath ): os.mkdir ( outPath )
            outFiles = glob.glob ( os.path.join ( outPath, "*.cif" ) )
            for outFile in outFiles: os.remove ( outFile )

        # . File names.
        cifFiles = glob.glob ( os.path.join ( sourcePath, "*.cif" ) )
        cifFiles.sort ( )

        # . Loop over the files.
        for cifFile in cifFiles:

            # . Process file name.
            ( head, fileName ) = os.path.split ( cifFile )

            # . Header.
            if log is not None:
                log.Paragraph ( "Processing " + fileName + ":" )

            try:

                # . Reading.
                system = mmCIFFile_ToSystem ( cifFile, log = log )
                system.Summary ( log = log )

                # . Writing.
                if outPath is not None:
                    mmCIFFile_FromSystem ( os.path.join ( outPath, fileName ), system )

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
    test = mmCIFFileReadWriteTest ( )
    test.run ( )
