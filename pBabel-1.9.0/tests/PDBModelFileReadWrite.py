"""Read and write PDB model files."""

import glob, os, os.path

from pBabel           import PDBFile_ToPDBModel, PDBModel_FromModelFile, PDBModel_ToModelFile
from pCore            import logFile, LogFileActive, TestCase
from pMoleculeScripts import PrintComponentFrequency

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class PDBModelFileReadWriteTest ( TestCase ):
    """A test case for reading and writing PDB model files."""

    def runTest ( self ):
        """The test."""

        # . Initialization.
        isOK         = True
        numberErrors = 0

        # . Output setup.
        sourcePath = os.path.join ( os.getenv ( "PDYNAMO_PBABEL" ), "data", "pdb" )
        outPath    = None
        if self.resultPath is not None:
            outPath = os.path.join ( self.resultPath, "pdbModel" )
        log = self.GetLog ( )

        # . Set up the output directory.
        if outPath is not None:
            if not os.path.exists ( outPath ): os.mkdir ( outPath )
            outFiles = glob.glob ( os.path.join ( outPath, "*.model" ) )
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
                model = PDBFile_ToPDBModel ( pdbFile, log = log )
                if log is not None:
                    log.LineBreak ( )
                    model.Summary ( log = log )

                # . Writing.
                if outPath is not None:
                    modelPath = os.path.join ( outPath, fileName[0:-4] + ".model" )
                    PDBModel_ToModelFile ( modelPath, model )

                    # . Reading again.
                    model = PDBModel_FromModelFile ( modelPath )
                    model.Summary ( log = log )

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
    test = PDBModelFileReadWriteTest ( )
    test.run ( )

