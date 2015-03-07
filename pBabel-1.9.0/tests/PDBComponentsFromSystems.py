"""Interconvert PDB components and systems."""

import glob, os, os.path

from pBabel import MOLFile_ToSystem, PDBComponent
from pCore  import logFile, LogFileActive, TestCase

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class PDBComponentSystemTest ( TestCase ):
    """A test case for reading EMSL basis set files in Gaussian 94 format."""

    def runTest ( self ):
        """The test."""

        # . Initialization.
        isOK         = True
        numberErrors = 0

        # . Output setup.
        dataPath = os.path.join ( os.getenv ( "PDYNAMO_PMOLECULE" ), "data", "mol" )
        log = self.GetLog ( )

        # . Get the files.
        molFiles = glob.glob ( os.path.join ( dataPath, "*.mol" ) )
        molFiles.sort ( )

        # . Loop over the files.
        for ( i, molFile ) in enumerate ( molFiles ):

            try:

                # . Get the system.
                system = MOLFile_ToSystem ( molFile )
                system.Summary ( log = log )

                # . Convert to PDB component.
                component = PDBComponent.FromSystem ( system, label = ( "{:03d}".format ( i ) ) )
                component.Summary ( log = log )

                # . Convert back to a system.
                system = component.ToSystem ( )
                system.Summary ( log = log )

                # . Finish up.
                log.Separator ( )

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
    test = PDBComponentSystemTest ( )
    test.run ( )
