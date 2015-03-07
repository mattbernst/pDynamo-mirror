"""Find CIP labels for several bALA configurations."""

import glob, os

from pBabel           import MOLFile_ToSystem, XYZFile_ToCoordinates3
from pCore            import TestCase
from pMoleculeScripts import CIPLabelFinder

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class CIPLabelTest ( TestCase ):
    """Find CIP labels."""

    def runTest ( self ):
        """The test."""

        # . Paths.
        dataPath = os.path.join ( os.getenv ( "PDYNAMO_PMOLECULE" ), "data" )
        molPath  = os.path.join ( dataPath, "mol" )
        xyzPath  = os.path.join ( dataPath, "bAlaConformations" )
        log = self.GetLog ( )

        # . Conformations.
        xyzFiles = glob.glob ( os.path.join ( xyzPath, "*.xyz" ) )
        xyzFiles.sort ( )

        # . Generate the molecule.
        molecule = MOLFile_ToSystem ( os.path.join ( molPath, "bAla_c7eq.mol" ) )
        molecule.Summary ( log = log )

        # . Initialization.
        isOK = True

        # . Loop over the structures in the xyz files.
        for xyzFile in xyzFiles:
            molecule.coordinates3 = XYZFile_ToCoordinates3 ( xyzFile )
            if log is not None:
                conformation = os.path.split ( xyzFile )[-1].split ( "_", 1 )[-1][0:-4]
                log.Heading ( "bALA Configuration " + conformation, QBLANKLINE = True )
            results = CIPLabelFinder ( molecule, log = log )
            if results is None:
                localIsOK = False
            else:
                ( ( tCenters, rtCenters, stCenters, utCenters ), ( dCenters, edCenters, zdCenters, udCenters ) ) = results
                localIsOK = ( len ( tCenters  ) == 4 ) and ( len ( dCenters  ) == 0 ) and \
                            ( len ( stCenters ) == 1 ) and ( len ( utCenters ) == 3 )
            isOK = ( isOK and localIsOK )

        # . Success/failure.
        self.assertTrue ( isOK )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":

    # . Run the test.
    test = CIPLabelTest ( )
    test.run ( )
