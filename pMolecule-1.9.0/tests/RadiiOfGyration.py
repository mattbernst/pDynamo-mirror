"""Radii of gyration of various bALA structures."""

import glob, os

from pBabel import MOLFile_ToSystem, XYZFile_ToCoordinates3
from pCore  import TestCase

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class RadiiOfGyrationTest ( TestCase ):
    """Determine radii of gyration."""

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

        # . Translate to principal axes using masses as weights.
        masses = molecule.atoms.GetItemAttributes ( "mass" )
        molecule.coordinates3.ToPrincipalAxes ( weights = masses )

        # . Loop over structures and output at the same time.
        if log is not None:
            table = log.GetTable ( columns = [ 20, 10 ] )
            table.Start ( )
            table.Title ( "Radii Of Gyration" )
            table.Heading ( "Conformation" )
            table.Heading ( "Value"        )
            for xyzFile in xyzFiles:
                conformation     = os.path.split ( xyzFile )[-1].split ( "_", 1 )[-1][0:-4]
                coordinates3     = XYZFile_ToCoordinates3 ( xyzFile )
                radiusOfGyration = coordinates3.RadiusOfGyration ( weights = masses )
                table.Entry ( conformation, alignment = "left" )
                table.Entry ( "{:.2f}".format ( radiusOfGyration ) )
            table.Stop ( )

        # . Success/failure.
        self.assertTrue ( True )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":

    # . Run the test.
    test = RadiiOfGyrationTest ( )
    test.run ( )
