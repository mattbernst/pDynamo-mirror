"""Test for hydrogen building."""

import glob, os

from pBabel           import MOLFile_ToSystem
from pCore            import RandomNumberGenerator, Selection, TestCase
from pMolecule        import MMModelOPLS, NBModelFull
from pMoleculeScripts import BuildHydrogenCoordinates3FromConnectivity

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class HydrogenBuildingTest ( TestCase ):
    """Test hydrogen building."""

    def runTest ( self ):
        """The test."""

        # . Paths.
        dataPath = os.path.join ( os.getenv ( "PDYNAMO_PMOLECULE" ), "data", "mol" )
        log = self.GetLog ( )

        # . Energy models.
        mmModel = MMModelOPLS ( "protein" )
        nbModel = NBModelFull ( )

        # . Get the files.
        molFiles = glob.glob ( os.path.join ( dataPath, "*.mol" ) )

        # . Initialization.
        numberErrors   = 0
        numberFailures = 0

        # . Loop over the files.
        for molFile in molFiles:

            try:

                # . Get the system.
                try:
                    system = MOLFile_ToSystem ( molFile )
                    system.DefineMMModel ( mmModel, log = log )
                    system.DefineNBModel ( nbModel )
                    system.Summary ( log = log )
                except:
                    continue

                # . Calculate an energy.
                eBefore = system.Energy ( log = log, doGradients = True )

                # . Define all hydrogen positions as undefined.
                for ( i, atom ) in enumerate ( system.atoms ):
                    if atom.atomicNumber == 1:
                        system.coordinates3.FlagCoordinateAsUndefined ( i )

                # . Build as many undefined coordinates as possible.
                randomNumberGenerator = RandomNumberGenerator.WithSeed ( 957197 )
                BuildHydrogenCoordinates3FromConnectivity ( system, log = log, randomNumberGenerator = randomNumberGenerator )

                # . Calculate an energy if all coordinates have been defined.
                if system.coordinates3.numberUndefined > 0:
                    numberFailures += 1
                    if log is not None: log.Paragraph ( "Not all hydrogens have been rebuilt." )
                else:
                    eAfter = system.Energy ( log = log, doGradients = True )
                    if log is not None: log.Paragraph ( "Energy difference after rebuilding = {:.1f}.".format ( eAfter - eBefore ) )

            except Exception as e:
                numberErrors += 1
                if log is not None: log.Text ( "\nError occurred> " +  e.args[0] + "\n" )

        # . Success/failure.
        self.assertTrue ( ( ( numberErrors == 0 ) and ( numberFailures == 0 ) ) )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":

    # . Run the test.
    test = HydrogenBuildingTest ( )
    test.run ( )
