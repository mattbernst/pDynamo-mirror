"""Test that rotation and translation is correctly removed during a dynamics simulation."""

import glob, os

from pBabel           import MOLFile_ToSystem
from pCore            import Clone, NormalDeviateGenerator, RandomNumberGenerator, TestCase, TestDataSet, TestReal
from pMolecule        import MMModelOPLS, NBModelFull
from pMoleculeScripts import LangevinDynamics_SystemGeometry

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Dynamics options.
_NSteps = 20000

# . Tolerances.
_RMSAbsoluteErrorTolerance = 1.0e-03

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class RotationTranslationRemovalTest ( TestCase ):
    """Check rotational and translational motion removal during a dynamics simulation."""

    def runTest ( self ):
        """The test."""

        # . Paths.
        dataPath = os.path.join ( os.getenv ( "PDYNAMO_PMOLECULE" ), "data", "mol" )
        log = self.GetLog ( )

        # . Get the system.
        molecule = MOLFile_ToSystem ( os.path.join ( dataPath, "tyrosineDipeptide.mol" ) )
        molecule.DefineMMModel ( MMModelOPLS ( "protein" ) )
        molecule.DefineNBModel ( NBModelFull ( ) )
        molecule.Summary ( log = log )
        molecule.Energy  ( log = log, doGradients = True )

        # . Save initial coordinates.
        reference3 = Clone ( molecule.coordinates3 )

        # . Do some dynamics.
        normalDeviateGenerator = NormalDeviateGenerator.WithRandomNumberGenerator ( RandomNumberGenerator.WithSeed ( 247171 ) )
        LangevinDynamics_SystemGeometry ( molecule                         ,
                                          collisionFrequency     =    25.0 ,
                                          log                    =     log ,
                                          logFrequency           =    1000 ,
                                          normalDeviateGenerator = normalDeviateGenerator ,
                                          steps                  = _NSteps ,
                                          temperature            =   300.0 ,
                                          timeStep               =   0.001 )

        # . Check RMSs which should be the same as rotation and translation are removed.
        masses = molecule.atoms.GetItemAttributes ( "mass" )
        rms0   = molecule.coordinates3.RMSDeviation ( reference3, weights = masses )
        molecule.coordinates3.Superimpose ( reference3, weights = masses )
        rms1 = molecule.coordinates3.RMSDeviation ( reference3, weights = masses )

        # . Get the observed and reference data.
        observed      = { "RMS Deviation" : rms1 }
        referenceData = TestDataSet ( "Rotation/Translation Removal" )
        referenceData.AddDatum ( TestReal ( "RMS Deviation", rms0, referenceData, absoluteErrorTolerance = _RMSAbsoluteErrorTolerance, toleranceFormat = "{:.3f}", valueFormat = "{:.3f}" ) )

        # . Check for success/failure.
        if len ( observed ) > 0:
            results = referenceData.VerifyAgainst ( observed )
            results.Summary ( log = log, fullSummary = self.fullVerificationSummary )
            isOK    = results.WasSuccessful ( )
        else:
            isOK    = True
        self.assertTrue ( isOK )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":

    # . Run the test.
    test = RotationTranslationRemovalTest ( )
    test.run ( )
