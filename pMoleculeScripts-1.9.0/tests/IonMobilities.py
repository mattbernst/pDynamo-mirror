"""Testing for ion mobilities."""

import glob, math, os

from pBabel           import XYZFile_ToSystem
from pCore            import logFile, LogFileActive, Pickle, RandomNumberGenerator, TestCase, TestDataSet, TestReal, Unpickle
from pMoleculeScripts import HardSphereIonMobilities

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Tolerance for acceptability (10%).
_percentErrorTolerance = 10.0

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class IonMobilityTest ( TestCase ):
    """An ion mobility test case."""

    def GenerateReferenceData ( self ):
        self.generateReferenceData = True
        return True

    def runTest ( self ):
        """The test."""

        # . Output setup.
        dataPath = os.path.join ( os.getenv ( "PDYNAMO_PMOLECULE" ), "data" )
        log = self.GetLog ( )

        # . Define the molecule.
        molecule = XYZFile_ToSystem ( os.path.join ( dataPath, "xyz", "serineOctamerZD4L4.xyz" ) )
        molecule.Summary ( log = log )

        # . Define the randomNumberGenerator so as to have reproducible results.
        randomNumberGenerator = RandomNumberGenerator.WithSeed ( 314159 )

        # . Do the test - 250000+ trajectories are generally necessary.
        observed = HardSphereIonMobilities ( molecule, nreflections = 30, ntrajectories = 100000, log = log, randomNumberGenerator = randomNumberGenerator )

        # . Remove non-real values.
        keys = observed.keys ( )
        for key in keys:
            if not isinstance ( observed[key], float ): del observed[key]

        # . Generate the reference data.
        if self.generateReferenceData:
            referenceData = TestDataSet ( "Ion Mobilities" )
            for ( key, value ) in observed.iteritems ( ):
                referenceData.AddDatum ( TestReal ( key, value, referenceData, percentErrorTolerance = _percentErrorTolerance ) )
            referenceData.Summary ( log = log )
            Pickle ( self.referenceDataPath, referenceData )
            isOK = True
        # . Verify the observed data against the reference data.
        else:
            referenceData = Unpickle ( self.referenceDataPath )
            results = referenceData.VerifyAgainst ( observed )
            results.Summary ( log = log, fullSummary = self.fullVerificationSummary )
            isOK    = results.WasSuccessful ( )

        # . Success/failure.
        self.assertTrue ( isOK )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":

    # . Run the test.
    test = IonMobilityTest ( )
    test.run ( )
