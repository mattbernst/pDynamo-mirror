"""MNDO CI QC/MM tests."""

import math, os.path

from pCore           import logFile, LogFileActive, TestCase, TestDataSet, TestReal
from pMolecule       import DIISSCFConverger, QCModelMNDO, SystemGeometryObjectiveFunction
from QCMMTestSystems import qcmmTestSystems

#===================================================================================================================================
# . Parameters for system definitions.
#===================================================================================================================================
# . QC models.
_qcModelKeywordLabels = ( "CIMethod", "activeElectrons", "activeOrbitals", "occupancyType", "requiredRoot", "rootMultiplicity" )
_qcModels = { "Alanine Dipeptide"    : ( ( "Singles/Doubles",  8, 6, "Cardinal"        , 0, 1 ), ), \
              "Cyclohexane 6"        : ( ( "Singles/Doubles",  4, 4, "Cardinal"        , 0, 1 ), ), \
              "Cyclohexane 9"        : ( ( "Singles/Doubles",  6, 6, "Cardinal"        , 1, 1 ), ), \
              "Tyrosine Dipeptide 1" : ( ( "Singles/Doubles",  6, 6, "Cardinal"        , 0, 3 ), ), \
              "Tyrosine Dipeptide 2" : ( ( "Singles/Doubles",  5, 6, "Fixed Fractional", 0, 2 ), ), \
              "Water Dimer 1"        : ( ( "Singles/Doubles",  8, 6, "Cardinal"        , 1, 1 ), ), \
              "Water Dimer 2"        : ( ( "Singles/Doubles",  8, 6, "Cardinal"        , 0, 3 ), )  }

#===================================================================================================================================
# . Parameters for test.
#===================================================================================================================================
# . Tolerances.
_EnergyTolerance                = 0.1
_GradientAbsoluteErrorTolerance = 1.0e-02

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class MNDOCIEnergiesQCMMTest ( TestCase ):
    """A test case for calculating MNDO CI QC/MM energies."""

    def __init__ ( self, *args ):
        """Constructor."""
        super ( MNDOCIEnergiesQCMMTest, self ).__init__ ( *args )
        # . Options.
        self.doPrinting    = True
        self.maximumAtoms  = 100
        self.testGradients = True

    def runTest ( self ):
        """The test."""

        # . Initialization.
        log          = self.GetLog ( )
        numberErrors = 0
        if self.testGradients:
            maximumGradientDeviation = 0.0

        # . Loop over systems and QC models.
        labels = self.qcmmTestSystems.keys ( )
        labels.sort ( )
        for label in labels:
            testSystem = self.qcmmTestSystems[label]
            for qcModel in self.qcModels[label]:

                # . Get the molecule.
                molecule = testSystem.GetSystem ( log = log, qcModel = qcModel )

                # . Energy.
                try:
                    energy  = molecule.Energy ( log = log, doGradients = True )

                    # . Charges.
                    charges = molecule.AtomicCharges ( )
                    charges.Print ( log = log, title = "Charges" )
                    if log is not None: log.Paragraph ( "Total Charge = {:.3f}".format ( sum ( charges ) ) )
                    charges = molecule.AtomicCharges ( spinDensities = True )
                    charges.Print ( log = log, title = "Spin Densities" )
                    if log is not None: log.Paragraph ( "Total Spin Density = {:.3f}".format ( sum ( charges ) ) )

                    # . Gradient testing.
                    if self.testGradients and ( len ( molecule.atoms ) < self.maximumAtoms ):
                        of = SystemGeometryObjectiveFunction.FromSystem ( molecule )
                        gradientDeviation = of.TestGradients ( delta = 1.0e-06, log = log, tolerance = 1.0e-03 )
                        maximumGradientDeviation = max ( maximumGradientDeviation, gradientDeviation )

                # . Error.
                except Exception as e:
                    numberErrors += 1
                    if log is not None: log.Text ( "\nError occurred> " +  e.args[0] + "\n" )

        # . Get the observed and reference data.
        observed      = {}
        referenceData = TestDataSet ( "MNDO CI QC/MM Energies" )
        if self.testGradients:
            observed["Gradient Error"] = maximumGradientDeviation
            referenceData.AddDatum ( TestReal ( "Gradient Error", 0.0, referenceData, absoluteErrorTolerance = _GradientAbsoluteErrorTolerance, toleranceFormat = "{:.3f}", valueFormat = "{:.3f}" ) )

        # . Check for success/failure.
        if len ( observed ) > 0:
            results = referenceData.VerifyAgainst ( observed )
            results.Summary ( log = log, fullSummary = self.fullVerificationSummary )
            isOK    = results.WasSuccessful ( )
        else:
            isOK    = True
        isOK = isOK and ( numberErrors == 0 )
        self.assertTrue ( isOK )

    def setUp ( self ):
        """Set up the calculation."""
        # . Set up the systems.
        self.qcmmTestSystems = qcmmTestSystems
        # . Define the QC models for each system.
        self.qcModels = {}
        for label in _qcModels:
            qcModels = []
            for qcModelOptions in _qcModels[label]:
                kwargs                    = { key : value for ( key, value ) in zip ( _qcModelKeywordLabels, qcModelOptions ) }
                kwargs["converger"      ] = DIISSCFConverger ( densityTolerance = 1.0e-11, maximumSCFCycles = 250 )
                kwargs["keepOrbitalData"] = True
                qcModels.append ( QCModelMNDO ( "am1", **kwargs ) )
            self.qcModels[label] = qcModels

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":

    # . Run the test.
    test = MNDOCIEnergiesQCMMTest ( )
    test.run ( )
