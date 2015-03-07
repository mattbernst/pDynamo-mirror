"""MNDO QC/MM tests."""

import math, os.path

from pCore           import logFile, LogFileActive, TestCase, TestDataSet, TestReal
from pMolecule       import DIISSCFConverger, QCModelMNDO, SystemGeometryObjectiveFunction
from QCMMTestSystems import qcmmTestSystems

#===================================================================================================================================
# . Parameters for test.
#===================================================================================================================================
# . Tolerances.
_EnergyTolerance                = 0.1
_GradientAbsoluteErrorTolerance = 1.0e-03

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class MNDOQCMMEnergiesTest ( TestCase ):
    """A test case for calculating MNDO QC/MM energies."""

    def __init__ ( self, *args ):
        """Constructor."""
        super ( MNDOQCMMEnergiesTest, self ).__init__ ( *args )
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
            if testSystem.multiplicity != 1: continue

            # . Get the QC model.
            converger = DIISSCFConverger ( densityTolerance = 1.0e-10, maximumSCFCycles = 250 )
            qcModel   = QCModelMNDO ( "am1", converger = converger, isSpinRestricted = ( testSystem.multiplicity == 1 ) )

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
                    gradientDeviation = of.TestGradients ( delta = 1.0e-04, log = log, tolerance = 1.0e-03 )
                    maximumGradientDeviation = max ( maximumGradientDeviation, gradientDeviation )

            # . Error.
            except Exception as e:
                numberErrors += 1
                if log is not None: log.Text ( "\nError occurred> " +  e.args[0] + "\n" )

        # . Get the observed and reference data.
        observed      = {}
        referenceData = TestDataSet ( "MNDO QC/MM Energies" )
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
        self.qcmmTestSystems = qcmmTestSystems

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":

    # . Run the test.
    test = MNDOQCMMEnergiesTest ( )
    test.run ( )
