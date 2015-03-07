"""Basic ORCA tests."""

# . Only limited accuracy can be expected in the gradients due to ORCA, not QC/MM.

import math, os.path

from pCore           import logFile, LogFileActive, TestCase, TestDataSet, TestReal
from pMolecule       import NBModelORCA, QCModelORCA, SystemGeometryObjectiveFunction
from QCMMTestSystems import qcmmTestSystems

#===================================================================================================================================
# . Parameters for system definitions.
#===================================================================================================================================
# . QC models.
_qcModels = { "Cyclohexane 9" : [ ( "HF:3-21G"    , "SCFCONV10", "EXTREMESCF" ), True ]  ,
              "Water Dimer 1" : [ ( "MP2:6-31G*"  , "SCFCONV10", "EXTREMESCF" ), True ]  ,
              "Water Dimer 1" : [ ( "HF:3-21G"    , "SCFCONV10", "EXTREMESCF" ),
                                  ( "DFT-ENERGY+" , "SCFCONV10", "EXTREMESCF" ), False ] }

#===================================================================================================================================
# . Parameters for test.
#===================================================================================================================================
# . Tolerances.
_EnergyTolerance                = 0.1
_GradientAbsoluteErrorTolerance = 0.1 # . Not very precise.

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class ORCAEnergiesTest ( TestCase ):
    """A test case for calculating MNDO QC/MM energies."""

    def __init__ ( self, *args ):
        """Constructor."""
        super ( ORCAEnergiesTest, self ).__init__ ( *args )
        # . Options.
        self.doLong        = True
        self.maximumAtoms  = 100
        self.testGradients = True

    def MakeShort ( self ):
        self.doLong        = False
        self.testGradients = False
        return True

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
            for ( doQCMM, qcModel ) in self.qcModels.get ( label, [] ):

                # . Get the molecule.
                molecule = testSystem.GetSystem ( doQCMM = doQCMM, log = log, nbModel = NBModelORCA ( ), qcModel = qcModel )

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
                        gradientDeviation = of.TestGradients ( delta = 5.0e-04, log = log, tolerance = 1.0e-02 )
                        maximumGradientDeviation = max ( maximumGradientDeviation, gradientDeviation )

                # . Error.
                except Exception as e:
                    numberErrors += 1
                    if log is not None: log.Text ( "\nError occurred> " +  e.args[0] + "\n" )

        # . Get the observed and reference data.
        observed      = {}
        referenceData = TestDataSet ( "ORCA Energies" )
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
            qcModelOptions = _qcModels[label]
            doQCMM         = qcModelOptions.pop ( -1 )
            for options in qcModelOptions:
                qcModels.append ( ( doQCMM, QCModelORCA ( *options, useRandomJob = True ) ) )
            self.qcModels[label] = qcModels

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":

    # . Run the test.
    test = ORCAEnergiesTest ( )
    test.run ( )
