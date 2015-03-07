"""Simple ONIOM tests with either a QC(MNDO) or a QC(MNDO)/MM lower layer and a QC(ab initio) upper layer."""

import math, os.path

from pCore           import logFile, LogFileActive, Selection, TestCase, TestDataSet, TestReal
from pMolecule       import DIISSCFConverger, ElectronicState, MultiLayerSystemGeometryObjectiveFunction, NBModelABFS, QCModelDFT, QCModelMNDO, QCModelORCA
from QCMMTestSystems import qcmmTestSystems

#===================================================================================================================================
# . Parameters for system definitions.
#===================================================================================================================================
# . Model definitions.
_converger   = DIISSCFConverger ( densityTolerance = 1.0e-8, maximumSCFCycles = 250 )
_qcModelDFT  = QCModelDFT  ( converger = _converger, densityBasis = "demon", orbitalBasis = "321g" )
_qcModelMNDO = QCModelMNDO ( converger = _converger )
_qcModelORCA = QCModelORCA ( "HF:6-31G*", "SCFCONV10", "EXTREMESCF" )
_nbModel     = NBModelABFS ( )

# . Job data.
# . Options: doQCMM for bottom layer, QC model for upper layer and QC selection for upper layer.
_jobData = { "Water Dimer 1"     : ( ( False, _qcModelDFT,  Selection.FromIterable ( range (  3     ) ) ), \
                                     ( False, _qcModelORCA, Selection.FromIterable ( range (  3,  6 ) ) ), \
                                     ( True,  _qcModelDFT,  Selection.FromIterable ( range (  3     ) ) )  ), \
             "Water Dimer 2"     : ( ( True,  _qcModelORCA, Selection.FromIterable ( range (  3,  6 ) ) ), ), \
             "Cyclohexane 9"     : ( ( False, _qcModelDFT,  Selection.FromIterable ( range (  6     ) ) ), \
                                     ( True,  _qcModelORCA, Selection.FromIterable ( range (  3,  6 ) ) )  ), \
             "Alanine Dipeptide" : ( ( False, _qcModelDFT,  Selection.FromIterable ( range ( 10, 14 ) ) ), \
                                     ( False, _qcModelORCA, Selection.FromIterable ( range ( 10, 14 ) ) ), \
                                     ( True,  _qcModelDFT,  Selection.FromIterable ( range ( 10, 14 ) ) ), \
                                     ( True,  _qcModelORCA, Selection.FromIterable ( range ( 10, 14 ) ) )  ) }

#===================================================================================================================================
# . Parameters for test.
#===================================================================================================================================
# . Tolerances.
_EnergyTolerance                = 0.1
_GradientAbsoluteErrorTolerance = 1.0e-02

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class ONIOMEnergiesTest ( TestCase ):
    """A test case for calculating MNDO QC/MM energies."""

    def __init__ ( self, *args ):
        """Constructor."""
        super ( ONIOMEnergiesTest, self ).__init__ ( *args )
        # . Options.
        self.doLong        = True
        self.maximumAtoms  = 100
        self.skipORCATests = False
        self.testGradients = True

    def MakeShort ( self ):
        self.doLong        = False
        self.maximumAtoms  = 10
        self.skipORCATests = True
        return True

    def runTest ( self ):
        """The test."""

        # . Initialization.
        log          = self.GetLog ( )
        numberErrors = 0
        if self.testGradients:
            maximumGradientDeviation = 0.0

        # . Loop over systems and jobs.
        labels = self.jobs.keys ( )
        labels.sort ( )
        for label in labels:
            testSystem = self.qcmmTestSystems[label]
            for ( doQCMM, qcModelU, qcSelectionU ) in self.jobs.get ( label, [] ):

                # . Check if ORCA tests are to be skipped.
                if self.skipORCATests and ( qcModelU is _qcModelORCA ): continue

                # . Get the molecule.
                molecule = testSystem.GetSystem ( doQCMM = doQCMM, log = log, nbModel = _nbModel, qcModel = _qcModelMNDO )

                # . Energy.
                try:
                    energy  = molecule.Energy ( log = log, doGradients = True )

                    # . Define the object function.
                    of = MultiLayerSystemGeometryObjectiveFunction.FromSystem ( molecule )

                    # . First layer.
                    of.DefineQCLayer ( qcSelectionU, qcModelU, electronicState = ElectronicState ( charge = 0 ) )
                    of.SubsystemSummary ( log = log )

                    # . Gradient testing.
                    if self.testGradients and ( len ( molecule.atoms ) < self.maximumAtoms ):
                        gradientDeviation = of.TestGradients ( delta = 5.0e-04, log = log, tolerance = 1.0e-02 )
                        maximumGradientDeviation = max ( maximumGradientDeviation, gradientDeviation )

                # . Error.
                except Exception as e:
                    numberErrors += 1
                    if log is not None: log.Text ( "\nError occurred> " +  e.args[0] + "\n" )

                # . Clear up.
                finally:
                    if qcModelU is _qcModelORCA: _qcModelORCA.DeleteJobFiles ( )

        # . Get the observed and reference data.
        observed      = {}
        referenceData = TestDataSet ( "ONIOM Energies" )
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
        self.jobs            = _jobData
        self.qcmmTestSystems = qcmmTestSystems

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":

    # . Run the test.
    test = ONIOMEnergiesTest ( )
    test.run ( )
