"""MNDO UHF tests."""

import math, os.path

from pCore         import Clone, logFile, LogFileActive, TestCase, TestDataSet, TestReal
from pMolecule     import QCModelMNDO, SystemGeometryObjectiveFunction
from QCTestSystems import GetRadicalMoleculeData, QCTestSystem

#===================================================================================================================================
# . Parameters for test.
#===================================================================================================================================
# . Tolerances.
_BondOrderTolerance             = 0.1
_EnergyTolerance                = 0.1
_GradientAbsoluteErrorTolerance = 1.0e-02

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class MNDOUHFEnergiesTest ( TestCase ):
    """A test case for calculating MNDO UHF energies."""

    def __init__ ( self, *args ):
        """Constructor."""
        super ( MNDOUHFEnergiesTest, self ).__init__ ( *args )
        # . Options.
        self.doPrinting           = True
        self.maximumEnergyAtoms   = 100
        self.maximumGradientAtoms = 100
        self.testGradients        = True

    def ConvergerKeywords ( self ): return { "densityTolerance" : 1.0e-8, "maximumSCFCycles" : 500 }
    def QCModelArguments  ( self ): return ( "am1", )
    def QCModelClass      ( self ): return QCModelMNDO
    def QCModelKeywords   ( self ): return {}

    def runTest ( self ):
        """The test."""

        # . Initialization.
        log          = self.GetLog ( )
        numberErrors = 0
        if self.testGradients:
            maximumGradientDeviation = 0.0

        # . Loop over systems.
        for testSystem in self.testSystems:

            # . Get the molecule.
            molecule = testSystem.GetSystem ( log = log, maximumAtoms = self.maximumEnergyAtoms )
            if molecule is None: continue

            # . Determine an energy.
            try:
                energy = molecule.Energy ( log = log, doGradients = True )

                # . Bond orders.
                labels = []
                for i in range ( len ( molecule.atoms ) ): labels.append ( molecule.atoms[i].path )
                ( bondOrders, charges, freevalence, totalvalence ) = molecule.energyModel.qcModel.BondOrders ( molecule.configuration )
                bondOrders.PrintWithCondition ( ( lambda x, i, j: ( i != j ) and ( math.fabs ( x ) >= _BondOrderTolerance ) ), itemFormat = "{:.3f}", labels = labels, log = log, title = "Bond Orders" )

                # . Charges.
                charges = molecule.AtomicCharges ( )
                charges.Print ( log = log, title = "Charges" )
                if log is not None: log.Paragraph ( "Total Charge = {:.3f}".format ( sum ( charges ) ) )
                charges = molecule.AtomicCharges ( spinDensities = True )
                charges.Print ( log = log, title = "Spin Densities" )
                if log is not None: log.Paragraph ( "Total Spin Density = {:.3f}".format ( sum ( charges ) ) )

                # . Dipole moment.
                dipole = molecule.DipoleMoment ( )
                dipole.Print ( log = log, title = "Dipole" )

                # . Gradient testing.
                if self.testGradients and ( len ( molecule.atoms ) < self.maximumGradientAtoms ):
                    of = SystemGeometryObjectiveFunction.FromSystem ( molecule )
                    gradientDeviation = of.TestGradients ( log = log )
                    maximumGradientDeviation = max ( maximumGradientDeviation, gradientDeviation )

            # . Error.
            except Exception as e:
                numberErrors += 1
                if log is not None: log.Text ( "\nError occurred> " +  e.args[0] + "\n" )

        # . Get the observed and reference data.
        observed      = {}
        referenceData = TestDataSet ( "MNDO UHF Energies" )
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
        self.testSystems = []
        for values in GetRadicalMoleculeData ( ):
            # . Basic keyword arguments.
            kwargs = Clone ( values )
            # . Specific keyword arguments.
            kwargs["convergerKeywords"] = self.ConvergerKeywords ( )
            kwargs["qcModelArguments" ] = self.QCModelArguments  ( )
            kwargs["qcModelClass"     ] = self.QCModelClass      ( )
            qcModelKeywords = kwargs.get ( "qcModelKeywords", {} )
            qcModelKeywords.update ( self.QCModelKeywords ( ) )
            kwargs["qcModelKeywords"  ] = qcModelKeywords
            # . Get the system.
            self.testSystems.append ( QCTestSystem ( **kwargs ) )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":

    # . Run the test.
    test = MNDOUHFEnergiesTest ( )
    test.run ( )
