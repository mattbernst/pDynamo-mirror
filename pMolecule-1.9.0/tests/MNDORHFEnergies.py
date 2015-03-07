"""MNDO RHF tests."""

import math, os.path

from pCore         import Clone, logFile, LogFileActive, Pickle, Statistics, TestCase, TestDataSet, TestReal, Unpickle
from pMolecule     import QCModelMNDO, SystemGeometryObjectiveFunction
from QCTestSystems import GetClosedShellMoleculeData, QCTestSystem

#===================================================================================================================================
# . Parameters for test.
#===================================================================================================================================
# . Parameters.
_TableModelWidth = 20

# . Tolerances.
_AbsoluteErrorTolerance         = 1.0
_BondOrderTolerance             = 0.1
_EnergyTolerance                = 0.1
_GradientAbsoluteErrorTolerance = 2.0e-02

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class MNDORHFEnergiesTest ( TestCase ):
    """A test case for calculating MNDO RHF energies."""

    def __init__ ( self, *args ):
        """Constructor."""
        super ( MNDORHFEnergiesTest, self ).__init__ ( *args )
        # . Options.
        self.doPrinting           = True
        self.maximumEnergyAtoms   =   30
        self.maximumEnergyTests   = 1000
        self.maximumGradientAtoms =   10
        self.maximumGradientTests =   20
        self.testGradients        = True

    def ConvergerKeywords ( self ):
        return { "densityTolerance" : 1.0e-8, "maximumSCFCycles" : 500 }

    def GenerateReferenceData ( self ):
        self.generateReferenceData = True
        self.testGradients         = False
        return True

    def QCModelClass   ( self ): return QCModelMNDO

    def QCModelOptions ( self ):
        self.modelLabels = ( "am1", "mndo", "pddgmndo", "pddgpm3", "pm3", "pm6", "rm1" )
        return ( ( ( "am1", ), {} ), ( ( "mndo", ), {} ), ( ( "pddgmndo", ), {} ), ( ( "pddgpm3", ), {} ), ( ( "pm3", ), {} ), ( ( "pm6", ), {} ), ( ( "rm1", ), {} ) )

    def runTest ( self ):
        """The test."""

        # . Initialization.
        log                 = self.GetLog ( )
        modelResults        = {}
        numberEnergyTests   = 0
        numberErrors        = 0
        numberGradientTests = 0
        if self.testGradients:
            maximumGradientDeviation = 0.0

        # . Energy data.
        if self.generateReferenceData:
            referenceData  = TestDataSet ( "MNDO RHF Energies" )
            energyDataSets = {}
            for model in self.modelLabels:
                energyDataSet         = TestDataSet ( model, parent = referenceData )
                energyDataSets[model] = energyDataSet
                referenceData.AddDatum ( energyDataSet )
        else:
            energyResults = {}

        # . Loop over systems.
        for testSystem in self.testSystems:

            # . Get results (if necessary).
            modelLabel = testSystem.modelLabel
            if not self.generateReferenceData:
                energies = energyResults.get ( modelLabel, None )
                if energies is None:
                    energies = {}
                    energyResults[modelLabel] = energies
            results = modelResults.get ( modelLabel, None )
            if results is None:
                results = {}
                modelResults[modelLabel] = results
            cycles  = results.get ( "Cycles", None )
            if cycles is None:
                cycles = []
                results["Cycles"] = cycles

            # . Get the molecule.
            try:
                molecule = testSystem.GetSystem ( log = log, maximumAtoms = self.maximumEnergyAtoms )
            except:
                molecule = None
                results["Undefined"] = results.get ( "Undefined", 0 ) + 1
            if molecule is None: continue

            # . Determine an energy.
            try:
                try:
                    e = molecule.Energy ( log = log, doGradients = True )
                    i = molecule.configuration.qcState.GetItem ( "SCF Cycles" )
                    results["Converged"]     = results.get ( "Converged", 0 ) + 1
                    cycles.append ( float ( i ) )
                    # . Generate reference data - converged only.
                    if self.generateReferenceData:
                        energyDataSet = energyDataSets[modelLabel]
                        energyDataSet.AddDatum ( TestReal ( molecule.label, e, energyDataSet, absoluteErrorTolerance = _AbsoluteErrorTolerance, toleranceFormat = "{:.3f}", valueFormat = "{:.3f}" ) )
                    # . Accumulate observed data.
                    else: energies[molecule.label] = e
                except Exception as e:
                    results["Not Converged"] = results.get ( "Not Converged", 0 ) + 1

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
                if self.testGradients and ( len ( molecule.atoms ) < self.maximumGradientAtoms ) and ( numberGradientTests < self.maximumGradientTests ):
                    of = SystemGeometryObjectiveFunction.FromSystem ( molecule )
                    gradientDeviation = of.TestGradients ( log = log )
                    maximumGradientDeviation = max ( maximumGradientDeviation, gradientDeviation )
                    numberGradientTests += 1

            # . Error.
            except Exception as e:
                numberErrors += 1
                if log is not None: log.Text ( "\nError occurred> " +  e.args[0] + "\n" )

            # . Check the number of tests.
            numberEnergyTests += 1
            if numberEnergyTests >= self.maximumEnergyTests: break

        # . Print model results.
        if log is not None:
            models = modelResults.keys ( )
            models.sort ( )
            table = log.GetTable ( columns = 7 * [ _TableModelWidth ] )
            table.Start   ( )
            table.Title   ( "Model Results" )
            table.Heading ( "Model"         )
            table.Heading ( "Undefined"     )
            table.Heading ( "Not Converged" )
            table.Heading ( "Converged"     )
            table.Heading ( "<Cycles>"      )
            table.Heading ( "Max. Cycles"   )
            table.Heading ( "Min. Cycles"   )
            for model in models:
                results = modelResults[model]
                cycles  = Statistics ( results["Cycles"] )
                table.Entry ( model, alignment = "left" )
                table.Entry ( "{:d}".format ( results.get ( "Undefined"    , 0 ) ) )
                table.Entry ( "{:d}".format ( results.get ( "Not Converged", 0 ) ) )
                table.Entry ( "{:d}".format ( results.get ( "Converged"    , 0 ) ) )
                if cycles.size > 0:
                    table.Entry ( "{:.1f}".format ( cycles.mean    ) )
                    table.Entry ( "{:.0f}".format ( cycles.maximum ) )
                    table.Entry ( "{:.0f}".format ( cycles.minimum ) )
                else:
                    table.Entry ( "-" )
                    table.Entry ( "-" )
                    table.Entry ( "-" )
            table.Stop ( )

        # . Energy data.
        # . Generate the reference data.
        if self.generateReferenceData:
            referenceData.Summary ( log = log )
            Pickle ( self.referenceDataPath, referenceData )
            isOK = True
        # . Verify the observed data against the reference data.
        else:
            referenceData = Unpickle ( self.referenceDataPath )
            results = referenceData.VerifyAgainst ( energyResults )
            results.Summary ( log = log, fullSummary = self.fullVerificationSummary )
            isOK    = results.WasSuccessful ( )

        # . Gradient data set.
        if self.testGradients:

            # . Get the observed and reference data.
            observed      = {}
            referenceData = TestDataSet ( "MNDO RHF Gradients" )
            observed["Gradient Error"] = maximumGradientDeviation
            referenceData.AddDatum ( TestReal ( "Gradient Error", 0.0, referenceData, absoluteErrorTolerance = _GradientAbsoluteErrorTolerance, toleranceFormat = "{:.3f}", valueFormat = "{:.3f}" ) )

            # . Check for success/failure.
            if len ( observed ) > 0:
                results = referenceData.VerifyAgainst ( observed )
                results.Summary ( log = log, fullSummary = self.fullVerificationSummary )
                isOK    = isOK and results.WasSuccessful ( )

        # . Finish up.
        isOK = isOK and ( numberErrors == 0 )
        self.assertTrue ( isOK )

    def setUp ( self ):
        """Set up the calculation."""
        # . Get basic options.
        convergerKeywords = self.ConvergerKeywords ( )
        qcModelClass      = self.QCModelClass      ( )
        qcModelOptions    = self.QCModelOptions    ( )
        # . Set up the systems for each set of options.
        self.testSystems  = []
        for values in GetClosedShellMoleculeData ( ):
            for ( i, ( arguments, keywords ) ) in enumerate ( qcModelOptions ):
                # . Basic keyword arguments.
                kwargs = Clone ( values )
                # . Specific keyword arguments.
                kwargs["convergerKeywords"] = convergerKeywords
                kwargs["qcModelArguments" ] = arguments
                kwargs["qcModelClass"     ] = qcModelClass
                qcModelKeywords = kwargs.get ( "qcModelKeywords", {} )
                qcModelKeywords.update ( keywords )
                kwargs["qcModelKeywords"  ] = qcModelKeywords
                kwargs["modelLabel"       ] = self.modelLabels[i]
                # . Get the system.
                self.testSystems.append ( QCTestSystem ( **kwargs ) )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":

    # . Run the test.
    test = MNDORHFEnergiesTest ( )
    test.run ( )
