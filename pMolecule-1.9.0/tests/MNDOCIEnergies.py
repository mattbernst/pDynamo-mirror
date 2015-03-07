"""MNDO CI tests."""

import math, os.path

from pBabel    import XYZFile_ToSystem
from pCore     import logFile, LogFileActive, TestCase, TestDataSet, TestReal
from pMolecule import DIISSCFConverger, ElectronicState, QCModelMNDO, SystemGeometryObjectiveFunction

#===================================================================================================================================
# . Parameters for system definitions.
#===================================================================================================================================
# . Paths.
_dataPath = os.path.join ( os.getenv ( "PDYNAMO_PMOLECULE" ), "data" )

# . User-specified sets of micro-states.
_AllylSet1 = ( "1111000011100000", "0111100011100000", "0111010011100000", "0111001011100000", "0111000111100000", "1011100011100000", "1011010011100000", "1011001011100000", \
               "1011000111100000", "1101100011100000", "1101010011100000", "1101001011100000", "1101000111100000", "1110100011100000", "1110010011100000", "1110001011100000", \
               "1110000111100000", "1111000001110000", "1111000001101000", "1111000001100100", "1111000001100010", "1111000001100001", "1111000010110000", "1111000010101000", \
               "1111000010100100", "1111000010100010", "1111000010100001", "1111000011010000", "1111000011001000", "1111000011000100", "1111000011000010", "1111000011000001", \
               "1110100001110000", "1110010001110000", "1110001001110000", "1110000101110000", "1110100010110000", "1110010010110000", "1110001010110000", "1110000110110000", \
               "1110100011010000", "1110010011010000", "1110001011010000", "1110000111010000" )
_AllylSet2 = ( "110100", "101100", "011100", "110010", "101010", "011010", "110001", "101001", "011001" )

# . Molecule data.
_keywordLabels        = ( "label", "directory", "charge", "multiplicity", "qcModelOptions" )
_qcModelKeywordLabels = ( "CIMethod", "activeElectrons", "activeOrbitals", "numberFractionalHOOs", "numberFractionalLUOs", "occupancyType", "requiredRoot", "rootMultiplicity", "minimalMultiplicity", "microStates" )
_moleculeData = ( ( "tyrosineDipeptide", "xyz"     , 0, 1, ( "Singles/Doubles", 10, 9, 0, 0, "Cardinal"        , 0, 1, 1, None       ) ), \
                  ( "allyl",             "radicals", 0, 2, ( "Singles"        ,  7, 8, 1, 0, "Fixed Fractional", 0, 2, 2, None       ) ), \
                  ( "allyl",             "radicals", 0, 2, ( "Doubles"        ,  7, 8, 1, 0, "Fixed Fractional", 0, 2, 2, None       ) ), \
                  ( "allyl",             "radicals", 0, 2, ( "Singles/Doubles",  7, 8, 1, 0, "Fixed Fractional", 0, 2, 2, None       ) ), \
                  ( "allyl",             "radicals", 0, 2, ( "User"           ,  7, 8, 1, 0, "Fixed Fractional", 0, 2, 2, _AllylSet1 ) ), \
                  ( "allyl",             "radicals", 0, 2, ( "User"           ,  3, 3, 1, 0, "Fixed Fractional", 0, 2, 2, _AllylSet2 ) ), \
                  ( "methylene",         "radicals", 0, 3, ( "Singles"        ,  6, 8, 1, 1, "Fixed Fractional", 0, 3, 1, None       ) ), \
                  ( "methylene",         "radicals", 0, 3, ( "Doubles"        ,  6, 8, 1, 1, "Fixed Fractional", 0, 3, 1, None       ) ), \
                  ( "methylene",         "radicals", 0, 3, ( "Full"           ,  6, 8, 1, 1, "Fixed Fractional", 0, 3, 1, None       ) ), \
                  ( "methylene",         "radicals", 0, 3, ( "Singles/Doubles",  6, 8, 1, 1, "Fixed Fractional", 0, 3, 1, None       ) ), \
                  ( "methylene",         "radicals", 0, 3, ( "Singles/Doubles",  6, 8, 1, 1, "Fixed Fractional", 0, 5, 1, None       ) ), \
                  ( "methylene",         "radicals", 0, 3, ( "Full"           ,  6, 8, 1, 1, "Fixed Fractional", 0, 5, 1, None       ) ), \
                  ( "methylene",         "radicals", 0, 3, ( "Full"           ,  6, 8, 1, 1, "Fixed Fractional", 0, 5, 5, None       ) ), \
                  ( "water",             "xyz"     , 0, 1, ( "Singles/Doubles",  8, 6, 0, 0, "Cardinal"        , 0, 1, 1, None       ) )  )

_AlgorithmKeywords = { "CIAlgorithm"    : "Full" ,
                       "checkAlgorithm" : False    ,
                       "numberOfStates" : 100      ,
                       "eigenvalueSolverIterations"      : 100000,
                       "eigenvalueSolverPreconditioning" :  False,
                       "eigenvalueSolverTolerance"       : 1.0e-12 }

#===================================================================================================================================
# . Parameters for test.
#===================================================================================================================================
# . Tolerances.
_EnergyTolerance                = 0.1
_GradientAbsoluteErrorTolerance = 1.0e-02

#===================================================================================================================================
# . Class for a CI test system.
#===================================================================================================================================
class CITestSystem ( object ):
    """CI test system."""

    def __init__ ( self, **kwargs ):
        """Constructor."""
        for ( attribute, value ) in kwargs.iteritems ( ):
            setattr ( self, attribute, value )

    def GetSystem ( self, log = logFile ):
        """Get the system with the energy model defined."""
        # . Define the QC model.
        kwargs                    = { key : value for ( key, value ) in zip ( _qcModelKeywordLabels, self.qcModelOptions ) }
        kwargs["converger"      ] = DIISSCFConverger ( densityTolerance = 1.0e-10, maximumSCFCycles = 250 )
        kwargs["keepOrbitalData"] = True
        kwargs.update ( _AlgorithmKeywords )
        qcModel = QCModelMNDO ( "am1", **kwargs )
        # . Read the molecule.
        molecule       = XYZFile_ToSystem  ( os.path.join ( _dataPath, self.directory, self.label + ".xyz" ) )
        molecule.label = self.label
        molecule.electronicState = ElectronicState ( charge = self.charge, multiplicity = self.multiplicity )
        molecule.DefineQCModel ( qcModel )
        molecule.Summary ( log = log )
        # . Finish up.
        return molecule

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class MNDOCIEnergiesTest ( TestCase ):
    """A test case for calculating MNDO CI energies."""

    def __init__ ( self, *args ):
        """Constructor."""
        super ( MNDOCIEnergiesTest, self ).__init__ ( *args )
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

        # . Loop over systems.
        for testSystem in self.ciTestSystems:

            # . Get the molecule.
            molecule = testSystem.GetSystem ( log = log )

            # . Determine an energy.
            try:
                energy  = molecule.Energy ( log = log, doGradients = True )

                # . Charges.
                charges = molecule.AtomicCharges ( )
                charges.Print ( log = log, title = "Charges" )
                if log is not None: log.Paragraph ( "Total Charge = {:.3f}".format ( sum ( charges ) ) )
                charges = molecule.AtomicCharges ( spinDensities = True )
                charges.Print ( log = log, title = "Spin Densities" )
                if log is not None: log.Paragraph ( "Total Spin Density = {:.3f}".format ( sum ( charges ) ) )

                # . Printing.
                if self.doPrinting:
                    molecule.configuration.qcState.CIWavefunctionSummary ( log = log )
                    molecule.configuration.qcState.CIVectorsTable ( log = log, nvectors = 20 )
                    molecule.energyModel.qcModel.PrintOnePDMArrays ( molecule.configuration, log = log )

                # . Gradient testing.
                if self.testGradients and ( len ( molecule.atoms ) < self.maximumAtoms ):
                    of = SystemGeometryObjectiveFunction.FromSystem ( molecule )
                    gradientDeviation = of.TestGradients ( delta = 5.0e-05, log = log, tolerance = 3.0e-04 )
                    maximumGradientDeviation = max ( maximumGradientDeviation, gradientDeviation )

            # . Error.
            except Exception as e:
                numberErrors += 1
                if log is not None: log.Text ( "\nError occurred> " +  e.args[0] + "\n" )

        # . Get the observed and reference data.
        observed      = {}
        referenceData = TestDataSet ( "MNDO CI Energies" )
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
        self.ciTestSystems = []
        for values in _moleculeData:
            kwargs = { key : value for ( key, value ) in zip ( _keywordLabels, values ) }
            self.ciTestSystems.append ( CITestSystem ( **kwargs ) )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":

    # . Run the test.
    test = MNDOCIEnergiesTest ( )
    test.run ( )
