"""Test for reading Amber Top and Crd Files."""

import glob, math, os

from pBabel           import AmberCrdFile_ToCoordinates3, AmberTopologyFile_ToSystem
from pCore            import logFile, LogFileActive, Pickle, TestCase, TestDataSet, TestReal, TextLogFileWriter, Unpickle
from pMolecule        import NBModelFull, SystemGeometryObjectiveFunction
from pMoleculeScripts import HardSphereIonMobilities

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Molecule names.
_SystemLabels = ( "dna", "glucose" )

# . Options.
_AbsoluteErrorTolerance         = 0.5
_GradientAbsoluteErrorTolerance = 1.0e-03
_MaximumAtoms                   = 100

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class AmberTopCrdReadTest ( TestCase ):
    """A test case for reading Amber Top and Crd files."""

    def GenerateReferenceData ( self ):
        self.generateReferenceData = True
        return True

    def runTest ( self ):
        """The test."""

        # . Initialization.
        if self.generateReferenceData: referenceData = TestDataSet ( "AMBER Test" )
        else:                          observed = {}

        # . Output setup.
        dataPath = os.path.join ( os.getenv ( "PDYNAMO_PBABEL" ), "data", "amber" )
        log = self.GetLog ( )

        # . Models.
        nbModel = NBModelFull ( )

        # . Loop over the molecules.
        for label in _SystemLabels:

            # . Read the data.
            molecule              = AmberTopologyFile_ToSystem  ( os.path.join ( dataPath, label + ".top" ), log = log )
            molecule.coordinates3 = AmberCrdFile_ToCoordinates3 ( os.path.join ( dataPath, label + ".crd" ), log = log )

            # . Calculation.
            molecule.DefineNBModel ( nbModel )
            molecule.Summary ( log = log )
            if log is not None: log.Text ( "\nFormula = " + molecule.atoms.FormulaString ( ) + ".\n" )
            energy = molecule.Energy ( doGradients = True, log = log )

            # . Get the dictionary of energies.
            localObserved = molecule.configuration.energyTerms.Terms ( asDictionary = True )
            localObserved["Potential Energy"] = energy

            # . Test the gradients.
            if len ( molecule.atoms ) <= _MaximumAtoms:
                of = SystemGeometryObjectiveFunction.FromSystem ( molecule )
                localObserved["Gradient Error"] = of.TestGradients ( log = log )

            # . Generate reference data.
            if self.generateReferenceData:
                localData = TestDataSet ( label, parent = referenceData )
                for ( key, value ) in localObserved.iteritems ( ):
                    if key == "Gradient Error": localData.AddDatum ( TestReal ( key, value, localData, absoluteErrorTolerance = _GradientAbsoluteErrorTolerance ) )
                    else:                       localData.AddDatum ( TestReal ( key, value, localData, absoluteErrorTolerance = _AbsoluteErrorTolerance, toleranceFormat = "{:.3f}", valueFormat = "{:.3f}" ) )
                referenceData.AddDatum ( localData )
            # . Accumulate observed data.
            else: observed[label] = localObserved

        # . Generate the reference data.
        if self.generateReferenceData:
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
    test = AmberTopCrdReadTest ( )
    test.run ( )
