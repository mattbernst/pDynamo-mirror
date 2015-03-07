"""Test for reading CHARMM Param and PSF files."""

import glob, math, os

from pBabel           import CHARMMCRDFile_ToCoordinates3, CHARMMParameterFiles_ToParameters, CHARMMPSFFileReader, CHARMMPSFFile_ToSystem, PDBFile_FromSystem
from pCore            import logFile, LogFileActive, Pickle, TestCase, TestDataSet, TestReal, TextLogFileWriter, UNITS_ENERGY_KILOCALORIES_PER_MOLE_TO_KILOJOULES_PER_MOLE, Unpickle
from pMolecule        import AtomSelection, NBModelFull, SystemGeometryObjectiveFunction
from pMoleculeScripts import HardSphereIonMobilities

# . The energy values should be the same as CHARMM values to within 1.0e-3 kcal/mole except for non-bonding where the deviations
# . for the larger systems may reach 0.1 or so.

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Molecule names.
_SystemLabels = ( "ava", "citrate_synthase" )

# . Options.
_AbsoluteErrorTolerance         = 0.5
_GradientAbsoluteErrorTolerance = 1.0e-03
_MaximumAtoms                   = 100

# . Parameter sets.
_ParameterPaths = ( "par_all27_prot_na", "par_coa", "par_oaa" )

# . Patterns.
_Patterns = { "ava"              : ( "*:*:C", "*:VAL.2:*", "AAAA:*:*" ), \
              "citrate_synthase" : ( "*:*:C", "*:PRO.*:*", "", "AABW:*:*" ) }

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class CharmmPSFReadTest ( TestCase ):
    """A test case for reading CHARMM Param and PSF files."""

    def GenerateReferenceData ( self ):
        self.generateReferenceData = True
        return True

    def runTest ( self ):
        """The test."""

        # . Initialization.
        if self.generateReferenceData: referenceData = TestDataSet ( "Charmm Test" )
        else:                          observed = {}

        # . Output setup.
        dataPath = os.path.join ( os.getenv ( "PDYNAMO_PBABEL" ), "data", "charmm" )
        outPath  = None
        if self.resultPath is not None:
            outPath = os.path.join ( self.resultPath, "pdb" )
            if not os.path.exists ( outPath ): os.mkdir ( outPath )
        log = self.GetLog ( )

        # . Get the parameters.
        parameterPaths = []
        for parameterPath in _ParameterPaths:
            parameterPaths.append ( os.path.join ( dataPath, parameterPath + ".prm" ) )
        parameters = CHARMMParameterFiles_ToParameters ( parameterPaths, log = log )

        # . Generate systems.
        for label in _SystemLabels:

            if log is not None:
                log.Text ( "\n" + ( 80 * "=" ) + "\n" )
                log.Text ( label + "\n" )
                log.Text ( 80 * "=" + "\n" )
            system              = CHARMMPSFFile_ToSystem ( os.path.join ( dataPath, label + ".psfx" ), isXPLOR = True, parameters = parameters, log = log )
            system.coordinates3 = CHARMMCRDFile_ToCoordinates3 ( os.path.join ( dataPath, label + ".chm" ), log = log )
            system.label        = label
            system.DefineNBModel ( NBModelFull ( ) )
            system.Summary ( log = log )
            energy = system.Energy  ( log = log )
            log.Text ( "\nEnergy (kcal/mole) = {:.4f}\n".format ( energy / UNITS_ENERGY_KILOCALORIES_PER_MOLE_TO_KILOJOULES_PER_MOLE ) )

            # . Get the dictionary of energies.
            localObserved = system.configuration.energyTerms.Terms ( asDictionary = True )
            localObserved["Potential Energy"] = energy

            # . Test gradients.
            if len ( system.atoms ) <= _MaximumAtoms:
                of = SystemGeometryObjectiveFunction.FromSystem ( system )
                of.TestGradients ( log = log )
                localObserved["Gradient Error"] = of.TestGradients ( log = log )

            # . Write PDB file and do various sequence tests.
            if not self.generateReferenceData:
                if outPath is not None:
                    PDBFile_FromSystem ( os.path.join ( outPath, label + ".pdb" ), system, useSegmentEntityLabels = True )
                system.sequence.PrintComponentSequence ( log = log )
                log.Text ( "\nSelections:\n\n" )
                for pattern in _Patterns[label]:
                    log.Text ( "{:<30s} {:5d}\n".format ( pattern, len ( AtomSelection.FromAtomPattern ( system, pattern ) ) ) )

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
    test = CharmmPSFReadTest ( )
    test.run ( )
