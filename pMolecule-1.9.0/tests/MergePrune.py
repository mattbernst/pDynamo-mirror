"""Testing merging and pruning."""

import glob, os

from pBabel  import MOLFile_ToSystem
from pCore   import Clone, Selection, TestCase, TestDataSet, TestReal, Vector3
from pMolecule import MMModelOPLS, NBModelABFS

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Models.
_mmModel = MMModelOPLS ( "bookSmallExamples" )
_nbModel = NBModelABFS ( )

# . The molecules to test.
_Molecules = ( "bAla_c7eq", "chlorideAnion", "water" )

# . The tests.
# . The number of molecules of each type to merge followed by the indices of the molecules to prune.
_Tests = ( ( ( 0, 2, 0 ), ( 0, ) ),
           ( ( 0, 2, 0 ), ( 1, ) ),
           ( ( 0, 4, 0 ), ( 1, 2 ) ),
           ( ( 0, 0, 2 ), ( 0, ) ),
           ( ( 0, 0, 2 ), ( 1, ) ),
           ( ( 0, 0, 4 ), ( 1, 2 ) ),
           ( ( 2, 0, 0 ), ( 0, ) ),
           ( ( 2, 0, 0 ), ( 1, ) ),
           ( ( 4, 0, 0 ), ( 1, 2 ) ),
           ( ( 1, 1, 1 ), ( 0, ) ),
           ( ( 1, 1, 1 ), ( 1, ) ),
           ( ( 1, 1, 1 ), ( 2, ) ),
           ( ( 2, 2, 2 ), ( 1, 3, 5 ) ) )

# . Tolerances.
_EnergyAbsoluteErrorTolerance = 1.0e-04

# . Translation.
_Displacement = 25.0

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class MergePruneTest ( TestCase ):
    """Test merge and prune."""

    def runTest ( self ):
        """The test."""

        # . Initialization.
        dataPath    = os.path.join ( os.getenv ( "PDYNAMO_PMOLECULE" ), "data", "mol" )
        log         = self.GetLog ( )
        translation = Vector3.Uninitialized ( )

        # . Get the individual systems.
        energies  = []
        molecules = []
        for ( i, label ) in enumerate ( _Molecules ):
            molecule       = MOLFile_ToSystem ( os.path.join ( dataPath, label + ".mol" ) )
            molecule.label = label
            molecule.DefineMMModel ( _mmModel )
            molecule.DefineNBModel ( _nbModel )
            molecule.Summary ( log = log )
            molecule.coordinates3.TranslateToCenter ( )
            molecule.configuration.Clear ( )
            energies.append  ( molecule.Energy  ( log = log, doGradients = True ) )
            molecules.append ( molecule )

        # . Data initialization.
        observed      = {}
        referenceData = TestDataSet ( "Merge/Prune Energies" )

        # . Loop over the tests.
        for ( testIndex, ( moleculeFrequencies, moleculePruneIndices ) ) in enumerate ( _Tests ):

            # . Heading.
            if log is not None: log.Heading ( "Merge/Prune Test {:d}".format ( testIndex ), QBLANKLINE = True )

            # . Initialization.
            mergedEnergy = 0.0
            prunedEnergy = 0.0
            translation.Set ( 0.0 )

            # . Gather items.
            index        = 0
            numberAtoms  = 0
            pruneIndices = []
            toMerge      = []
            for ( i, frequency ) in enumerate ( moleculeFrequencies ):
                molecule     = molecules[i]
                mergedEnergy += energies[i] * frequency
                for f in range ( frequency ):
                    cloned = Clone ( molecule )
                    translation[0] += _Displacement
                    cloned.coordinates3.Translate ( translation )
                    toMerge.append ( cloned )
                    if index in moleculePruneIndices:
                        pruneIndices.extend ( range ( numberAtoms, numberAtoms + len ( cloned.atoms ) ) )
                        prunedEnergy += energies[i]
                    index       += 1
                    numberAtoms += len ( cloned.atoms )

            # . Merging.
            merged = toMerge[0].Merge ( toMerge[1:] )
            merged.Summary    ( log = log )
            eMerged = merged.Energy ( log = log )

            # . Pruning.
            pruned = merged.Prune ( Selection.FromIterable ( pruneIndices ) )
            pruned.Summary    ( log = log )
            ePruned = pruned.Energy ( log = log )

            # . Get the observed and reference data.
            for ( tag, eObserved, eReference ) in ( ( "Merged Energy {:d}".format ( testIndex ), eMerged, mergedEnergy ), ( "Pruned Energy {:d}".format ( testIndex ), ePruned, prunedEnergy ) ):
                observed[tag] = eObserved
                referenceData.AddDatum ( TestReal ( tag, eReference, referenceData, absoluteErrorTolerance = _EnergyAbsoluteErrorTolerance, toleranceFormat = "{:.3f}", valueFormat = "{:.3f}" ) )

        # . Finish up.
        if log is not None: log.Separator ( )

        # . Check for success/failure.
        if len ( observed ) > 0:
            results = referenceData.VerifyAgainst ( observed )
            results.Summary ( log = log, fullSummary = self.fullVerificationSummary )
            isOK = results.WasSuccessful ( )
        else:
            isOK = True

        # . Success/failure.
        self.assertTrue ( isOK )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":

    # . Run the test.
    test = MergePruneTest ( )
    test.run ( )
