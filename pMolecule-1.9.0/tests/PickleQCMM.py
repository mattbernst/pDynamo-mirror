"""Test pickling and unpickling."""

import glob, math, os, os.path

from pBabel    import PDBFile_ToSystem
from pCore     import Pickle, Selection, TestCase, Unpickle, YAMLPickle, YAMLUnpickle
from pMolecule import AtomSelection, CrystalClassCubic, ElectronicState, MMModelCHARMM, NBModelABFS, QCModelMNDO

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Scalars.
_BoxSize     = 28.0
_Destination = "picklingQCMM"
_Tolerance   = 0.1

# . Atom names.
_Tags = ( "CB",  "CG",  "CD1",  "CD2",  "CE1",  "CE2",  "CZ",  "OH",  "HB2",  "HB3",  "HD1",  "HD2",  "HE1",  "HE2",  "HH" )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class PickleQCMMTest ( TestCase ):
    """Test pickling and unpickling."""

    def runTest ( self ):
        """The test."""

        # . Initialization.
        failures      = 0
        otherFailures = 0
        successes     = 0

        # . Paths.
        dataPath = os.path.join ( os.getenv ( "PDYNAMO_PMOLECULE" ), "data", "pdb" )
        if self.resultPath is None: outPath  = os.path.join ( os.getenv ( "PDYNAMO_SCRATCH" ), _Destination )
        else:                       outPath  = os.path.join ( self.resultPath, _Destination )
        if not os.path.exists ( outPath ): os.mkdir ( outPath )
        log = self.GetLog ( )

        # . Models.
        mmModel = MMModelCHARMM ( "c36a2" )
        nbModel = NBModelABFS ( )
        qcModel = QCModelMNDO ( )

        # . Get the file.
        pdbFile = os.path.join ( dataPath, "2E4E_folded_solvated.pdb" )
        ( head, tail ) = os.path.split ( pdbFile )
        tag = tail[0:-4]

        if log is not None: log.Text ( "\nProcessing " + pdbFile + ":\n" )
        system = PDBFile_ToSystem ( pdbFile, log = log, useComponentLibrary = True )
        try:

            # . Fixed atoms.
            fixedAtoms = Selection ( ~ AtomSelection.FromAtomPattern ( system, "A:*:*" ) )

            # . QC selection.
            indices = set ( )
            for atomTag in _Tags:
                indices.add ( system.sequence.AtomIndex ( "A:TYR.2:" + atomTag ) )
            tyrosine = Selection ( indices )

            # . Setup.
            system.electronicState = ElectronicState ( charge = 0, multiplicity = 1 )
            system.DefineFixedAtoms ( fixedAtoms )
            system.DefineSymmetry   ( crystalClass = CrystalClassCubic ( ), a = _BoxSize )
            system.DefineMMModel    ( mmModel, log = log )
            system.DefineQCModel    ( qcModel, qcSelection = tyrosine )
            system.DefineNBModel    ( nbModel )
            system.Summary ( log = log )
            referenceEnergy = system.Energy  ( log = log, doGradients = True )

            # . Pickling.
            pklFile  = os.path.join ( outPath, tag + ".pkl"  )
            yamlFile = os.path.join ( outPath, tag + ".yaml" )
            Pickle     ( pklFile , system )
            YAMLPickle ( yamlFile, system )

            # . Unpickling.
            pklSystem = Unpickle ( pklFile )
            pklSystem.label += " (Pickled)"
            pklSystem.Summary ( log = log )
            pklEnergy = pklSystem.Energy  ( log = log, doGradients = True )
            if math.fabs ( referenceEnergy - pklEnergy  <= _Tolerance ): successes += 1
            else:                                                        failures  += 1

            yamlSystem = YAMLUnpickle ( yamlFile )
            yamlSystem.label += " (YAMLPickled)"
            yamlSystem.Summary ( log = log )
            yamlEnergy = yamlSystem.Energy  ( log = log, doGradients = True )
            if math.fabs ( referenceEnergy - yamlEnergy  <= _Tolerance ): successes += 1
            else:                                                         failures  += 1

        except Exception as e:
            otherFailures += 1
            if log is not None: log.Text ( "\nError occurred> " +  e.args[0] + "\n" )

        # . Summary of results.
        if log is not None:
            summary = log.GetSummary ( )
            summary.Start ( "Pickle Tests" )
            summary.Entry ( "Pickle Successes", "{:d}".format ( successes            ) )
            summary.Entry ( "Pickle Failures" , "{:d}".format ( failures             ) )
            summary.Entry ( "Total Tests"     , "2"                                    )
            summary.Entry ( "Loop Failures"   , "{:d}".format ( otherFailures        ) )
            summary.Stop  ( )

        # . Success/failure.
        self.assertTrue ( ( failures == 0 ) and ( otherFailures == 0 ) )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":

    # . Run the test.
    test = PickleQCMMTest ( )
    test.run ( )
