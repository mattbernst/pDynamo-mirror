"""Test pickling and unpickling."""

import glob, math, os, os.path

from pBabel           import PDBFile_ToSystem
from pCore            import Pickle, TestCase, Unpickle, YAMLPickle, YAMLUnpickle
from pMolecule        import MMModelCHARMM, MMModelOPLS, NBModelABFS
from pMoleculeScripts import BuildHydrogenCoordinates3FromConnectivity

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
_Destination = "picklingMM"
_Tolerance   = 0.1

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class PickleMMTest ( TestCase ):
    """Test pickling and unpickling."""

    def runTest ( self ):
        """The test."""

        # . Paths.
        dataPath = os.path.join ( os.getenv ( "PDYNAMO_PMOLECULE" ), "data", "pdb" )
        if self.resultPath is None: outPath  = os.path.join ( os.getenv ( "PDYNAMO_SCRATCH" ), _Destination )
        else:                       outPath  = os.path.join ( self.resultPath, _Destination )
        if not os.path.exists ( outPath ): os.mkdir ( outPath )
        log = self.GetLog ( )

        # . Models.
        mmModel = MMModelCHARMM ( "c36a2" )
#        mmModel = MMModelOPLS ( "protein" )
        nbModel = NBModelABFS ( )

        # . Get all files.
        pdbFiles = glob.glob ( os.path.join ( dataPath, "*.pdb" ) )
        pdbFiles.sort ( )

        # . Read all PDB files.
        failures      = 0
        otherFailures = 0
        successes     = 0
        for pdbFile in pdbFiles:

            ( head, tail ) = os.path.split ( pdbFile )
            tag = tail[0:-4]

            if log is not None: log.Text ( "\nProcessing " + pdbFile + ":\n" )
            system = PDBFile_ToSystem ( pdbFile, log = log, useComponentLibrary = True )
            BuildHydrogenCoordinates3FromConnectivity ( system )
            try:

                # . Setup.
                if system.coordinates3.numberUndefined > 0: raise
                system.DefineMMModel ( mmModel, log = log )
                system.DefineNBModel ( nbModel )
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
            summary.Entry ( "Total Tests"     , "{:d}".format ( 2 * len ( pdbFiles ) ) )
            summary.Entry ( "Loop Failures"   , "{:d}".format ( otherFailures        ) )
            summary.Stop  ( )

        # . Success/failure.
        self.assertTrue ( ( failures == 0 ) and ( otherFailures == 0 ) )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":

    # . Run the test.
    test = PickleMMTest ( )
    test.run ( )
