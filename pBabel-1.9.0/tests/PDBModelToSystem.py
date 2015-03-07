"""Read PDB model files and convert them to systems using atom data from original PDB files."""

import glob, os, os.path

from pBabel           import PDBModel_FromModelFile, PDBFile_ToPDBModel, XYZFile_FromSystem
from pCore            import RandomNumberGenerator, TestCase
from pMolecule        import MMModelOPLS, NBModelABFS
from pMoleculeScripts import BuildHydrogenCoordinates3FromConnectivity

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class PDBModelToSystemTest ( TestCase ):
    """A test case for reading and writing PDB model files."""

    def runTest ( self ):
        """The test."""

        # . Initialization.
        isOK         = True
        numberErrors = 0

        # . Energy models.
        mmModel = MMModelOPLS ( "protein" )
        nbModel = NBModelABFS ( )

        # . Paths.
        dataPath  = os.path.join ( os.getenv ( "PDYNAMO_PBABEL" ), "data" )
        modelPath = os.path.join ( dataPath, "pdbModel" )
        pdbPath   = os.path.join ( dataPath, "pdb"      )
        outPath   = None
        if self.resultPath is not None:
            outPath = os.path.join ( self.resultPath, "xyz" )
        log = self.GetLog ( )

        # . Set up the output directory.
        if outPath is not None:
            if not os.path.exists ( outPath ): os.mkdir ( outPath )
            outFiles = glob.glob ( os.path.join ( outPath, "*.xyz" ) )
            for outFile in outFiles: os.remove ( outFile )

        # . Get the files to process.
        modelFiles = glob.glob ( os.path.join ( modelPath, "*.model" ) )

        # . Get the model file names.
        pdbNames = set ( )
        for modelFile in modelFiles:
            ( head, tail ) = os.path.split ( modelFile )
            pdbNames.add ( tail[0:-6] )
        pdbNames = list ( pdbNames )
        pdbNames.sort ( )

        # . Loop over the files.
        for pdbName in pdbNames:

            # . Check file names.
            modelFile = os.path.join ( modelPath, pdbName + ".model" )
            pdbFile   = os.path.join ( pdbPath,   pdbName + ".pdb"   )
            if os.path.exists ( modelFile ) and os.path.exists ( pdbFile ):

                # . Get the model and its raw counterpart.
                model1   = PDBModel_FromModelFile ( modelFile, log = log )
                rawModel = PDBFile_ToPDBModel     ( pdbFile  , log = log )

                # . Route 1 - make an atomic model.
                model1.Summary ( log = log )
                try:

                    # . Make the atomic model.
                    model1.MakeAtomicModelFromComponentLibrary ( log = log )
                    model1.ExtractAtomData ( rawModel, log = log )
                    model1.Summary ( log = log )

                    # . Make a system.
                    system1 = model1.MakeSystem ( )
                    system1.Summary ( log = log )

                    # . Add energy models.
                    system1.DefineMMModel ( mmModel, log = log )
                    system1.DefineNBModel ( nbModel )

                    # . Build as many undefined coordinates as possible.
                    if system1.coordinates3.numberUndefined > 0:
                        rng = RandomNumberGenerator.WithSeed ( 117513 )
                        BuildHydrogenCoordinates3FromConnectivity ( system1, log = log, randomNumberGenerator = rng )

                    # . Calculate an energy if all coordinates have been defined.
                    if system1.coordinates3.numberUndefined <= 0:
                        system1.Energy ( log = log, doGradients = True )

                # . Error.
                except Exception as e:
                    numberErrors += 1
                    if log is not None: log.Text ( "\nError occurred> " +  e.args[0] + "\n" )

                # . Route 2 - extract atoms.
                model2  = PDBModel_FromModelFile ( modelFile, log = log )
                model2.ExtractAtoms ( rawModel, log = log )
                model2.Summary ( log = log )
                system2 = model2.MakeSystem ( )
                system2.Summary ( log = log )

                # . Output the xyz file if there are no undefined coordinates.
                n = system2.coordinates3.numberUndefined
                if n > 0:
                    if log is not None:
                        log.Paragraph ( "System has {:d} undefined coordinates.".format ( n ) )
                elif outPath is not None:
                    XYZFile_FromSystem ( os.path.join ( outPath, pdbName + ".xyz" ), system2 )

                # . Separator.
                if log is not None: log.Separator ( )

        # . Success/failure.
        self.assertTrue ( isOK and ( numberErrors == 0 ) )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":

    # . Run the test.
    test = PDBModelToSystemTest ( )
    test.run ( )

