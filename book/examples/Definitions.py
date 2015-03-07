"""Some definitions needed by the testsuite examples."""

import glob, math, os.path

# . Get the root directory.
rootDirectory = os.getenv ( "PDYNAMO_ROOT" )
if rootDirectory is None:
    # . Go two directories down from the current directory.
    rootDirectory = os.path.split ( os.path.split ( os.getcwd ( ) )[0] )[0]

# . The pdynamo data path.
dataPath = os.path.join ( rootDirectory, "book", "data" )

# . Paths for finding data files of various types.
molPath = os.path.join ( dataPath, "mol" )
pdbPath = os.path.join ( dataPath, "pdb" )
pklPath = os.path.join ( dataPath, "pkl" )
xyzPath = os.path.join ( dataPath, "xyz" )

# . Paths for generated files and logs.
# . The scratch directory.
pDynamoScratch = os.getenv ( "PDYNAMO_SCRATCH" )
if pDynamoScratch is None: pDynamoScratch = os.path.join ( rootDirectory, "scratch" )
if not os.path.exists ( pDynamoScratch ): os.mkdir ( pDynamoScratch )

# . The directory for files and logs.
bookExamples = os.path.join ( pDynamoScratch, "bookExamples" )
if not os.path.exists ( bookExamples ): os.mkdir ( bookExamples )

# . The directory for generated files.
scratchPath = os.path.join ( bookExamples, "generatedFiles" )
if not os.path.exists ( scratchPath ): os.mkdir ( scratchPath )

# . The directory for errors and logs.
errorPath = os.path.join ( bookExamples, "errors" )
logPath   = os.path.join ( bookExamples, "logs"   )

# . Output options.
errorExtension = ".err"
logExtension   = ".log"
xhtmlExtension = ".html"

# . pDynamo package imports.
from pBabel           import MOLFile_ToCoordinates3                       , \
                             MOLFile_ToSystem                             , \
                             PDBFile_ToSystem                             , \
                             SMILES_FromSystem                            , \
                             SMILES_ToSystem                              , \
                             SystemGeometryTrajectory                     , \
                             SystemSoftConstraintTrajectory               , \
                             XYZFile_FromSystem                           , \
                             XYZFile_ToCoordinates3                       , \
                             XYZFile_ToSystem
from pCore            import Clone                                        , \
                             ConjugateGradientMinimizer                   , \
                             CONSTANT_MOLAR_GAS                           , \
                             logFile                                      , \
                             NormalDeviateGenerator                       , \
                             Pickle                                       , \
                             RandomNumberGenerator                        , \
                             Selection                                    , \
                             Statistics                                   , \
                             UNITS_MASS_AMU_TO_KG                         , \
                             Unpickle                                     , \
                             Vector3                                      , \
                             XHTMLLogFileWriter
from pMolecule        import CrystalClassCubic                            , \
                             DIISSCFConverger                             , \
                             MMModelOPLS                                  , \
                             MonteCarlo_IsolateInteractionEnergy          , \
                             MonteCarlo_ScaleIsolateInteractionParameters , \
                             MonteCarlo_SystemGeometry                    , \
                             NBModelABFS                                  , \
                             NBModelFull                                  , \
                             NBModelMonteCarlo                            , \
                             NBModelORCA                                  , \
                             QCModelDFT                                   , \
                             QCModelMNDO                                  , \
                             QCModelORCA                                  , \
                             SoftConstraintContainer                      , \
                             SoftConstraintDistance                       , \
                             SoftConstraintEnergyModelHarmonic            , \
                             SoftConstraintEnergyModelHarmonicRange       , \
                             SoftConstraintTether                         , \
                             SymmetryParameters                           , \
                             System                                       , \
                             SystemGeometryObjectiveFunction
from pMoleculeScripts import BakerSaddleOptimize_SystemGeometry           , \
                             BuildCubicSolventBox                         , \
                             BuildHydrogenCoordinates3FromConnectivity    , \
                             BuildSolventBox                              , \
                             ChainOfStatesOptimizePath_SystemGeometry     , \
                             ConjugateGradientMinimize_SystemGeometry     , \
                             GrowingStringInitialPath                     , \
                             IdentifyUndefinedCoordinates3                , \
                             LeapFrogDynamics_SystemGeometry              , \
                             MergeByAtom                                  , \
                             NormalModes_SystemGeometry                   , \
                             NormalModesTrajectory_SystemGeometry         , \
                             PruneByAtom                                  , \
                             RadialDistributionFunction                   , \
                             SelfDiffusionFunction                        , \
                             SolventCubicBoxDimensions                    , \
                             SolvateSystemBySuperposition                 , \
                             SteepestDescentPath_SystemGeometry           , \
                             SystemDensity                                , \
                             ThermodynamicsRRHO_SystemGeometry            , \
                             VelocityVerletDynamics_SystemGeometry        , \
                             WHAM_ConjugateGradientMinimize
