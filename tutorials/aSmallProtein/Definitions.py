"""Some definitions."""

import glob, math, os.path

# . Output directory for large files.
scratchPath = os.getenv ( "PDYNAMO_SCRATCH" )
if scratchPath is None:
    outPath = ""
else:
    outPath = os.path.join ( scratchPath, "tutorials"      )
    if not os.path.exists  ( outPath ): os.mkdir ( outPath )
    outPath = os.path.join ( outPath    , "aSmallProtein"  )
    if not os.path.exists  ( outPath ): os.mkdir ( outPath )
    outPath = os.path.join ( outPath    , "generatedFiles" )
    if not os.path.exists  ( outPath ): os.mkdir ( outPath )

# . Data path.
dataPath    = os.path.join ( os.getenv ( "PDYNAMO_ROOT" ), "tutorials", "aSmallProtein", "data" )

# . pBabel, pCore and pMolecule imports.
from pBabel           import AmberTrajectoryFileReader                 , \
                             AmberTrajectoryFileWriter                 , \
                             MOLFile_ToSystem                          , \
                             PDBFile_FromSystem                        , \
                             PDBFile_ToSystem                          , \
                             XYZFile_FromSystem                        , \
                             XYZFile_ToCoordinates3
from pCore            import Clone                                     , \
                             logFile                                   , \
                             NormalDeviateGenerator                    , \
                             Pickle                                    , \
                             RandomNumberGenerator                     , \
                             Selection                                 , \
                             Statistics                                , \
                             Unpickle
from pMolecule        import AtomSelection                             , \
                             CrystalClassCubic                         , \
                             MMModelOPLS                               , \
                             NBModelABFS                               , \
                             SoftConstraintContainer                   , \
                             SoftConstraintDistance                    , \
                             SoftConstraintEnergyModelHarmonic         , \
                             SoftConstraintMultipleTether              , \
                             SymmetryParameters
from pMoleculeScripts import AddCounterIons                            , \
                             BuildHydrogenCoordinates3FromConnectivity , \
                             BuildSolventBox                           , \
                             ConjugateGradientMinimize_SystemGeometry  , \
                             DetermineSolvationParameters              , \
                             LangevinDynamics_SystemGeometry           , \
                             SolvateSystemBySuperposition              , \
                             SystemDensity
