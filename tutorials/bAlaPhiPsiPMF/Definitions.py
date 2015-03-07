"""Some definitions."""

import glob, math, os.path

# . Output directory for large files.
scratchPath = os.getenv ( "PDYNAMO_SCRATCH" )
if scratchPath is None:
    outPath = ""
else:
    outPath = os.path.join ( scratchPath, "tutorials"      )
    if not os.path.exists  ( outPath ): os.mkdir ( outPath )
    outPath = os.path.join ( outPath    , "bAlaPhiPsiPMF"  )
    if not os.path.exists  ( outPath ): os.mkdir ( outPath )
    outPath = os.path.join ( outPath    , "generatedFiles" )
    if not os.path.exists  ( outPath ): os.mkdir ( outPath )

# . Data path.
dataPath = os.path.join ( os.getenv ( "PDYNAMO_ROOT" ), "tutorials", "bAlaPhiPsiPMF", "data" )

# . Phi/psi atom indices.
phiAtomIndices = ( 4, 6,  8, 14 )
psiAtomIndices = ( 6, 8, 14, 16 )

# . pDynamo imports.
from pBabel           import DCDTrajectoryFileReader                   , \
                             DCDTrajectoryFileWriter                   , \
                             MOLFile_ToSystem                          , \
                             SystemSoftConstraintTrajectory            , \
                             SystemSoftConstraintTrajectoryDataHandler , \
                             XYZFile_FromSystem                        , \
                             XYZFile_ToCoordinates3
from pCore            import logFile                                   , \
                             NormalDeviateGenerator                    , \
                             Pickle                                    , \
                             RandomNumberGenerator                     , \
                             RegularHistogram                          , \
                             RegularHistogramDimension                 , \
                             Statistics                                , \
                             Unpickle
from pMolecule        import MMModelOPLS                               , \
                             NBModelFull                               , \
                             SoftConstraintContainer                   , \
                             SoftConstraintDihedral                    , \
                             SoftConstraintEnergyModelHarmonic
from pMoleculeScripts import ConjugateGradientMinimize_SystemGeometry  , \
                             LangevinDynamics_SystemGeometry           , \
                             WHAM_Bootstrapping                        , \
                             WHAM_ConjugateGradientMinimize
