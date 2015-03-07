#-------------------------------------------------------------------------------
# . File      : __init__.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""pBabel handles the I/O of molecular data in various formats.

   The name is derived from the well-known series of programs
   which transform between various molecular data formats.
   However, the modules here have by no means the same range
   of functionality.
"""

# . The version number of the package.
PBABEL_VERSION = "1.9.0"

# . Import various modules.
from AmberCrdFileReader             import AmberCrdFileReader                           , \
                                           AmberCrdFile_ToCoordinates3
from AmberCrdFileWriter             import AmberCrdFileWriter                           , \
                                           AmberCrdFile_FromSystem
from AmberTopologyFileReader        import AmberTopologyFileReader                      , \
                                           AmberTopologyFile_ToSystem
from AmberTrajectoryFileReader      import AmberTrajectory_ToSystemGeometryTrajectory   , \
                                           AmberTrajectoryFileReader
from AmberTrajectoryFileWriter      import AmberTrajectory_FromSystemGeometryTrajectory , \
                                           AmberTrajectoryFileWriter
from CHARMMCRDFileReader            import CHARMMCRDFileReader                          , \
                                           CHARMMCRDFile_ToCoordinates3                 , \
                                           CHARMMCRDFile_ToSequence                     , \
                                           CHARMMCRDFile_ToSystem
from CHARMMParameterFileReader      import CHARMMParameterContainer                     , \
                                           CHARMMParameterFileReader                    , \
                                           CHARMMParameterFiles_ToParameters
from CHARMMPSFFileReader            import CHARMMPSFFileReader                          , \
                                           CHARMMPSFFile_ToSystem
from CHARMMPSFFileWriter            import CHARMMPSFFileWriter                          , \
                                           CHARMMPSFFile_FromSystem
from CHARMMTopologyFileReader       import CHARMMTopologyFileReader
from CIFFileReader                  import CIFFileReader                                , \
                                           CIFFile_ToSystem                             , \
                                           CIFFile_ToSystems
from DCDTrajectoryFileReader        import DCDTrajectory_ToSystemGeometryTrajectory     , \
                                           DCDTrajectoryFileReader
from DCDTrajectoryFileWriter        import DCDTrajectory_FromSystemGeometryTrajectory   , \
                                           DCDTrajectoryFileWriter
from EMSLFileReader                 import EMSLG94File_ToGaussianBases
from ExportImport                   import ExportFileFormats                            , \
                                           ExportOptions                                , \
                                           ExportSystem                                 , \
                                           ImportCoordinates3                           , \
                                           ImportFileFormats                            , \
                                           ImportObjects                                , \
                                           ImportOptions                                , \
                                           ImportSystem
from fDynamoCRDFileReader           import fDynamoCRDFileReader                         , \
                                           fDynamoCRDFile_ToCoordinates3                , \
                                           fDynamoCRDFile_ToSequence                    , \
                                           fDynamoCRDFile_ToSystem
from GaussianCubeFileReader         import GaussianCubeFile_ToCoordinates3              , \
                                           GaussianCubeFile_ToSystem                    , \
                                           GaussianCubeFileReader
from GaussianCubeFileWriter         import GaussianCubeFile_FromSystemDensity           , \
                                           GaussianCubeFile_FromSystemOrbitals          , \
                                           GaussianCubeFile_FromSystemPotential
from GromacsCrdFileReader           import GromacsCrdFileReader                         , \
                                           GromacsCrdFile_Process                       , \
                                           GromacsCrdFile_ToCoordinates3                , \
                                           GromacsCrdFile_ToSymmetry
from GromacsTopologyFileReader      import GromacsParameterContainer                    , \
                                           GromacsFileReader                            , \
                                           GromacsParameterReader                       , \
                                           GromacsDefinitionsReader                     , \
                                           GromacsDefinitions_ToSystem                  , \
                                           GromacsParameters_ToParameters
from JaguarInputFileReader          import JaguarInputFileReader                        , \
                                           JaguarInputFile_ToCoordinates3               , \
                                           JaguarInputFile_ToSystem
from JaguarInputFileWriter          import JaguarInputFileWriter                        , \
                                           JaguarInputFile_FromSystem
from JaguarOutputFileReader         import JaguarOutputFileReader                       , \
                                           JaguarOutputFile_ToCoordinates3              , \
                                           JaguarOutputFile_ToSystem
from mmCIFFileReader                import mmCIFFileReader                              , \
                                           mmCIFFile_ToSystem
from mmCIFFileWriter                import mmCIFFileWriter                              , \
                                           mmCIFFile_FromSystem
from MOL2FileReader                 import MOL2FileReader                               , \
                                           MOL2File_ToAtomNames                         , \
                                           MOL2File_ToCharges                           , \
                                           MOL2File_ToCoordinates3                      , \
                                           MOL2File_ToSystem
from MOL2FileWriter                 import MOL2FileWriter                               , \
                                           MOL2File_FromSystem
from MOLFileReader                  import MOLFileReader                                , \
                                           MOLFile_ToCoordinates3                       , \
                                           MOLFile_ToSystem                             , \
                                           SDFFile_ToSystems
from MOLFileWriter                  import MOLFileWriter                                , \
                                           MOLFile_FromSystem
from MopacInputFileReader           import MopacInputFileReader                         , \
                                           MopacInputFile_ToCoordinates3                , \
                                           MopacInputFile_ToSystem
from OOGLOffFileReader              import OOGLOffFileReader                            , \
                                           OOGLOffFile_ToPolygonalSurface
from OOGLOffFileWriter              import OOGLOffFileWriter                            , \
                                           OOGLOffFile_FromPolygonalSurface
from ORCAOutputFileReader           import ORCAOutputFileReader                         , \
                                           ORCAOutputFile_ToAbsorptionSpectrum          , \
                                           ORCAOutputFile_ToCoordinates3                , \
                                           ORCAOutputFile_ToSystem
from PDBComponent                   import PDBComponent                                 , \
                                           PDBComponentAtom                             , \
                                           PDBComponentBond                             , \
                                           PDBComponentLink                             , \
                                           PDBComponentVariant
from PDBComponentCIFFileReader      import PDBComponentCIFFile_ToComponents
from PDBComponentLibrary            import MakeDefaultPDBComponentLibrary               , \
                                           PDBComponentLibrary
from PDBFileReader                  import PDBFileReader                                , \
                                           PDBFile_ToCoordinates3                       , \
                                           PDBFile_ToPDBModel                           , \
                                           PDBFile_ToSystem                             , \
                                           PQRFile_ToCoordinates3                       , \
                                           PQRFile_ToSystem
from PDBFileWriter                  import PDBFileWriter                                , \
                                           PDBFile_FromSystem
from PDBModel                       import PDBModel_FromModelFile                       , \
                                           PDBModel_ToModelFile                         , \
                                           PDBModelLink                                 , \
                                           PDBModelVariant
from SMILESReader                   import SMILESReaderError                            , \
                                           SMILES_ToSystem
from SMILESUtilities                import SMILESConnectivity                           , \
                                           SMILESConnectivityError
from SMILESWriter                   import SMILES_FromSystem
from SystemGeometryTrajectory       import SystemGeometryTrajectory
from SystemSoftConstraintTrajectory import SystemSoftConstraintTrajectory               , \
                                           SystemSoftConstraintTrajectoryDataHandler
from XYZFileReader                  import XYZFileReader                                , \
                                           XYZFile_ToCoordinates3                       , \
                                           XYZFile_ToSystem                             , \
                                           XYZFiles_ToSystemGeometryTrajectory
from XYZFileWriter                  import XYZFileWriter                                , \
                                           XYZFile_FromSystem                           , \
                                           XYZFiles_FromSystemGeometryTrajectory
from XYZTrajectoryFileReader        import XYZTrajectory_ToSystemGeometryTrajectory     , \
                                           XYZTrajectoryFileReader
