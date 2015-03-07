#-------------------------------------------------------------------------------
# . File      : __init__.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""pCore contains foundation classes and functions for general Python/C programs."""

# . Add the root directory to the python path list.
import os, sys
sys.path.append ( os.path.split ( os.getcwd ( ) )[0] )
del os, sys

# . The version number of the package.
PCORE_VERSION = "1.9.0"

# . Import various modules.
# . Python.
from BakerOptimizer                   import BakerOptimizer
from Clone                            import Clone, DeepClone, ShallowClone
from CoreObjects                      import CLibraryError, pObjectProperty, SingleObjectContainer
from ConjugateGradientMinimizer       import ConjugateGradientMinimizer
from FIREMinimizer                    import FIREMinimizer
from Histogram                        import RegularHistogram, RegularHistogramBinMidPointIterator, RegularHistogramDimension
from LangevinVelocityVerletIntegrator import LangevinVelocityVerletIntegrator
from LBFGSMinimizer                   import LBFGSMinimizer
from LeapFrogIntegrator               import LeapFrogIntegrator
from LogFileWriter                    import logFile, LogFileActive, TextLogFileWriter, XHTMLLogFileWriter, PrintPriority_None, \
                                             PrintPriority_Verylow, PrintPriority_Low, PrintPriority_Medium,                    \
                                             PrintPriority_High, PrintPriority_Veryhigh, PrintPriority_Debug
from MoreThuenteLineSearch            import MoreThuenteLineSearcher
from MultiCubicSpline                 import MultiCubicSpline
from MultiDimensionalMinimizer        import MultiDimensionalMinimizer, MultiDimensionalMinimizerState
from NumericalIntegration             import AdaptiveSimpsonsRule, TrapezoidalRule
from ObjectiveFunction                import ObjectiveFunction, UniDimensionalObjectiveFunction
from QuasiNewtonMinimizer             import QuasiNewtonMinimizer
from RingFinder                       import FindRings, FindRingSets
from Serialization                    import Pickle, PickleFileExtension, RawObjectConstructor, Unpickle, \
                                             YAMLMappingFile_FromObject, YAMLMappingFile_ToObject,        \
                                             YAMLPickle, YAMLPickleFileExtension, YAMLUnpickle
from Singleton                        import Singleton
from SpectraHandling                  import BlackBodySpectrum   , \
                                             ContinuousSpectrum  , \
                                             DataRange           , \
                                             DataSet             , \
                                             GaussianLineShape   , \
                                             LorentzianLineShape , \
                                             RayleighSpectrum    , \
                                             StickSpectrum
from Statistics                       import Statistics, StatisticsAccumulator
from SteepestDescentMinimizer         import SteepestDescentMinimizer
from SteepestDescentPathFinder        import SteepestDescentPathFinder
from Testing                          import TestCase, TestDataResult, TestDataSet, TestReal, TestResult
from TextFile                         import TextFile, TextFileReader, TextFileReaderError, TextFileWriter, TextFileWriterError
from Time                             import CPUTime
from TimeAnalysis                     import Timings, TimingsAverager
from Tree                             import TreeBranchNode, TreeLeafNode, TreeRootNode
from Units                            import *
from VelocityVerletIntegrator         import VelocityVerletIntegrator

# . Extension types.
from Coordinates3                     import Coordinates3, Coordinates3_FromGrid
from Correlation                      import Correlation_AutoCorrelation, Correlation_CrossCorrelation, Correlation_DotProductAutoCorrelation, Correlation_DotProductCrossCorrelation
from CubicSpline                      import CubicSpline
from Integer1DArray                   import Integer1DArray
from Integer2DArray                   import Integer2DArray
from Matrix33                         import Matrix33, Matrix33_RotationAboutAxis
from PairList                         import CrossPairList, CrossPairListIterator, CrossPairList_FromIntegerPairs, SelfPairList, SelfPairListIterator, SelfPairList_FromIntegerPairs
from PairListGenerator                import CrossPairList_FromDoubleCoordinates3, CrossPairList_FromSingleCoordinates3, SelfPairList_FromCoordinates3, PairListGenerator
from PolygonalSurface                 import Isosurface_FromRegularGridData, PolygonalSurface
from RandomDeviateGenerator           import NormalDeviateGenerator
from RandomNumberGenerator            import RandomNumberGenerator
from Real1DArray                      import Real1DArray
from Real2DArray                      import Real2DArray
from RealNDArray                      import RealNDArray
from RegularGrid                      import RegularGrid, RegularGrid_FromDimensionData
from Selection                        import Selection, Selection_SortComparison
from SelectionContainer               import SelectionContainer
from SymmetricMatrix                  import SymmetricMatrix
from Transformation3                  import Transformation3, Transformation3_FromSymmetryOperationString
from Transformation3Container         import Transformation3Container, Transformation3Container_Identity
from Vector3                          import Vector3
