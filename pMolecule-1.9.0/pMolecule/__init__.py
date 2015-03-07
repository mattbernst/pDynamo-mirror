#-------------------------------------------------------------------------------
# . File      : __init__.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""pMolecule is a Python version of the Dynamo modeling program."""

# . Add the root directory to the python path list.
import os, sys
sys.path.append ( os.path.split ( os.getcwd ( ) )[0] )
del os, sys

# . The version number of the package.
PDYNAMO_VERSION = "1.9.0"

# . Python.
from Angle                                     import AngleContainer
from Atom                                      import Atom, AtomContainer
from AtomSelection                             import AtomSelection, AtomSelectionError
from AtomSelector                              import SQLAtomSelector
from Bond                                      import Bond, BondContainer, BondDefinition, \
                                                      AromaticDoubleBond, AromaticSingleBond, DoubleBond, NullBond, SingleBond, TripleBond, UndefinedBond
from Configuration                             import Configuration
from Connectivity                              import Connectivity
from ConnectivityPattern                       import ConnectivityPattern
from CrystalClass                              import CrystalClass             , \
                                                      CrystalClassCubic        , \
                                                      CrystalClassHexagonal    , \
                                                      CrystalClassMonoclinic   , \
                                                      CrystalClassOrthorhombic , \
                                                      CrystalClassRhombohedral , \
                                                      CrystalClassTetragonal   , \
                                                      CrystalClassTriclinic    , \
                                                      CrystalClass_FromSpaceGroupNumber
from Dihedral                                  import DihedralContainer
from ElectronicState                           import ElectronicState
from Element                                   import PeriodicTable
from EnergyTerms                               import EnergyTerms
from EnergyModel                               import EnergyModel
from HardConstraintContainer                   import HardConstraintContainer
from MMModel                                   import MMModel, MMModelAMBER, MMModelCHARMM, MMModelOPLS
from MMModelError                              import MMModelError
from MMSequence                                import MMSequenceAtom, MMSequenceComponent, MMSequenceLink, MMSequenceVariant
from MMTermContainer                           import MMTermContainer
from MultiLayerSystemGeometryObjectiveFunction import MultiLayerSystemGeometryObjectiveFunction
from SEAMObjectiveFunction                     import SEAMObjectiveFunction
from Sequence                                  import Sequence, SequenceComponent, SequenceEntity, SequenceLinearPolymer, SequenceLink, SequenceVariant
from SoftConstraint                            import SoftConstraintDihedral, SoftConstraintDistance, SoftConstraintAngleDotProduct, SoftConstraintMultipleDistance, SoftConstraintMultipleTether, SoftConstraintTether, \
                                                      SoftConstraintEnergyModel, SoftConstraintEnergyModelHarmonic, SoftConstraintEnergyModelHarmonicRange
from SoftConstraintContainer                   import SoftConstraintContainer
from Symmetry                                  import Symmetry
from System                                    import System
from SystemWithTimings                         import SystemWithTimings
from SystemGeometryObjectiveFunction           import SystemGeometryObjectiveFunction

# . Extension types.
from ADIISSCFConverger                         import ADIISSCFConverger
from ChargeConstraintContainer                 import ChargeConstraintContainer
from DIISSCFConverger                          import DIISSCFConverger
from CMAPDihedralContainer                     import CMAPDihedralContainer
from FourierDihedralContainer                  import FourierDihedralContainer
from GaussianBasis                             import GaussianBasis
from HarmonicAngleContainer                    import HarmonicAngleContainer
from HarmonicBondContainer                     import HarmonicBondContainer
from HarmonicImproperContainer                 import HarmonicImproperContainer
from LebedevLaikovGrid                         import LebedevLaikovGrid_GetGridPoints
from LJParameterContainer                      import LJParameterContainer, LJParameterContainer_FromEpsilonSigma, LJParameterContainer_FromTableCoefficients
from MMAtomContainer                           import MMAtomContainer
from MMTerm                                    import MMTerm
from MNDOParameters                            import MNDOParameters
from MonteCarloSystemGeometry                  import MonteCarlo_IsolateInteractionEnergy, MonteCarlo_ScaleIsolateInteractionParameters, MonteCarlo_SystemGeometry
from NBModelABFS                               import NBModelABFS
from NBModelFull                               import NBModelFull
from NBModelMonteCarlo                         import NBModelMonteCarlo
from NBModelMonteCarloState                    import NBModelMonteCarloState
from NBModelORCA                               import NBModelORCA
from NBModelSSBP                               import NBModelSSBP
from PairwiseInteraction                       import PairwiseInteractionABFS
from QCAtomContainer                           import QCAtomContainer, QCAtomContainer_FromAtomContainer
from QCMMLinkAtomCouplingOptions               import QCMMLinkAtomCoupling_ToEnum, QCMMLinkAtomCoupling_ToString
from QCMMInteractionState                      import QCMMInteractionState
from QCModel                                   import QCModel
from QCModelDFT                                import QCModelDFT
from QCModelMNDO                               import QCModelMNDO
from QCModelORCA                               import QCModelORCA, QCModelORCAState
from QCOnePDM                                  import QCOnePDM
from QCParameters                              import QCParameters, QCParameters_Define
from SymmetryParameters                        import SymmetryParameters
from SymmetryParameterGradients                import SymmetryParameterGradients
