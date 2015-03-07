#-------------------------------------------------------------------------------
# . File      : pMolecule.QCModelDFT.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions             cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Coordinates3             cimport CCoordinates3, Coordinates3
from pCore.Integer1DArray           cimport Integer1DArray
from pCore.Memory                   cimport Memory_Allocate_Array_Integer, Memory_Deallocate_Integer
from pCore.Real1DArray              cimport CReal1DArray, Real1DArray
from pCore.Real2DArray              cimport CReal2DArray, Real2DArray
from pCore.Status                   cimport Status, Status_Continue
from pCore.SymmetricMatrix          cimport SymmetricMatrix, CSymmetricMatrix
from pCore.Vector3                  cimport CVector3, Vector3
from pMolecule.DFTFunctionalModel   cimport CDFTFunctionalModel, DFTFunctionalModel_Deallocate, DFTFunctionalModel_MakeFromIDs
from pMolecule.DFTGrid              cimport DFTGridAccuracy, DFTGridAccuracy_VeryLow, DFTGridAccuracy_Low, DFTGridAccuracy_Medium, DFTGridAccuracy_High, DFTGridAccuracy_VeryHigh, DFTGrid_Construct
from pMolecule.QCAtomContainer      cimport QCAtomContainer, CQCAtomContainer
from pMolecule.QCChargeModelOptions cimport QCChargeModel
from pMolecule.QCMMInteractionState cimport QCMMInteractionState, CQCMMInteractionState
from pMolecule.QCModel              cimport QCModel
from pMolecule.QCModelDFTState      cimport QCModelDFTState, QCModelDFTState_Finalize, QCModelDFTState_Initialize, QCModelDFTState_Setup, CQCModelDFTState
from pMolecule.QCOnePDM             cimport QCOnePDM, QCOnePDM_MakeWeighted, QCOnePDM_Reorthogonalize, QCOnePDMOccupancyType
from pMolecule.QCParameters         cimport QCParameters, CQCParameter

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "QCModelDFT.h":

    ctypedef struct CQCModelDFT "QCModelDFT":
        Boolean               inCore
        Boolean               isSpinRestricted
        Boolean               keepOrbitalData
        Boolean               linkAtomRatio
        Real                  fermiBroadening
        Integer               numberFractionalAlphaHOOs
        Integer               numberFractionalAlphaLUOs
        Integer               numberFractionalBetaHOOs
        Integer               numberFractionalBetaLUOs
        Integer               numberFractionalHOOs
        Integer               numberFractionalLUOs
        CDFTFunctionalModel  *functionalModel
        DFTGridAccuracy       accuracy
        QCChargeModel         qcChargeModel
        QCOnePDMOccupancyType occupancyType

    cdef CQCModelDFT  *QCModelDFT_Allocate                ( )
    cdef void          QCModelDFT_AtomicCharges           ( CQCModelDFT  *self, CQCModelDFTState *qcState, QCChargeModel *qcChargeModel, Boolean spinDensities, CReal1DArray *qcCharges )
    cdef CQCModelDFT  *QCModelDFT_Clone                   ( CQCModelDFT  *self )
    cdef void          QCModelDFT_Deallocate              ( CQCModelDFT **self )
    cdef CVector3     *QCModelDFT_DipoleMoment            ( CQCModelDFT  *self, CQCModelDFTState *qcState, CCoordinates3 *coordinates3, CVector3 *center )
    cdef void          QCModelDFT_Fock                    ( CQCModelDFT  *self, CQCModelDFTState *qcState, Real *eelectronic )
    cdef void          QCModelDFT_Gradients               ( CQCModelDFT  *self, CQCModelDFTState *qcState )
    cdef void          QCModelDFT_GridPointDensities      ( CQCModelDFT  *self, CQCModelDFTState *qcState, CCoordinates3 *coordinates3, CCoordinates3 *points, CReal1DArray *data, Boolean spinDensities )
    cdef void          QCModelDFT_GridPointOrbitals       ( CQCModelDFT  *self, CQCModelDFTState *qcState, CCoordinates3 *coordinates3, CCoordinates3 *points, CReal1DArray *data, Integer norbitals, Integer *orbitals, Boolean useDensityP )
    cdef void          QCModelDFT_GridPointPotentials     ( CQCModelDFT  *self, CQCModelDFTState *qcState, CCoordinates3 *coordinates3, CCoordinates3 *points, CReal1DArray *data )
    cdef void          QCModelDFT_Integrals               ( CQCModelDFT  *self, CQCModelDFTState *qcState )
    cdef void          QCModelDFT_MayerBondOrders         ( CQCModelDFT  *self, CQCModelDFTState *qcState, QCChargeModel *qcChargeModel, CSymmetricMatrix *bondOrders, CReal1DArray *charges, CReal1DArray *freevalence, CReal1DArray *totalvalence )
    cdef CReal1DArray *QCModelDFT_OrbitalEnergies         ( CQCModelDFT  *self, CQCModelDFTState *qcState, Boolean useDensityPORBITALS, Integer *homo, Integer *lumo )
    cdef CReal2DArray *QCModelDFT_Orbitals                ( CQCModelDFT  *self, CQCModelDFTState *qcState, Boolean useDensityP )
    cdef void          QCModelDFT_QCMMFock                ( CQCModelDFT  *self, CQCModelDFTState *qcState, Real *eelectronic )
    cdef void          QCModelDFT_QCMMMakeWeightedDensity ( CQCModelDFT  *self, CQCModelDFTState *qcState )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class QCModelDFT ( QCModel ):

    cdef CQCModelDFT  *cObject
    cdef public object isOwner
    cdef public object densityBasis
    cdef public object functional
    cdef public object orbitalBasis

