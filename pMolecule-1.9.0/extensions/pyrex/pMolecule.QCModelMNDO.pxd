#-------------------------------------------------------------------------------
# . File      : pMolecule.QCModelMNDO.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions                  cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Coordinates3                  cimport CCoordinates3, Coordinates3
from pCore.Integer1DArray                cimport CInteger1DArray, Integer1DArray, Integer1DArray_Allocate, Integer1DArray_GetItem, Integer1DArray_Set, Integer1DArray_SetItem
from pCore.Integer2DArray                cimport CInteger2DArray, Integer2DArray_Allocate, Integer2DArray_GetItem, Integer2DArray_Set, Integer2DArray_SetItem
from pCore.Matrix33                      cimport Matrix33, CMatrix33
from pCore.Memory                        cimport Memory_Allocate_Array_Integer, Memory_Deallocate_Integer
from pCore.Real1DArray                   cimport CReal1DArray, Real1DArray
from pCore.Real2DArray                   cimport CReal2DArray, Real2DArray
from pCore.Selection                     cimport Selection, CSelection
from pCore.Status                        cimport Status
from pCore.SymmetricMatrix               cimport CSymmetricMatrix, SymmetricMatrix
from pCore.Vector3                       cimport Vector3, CVector3
from pMolecule.ChargeConstraintContainer cimport CChargeConstraintContainer, ChargeConstraintContainer
from pMolecule.QCAtomContainer           cimport QCAtomContainer, CQCAtomContainer
from pMolecule.QCChargeModelOptions      cimport QCChargeModel
from pMolecule.QCMMInteractionState      cimport QCMMInteractionState, CQCMMInteractionState
from pMolecule.QCModel                   cimport QCModel
# #ifdef MNDOCI
from pMolecule.QCModelMNDOState          cimport MNDOCIState, QCModelMNDOState, QCModelMNDOState_Finalize, QCModelMNDOState_Initialize, QCModelMNDOState_Setup, CQCModelMNDOState
# else /*MNDOCI*/
#from pMolecule.QCModelMNDOState          cimport QCModelMNDOState, QCModelMNDOState_Finalize, QCModelMNDOState_Initialize, QCModelMNDOState_Setup, CQCModelMNDOState
# #endif /*MNDOCI*/
from pMolecule.QCOnePDM                   cimport QCOnePDMOccupancyType
from pMolecule.QCParameters               cimport QCParameters, CQCParameter

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
# #ifdef MNDOCI
cdef extern from "JDEigenvalueSolver.h":

    ctypedef struct JDEigenvalueSolver:
        Boolean usePreconditioning
        Integer maximumMatrixVectorMultiplications
        Integer printLevel
        Real    errorTolerance

cdef extern from "MNDOCIModel.h":

    ctypedef enum MNDOCIAlgorithm:
        MNDOCIAlgorithm_Direct = 0,
        MNDOCIAlgorithm_Full   = 1,
        MNDOCIAlgorithm_Sparse = 2

    ctypedef enum MNDOCIMethod:
        MNDOCIMethod_Doubles        = 0,
        MNDOCIMethod_Full           = 1,
        MNDOCIMethod_Singles        = 2,
        MNDOCIMethod_SinglesDoubles = 3,
        MNDOCIMethod_UserSpecified  = 4

    ctypedef struct MNDOCIModel:
        Boolean            checkAlgorithm
        Boolean            doAllStates
        Boolean            identifyRootSpin
        Boolean            localizeOrbitals
        Integer            activeElectrons
        Integer            activeOrbitals
        Integer            localizeStart
        Integer            localizeStop
        Integer            minimalMultiplicity
        Integer            numberOfStates
        Integer            requiredRoot
        Integer            rootMultiplicity
        MNDOCIAlgorithm    algorithm
        MNDOCIMethod       method
        CInteger2DArray   *microstates
        JDEigenvalueSolver eigenvalueSolver

    cdef  MNDOCIModel *MNDOCIModel_Allocate   ( Status *status )
    cdef  MNDOCIModel *MNDOCIModel_Clone      ( MNDOCIModel  *self, Status *status )
    cdef  void         MNDOCIModel_Deallocate ( MNDOCIModel **self )
    cdef  MNDOCIState *MNDOCIModel_MakeState  ( MNDOCIModel  *self, Integer multiplicity, Integer nalpha, Integer nbeta, Integer norbitals, Status *status )
# #endif /*MNDOCI*/

cdef extern from "QCModelMNDO.h":

    ctypedef struct CQCModelMNDO "QCModelMNDO":
        Boolean                  isSpinRestricted
        Boolean                  keepOrbitalData
        Boolean                  linkAtomRatio
        Real                fermiBroadening
        Integer                   numberFractionalAlphaHOOs
        Integer                   numberFractionalAlphaLUOs
        Integer                   numberFractionalBetaHOOs
        Integer                   numberFractionalBetaLUOs
        Integer                   numberFractionalHOOs
        Integer                   numberFractionalLUOs
# #ifdef MNDOCI
        MNDOCIModel          *ciModel
# #endif /*MNDOCI*/
        QCChargeModel         qcChargeModel
        QCOnePDMOccupancyType occupancyType

    cdef CQCModelMNDO *QCModelMNDO_Allocate            ( )
    cdef void          QCModelMNDO_CCFockLowdin        ( CQCModelMNDO  *self, CQCModelMNDOState *qcState, CChargeConstraintContainer *chargeConstraints, CReal1DArray *lambdas, Status *status )
    cdef Real          QCModelMNDO_CCFunction          ( CQCModelMNDO               *self              ,
                                                         CQCModelMNDOState          *qcState           ,
                                                         CChargeConstraintContainer *chargeConstraints ,
                                                         CReal1DArray               *lambdas           ,
                                                         Boolean                     buildFock         ,
                                                         CReal1DArray               *gradients         ,
                                                         CSymmetricMatrix           *hessian           ,
                                                         Status                     *status            )
    cdef void          QCModelMNDO_CCHessian           ( CQCModelMNDO  *self, CQCModelMNDOState *qcState, CChargeConstraintContainer *chargeConstraints, CSymmetricMatrix *hessian, Status *status )
# #ifdef MNDOCI
    cdef void          QCModelMNDO_CIEnergy            ( CQCModelMNDO  *self, CQCModelMNDOState *state )
    cdef CReal1DArray *QCModelMNDO_CIStateCharacters   ( CQCModelMNDO  *self, CQCModelMNDOState *state, CMatrix33 *rotation, CInteger1DArray *mapping,
                                                                               CSelection *stateIndices, Boolean includeCoreOrbitals, Status *status )
# #endif /*MNDOCI*/
    cdef CQCModelMNDO *QCModelMNDO_Clone               ( CQCModelMNDO  *self )
    cdef void          QCModelMNDO_Deallocate          ( CQCModelMNDO **self )
    cdef CVector3     *QCModelMNDO_DipoleMoment        ( CQCModelMNDO  *self, CQCModelMNDOState *qcState, CCoordinates3 *coordinates3, CVector3 *center )
    cdef void          QCModelMNDO_Fock                ( CQCModelMNDO  *self, CQCModelMNDOState *qcState, Real *eelectronic )
    cdef void          QCModelMNDO_Gradients           ( CQCModelMNDO  *self, CQCModelMNDOState *qcState )
    cdef void          QCModelMNDO_GridPointDensities  ( CQCModelMNDO  *self, CQCModelMNDOState *qcState, CCoordinates3 *coordinates3, CCoordinates3 *points, CReal1DArray *data, Boolean spinDensities )
    cdef void          QCModelMNDO_GridPointOrbitals   ( CQCModelMNDO  *self, CQCModelMNDOState *qcState, CCoordinates3 *coordinates3, CCoordinates3 *points, CReal1DArray *data, Integer norbitals, Integer *orbitals, Boolean useDensityP )
    cdef void          QCModelMNDO_GridPointPotentials ( CQCModelMNDO  *self, CQCModelMNDOState *qcState, CCoordinates3 *coordinates3, CCoordinates3 *points, CReal1DArray *data )
    cdef void          QCModelMNDO_Integrals           ( CQCModelMNDO  *self, CQCModelMNDOState *qcState )
    cdef void          QCModelMNDO_LocalizeOrbitals    ( CQCModelMNDO  *self, CQCModelMNDOState *state, Boolean doQ, Integer startOrbital, Integer stopOrbital, Status *status )
    cdef void          QCModelMNDO_LowdinCharges       ( CQCModelMNDO  *self, CQCModelMNDOState *qcState, Boolean spinDensities, CReal1DArray *qcCharges )
    cdef void          QCModelMNDO_MayerBondOrders     ( CQCModelMNDO  *self, CQCModelMNDOState *qcState, CSymmetricMatrix *bondOrders, CReal1DArray *charges, CReal1DArray *freevalence, CReal1DArray *totalvalence )
    cdef CReal1DArray *QCModelMNDO_OrbitalCharacters   ( CQCModelMNDO  *self, CQCModelMNDOState *state, CMatrix33 *rotation, CInteger1DArray *mapping,
                                                                                     Boolean useDensityP, CSelection *orbitalIndices, Status *status )
    cdef CReal1DArray *QCModelMNDO_OrbitalEnergies     ( CQCModelMNDO  *self, CQCModelMNDOState *qcState, Boolean useDensityPORBITALS, Integer *homo, Integer *lumo )
    cdef CReal2DArray *QCModelMNDO_Orbitals            ( CQCModelMNDO  *self, CQCModelMNDOState *qcState, Boolean useDensityP )
    cdef void          QCModelMNDO_QCMMFockLowdin      ( CQCModelMNDO  *self, CQCModelMNDOState *qcState, Real *eelectronic )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class QCModelMNDO ( QCModel ):

    cdef CQCModelMNDO *cObject
    cdef public object     isOwner
    cdef public object     mndoModel
    cdef public object     orbitalBasis
