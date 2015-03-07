#-------------------------------------------------------------------------------
# . File      : pMolecule.QCModelMNDOState.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.BlockStorage             cimport CBlockStorage
from pCore.cDefinitions             cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Coordinates3             cimport CCoordinates3, Coordinates3
from pCore.Integer1DArray           cimport CInteger1DArray, Integer1DArray_GetItem
from pCore.Real1DArray              cimport CReal1DArray, Real1DArray, Real1DArray_GetItem
from pCore.Real2DArray              cimport CReal2DArray, Real2DArray, Real2DArray_GetItem
from pCore.Status                   cimport Status, Status_Continue
from pCore.SymmetricMatrix          cimport SymmetricMatrix, CSymmetricMatrix
from pMolecule.QCAtomContainer      cimport QCAtomContainer, QCAtomContainer_Size, CQCAtomContainer
from pMolecule.QCMMInteractionState cimport QCMMInteractionState, CQCMMInteractionState
from pMolecule.QCOnePDM             cimport CQCOnePDM, QCOnePDM, QCOnePDM_SpinExpectationValues, QCOnePDMOccupancyType_FractionalVariable, QCOnePDMOccupancyType
from pMolecule.QCParameters         cimport QCParameters, CQCParameter

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
# #ifdef MNDOCI
cdef extern from "JDEigenvalueSolver.h":

    ctypedef struct JDEigenvalueSolverReport:
        Boolean isConverged
        Boolean solutionChecked
        Integer convergedPairs
        Integer numberMatrixVectorMultiplications
        Integer returnCode
        Real    eigenvalueError
        Real    eigenvectorError
        Real    normalizationError

cdef extern from "MNDOCIState.h":

    ctypedef struct MNDOCIConfiguration:
        Integer          nalphas
        Integer          nspqr
        Real             spin
        CInteger1DArray *alphas
        CInteger1DArray *betas
        CInteger1DArray *spqr

    ctypedef struct MNDOCIState:
        Boolean                 doFull
        Boolean                 doSparse
        Boolean                 fractionallyOccupiedInactive
        Boolean                 orbitalDegeneracies
        Boolean                 rootNotFound
        Boolean                 usePreconditioning
        Integer                 ciMatrixNonZero
        Integer                 localizeStart
        Integer                 localizeStop
        Integer                 nactive
        Integer                 nconfigurations
        Integer                 ncore
        Integer                 nelectrons
        Integer                 norbitals
        Integer                 numberOfStates
        Integer                 root
        Real                    baseline
        Real                    ciEnergy
        Real                    ciMatrixSparsity
        Real                    ecore
        Real                    requiredSpin
        Real                    rootenergy
        MNDOCIConfiguration    *configurations
#        CBlockStorage    *twoelectronintegrals
#        CSymmetricMatrix *oneelectronmatrix
#        DoubleSymmetricMatrix *moteis
#        DoubleSymmetricMatrix *twopdm
        CReal1DArray           *ciEnergies
        CReal1DArray           *ciVector
#        Real1DArray           *energies
#        Real1DArray           *ocore
        CReal1DArray           *spins
        CReal2DArray           *ciVectors
#        Real2DArray           *kpa
#        Real2DArray           *kpaMO
#        Real2DArray           *motei34
#        Real2DArray           *orbitals
#        RealNDArray           *motei234
#        CSymmetricMatrix  *ciMatrixFull
#        CSymmetricMatrix  *fcore
#        CSymmetricMatrix  *fcoreMO
#        CSymmetricMatrix  *onepdma
#        CSymmetricMatrix  *onepdmb
#        CSymmetricMatrix  *onepdmMOa
#        CSymmetricMatrix  *onepdmMOb
#        CSymmetricMatrix  *pcore
#        Real2DArray            activemos
        JDEigenvalueSolverReport eigenvalueSolverReport
        # . CI gradient.
        Boolean                  ciGradientError
        Boolean                  doGradients
        Integer                  cphfErrorFlag
        Integer                  cphfIterations
        Integer                  numberDegenerateRedundant
        Integer                  numberNonRedundant
        Integer                  numberRedundant
# #endif /*MNDOCI*/

cdef extern from "QCModelMNDOState.h":

    ctypedef struct CQCModelMNDOState "QCModelMNDOState":
        Boolean                   isConverged
        Real                 eelectronic
        Real                 enuclear
        Real                 eocc
        Real                 eoei
        Real                 eqcmm
        Real                 eqcqc
        Real                 etei
        Integer                    ncycles
        CQCAtomContainer      *qcAtoms
        CQCParameter          *qcParameters
        CCoordinates3         *coordinates3
        CCoordinates3         *gradients3
        CQCMMInteractionState *qcmmstate
        CQCOnePDM             *densityp
        CQCOnePDM             *densityq
        CCoordinates3         *qccoordinates3
        CSymmetricMatrix      *oneelectronmatrix
        CBlockStorage         *twoelectronintegrals
# #ifdef MNDOCI
        MNDOCIState           *cistate
# #endif /*MNDOCI*/

    cdef CQCModelMNDOState *QCModelMNDOState_Allocate       ( )                                          
    cdef void               QCModelMNDOState_BackupFock     ( CQCModelMNDOState  *self, Status *status ) 
    cdef void               QCModelMNDOState_Deallocate     ( CQCModelMNDOState **self )                 
    cdef Real               QCModelMNDOState_EnergyFromFock ( CQCModelMNDOState *self )
    cdef void               QCModelMNDOState_Finalize       ( CQCModelMNDOState  *self, Boolean QKEEPDATA )
    cdef void               QCModelMNDOState_Initialize     ( CQCModelMNDOState  *self, CCoordinates3 *coordinates3, CQCMMInteractionState *qcmmstate, CCoordinates3 *gradients3 )
    cdef Real               QCModelMNDOState_MakeDensities  ( CQCModelMNDOState  *self )
    cdef void               QCModelMNDOState_RestoreFock    ( CQCModelMNDOState  *self )
    cdef CQCModelMNDOState *QCModelMNDOState_Setup          ( CQCAtomContainer     *qcAtoms                   ,
                                                              CQCParameter         *qcParameters              ,
                                                              Real                  alphaCharge               ,
                                                              Real                  betaCharge                ,
                                                              QCOnePDMOccupancyType occupancyType             ,
                                                              Boolean               isSpinRestricted          ,
                                                              Integer               numberFractionalHOOs      ,
                                                              Integer               numberFractionalLUOs      ,
                                                              Integer               numberFractionalAlphaHOOs ,
                                                              Integer               numberFractionalAlphaLUOs ,
                                                              Integer               numberFractionalBetaHOOs  ,
                                                              Integer               numberFractionalBetaLUOs  ,
                                                              Real                  fermiBroadening           )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class QCModelMNDOState:

    cdef CQCModelMNDOState *cObject
    cdef public object      isOwner
