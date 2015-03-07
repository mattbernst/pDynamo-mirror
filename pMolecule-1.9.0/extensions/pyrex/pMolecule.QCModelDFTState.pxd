#-------------------------------------------------------------------------------
# . File      : pMolecule.QCModelDFTState.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.BlockStorage             cimport BlockStorage_Size, CBlockStorage
from pCore.cDefinitions             cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Coordinates3             cimport CCoordinates3, Coordinates3
from pCore.Real1DArray              cimport CReal1DArray, Real1DArray
from pCore.Real2DArray              cimport CReal2DArray, Real2DArray, Real2DArray_Length
from pCore.SymmetricMatrix          cimport SymmetricMatrix, CSymmetricMatrix
from pMolecule.DFTGrid              cimport CDFTGrid, DFTGrid_FunctionDataSize, DFTGrid_HasFunctionData, DFTGrid_NumberOfFunctionValues
from pMolecule.QCAtomContainer      cimport QCAtomContainer, QCAtomContainer_Size, CQCAtomContainer
from pMolecule.QCMMInteractionState cimport QCMMInteractionState, CQCMMInteractionState
from pMolecule.QCOnePDM             cimport CQCOnePDM, QCOnePDM, QCOnePDM_SpinExpectationValues, QCOnePDMOccupancyType_FractionalVariable, QCOnePDMOccupancyType
from pMolecule.QCParameters         cimport QCParameters, CQCParameter

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "QCModelDFTState.h":

    ctypedef struct CQCModelDFTState "QCModelDFTState":
        Boolean               isConverged
        Real                  eelectronic
        Real                  enuclear
        Real                  eocc
        Real                  eoei
        Real                  eqcmm
        Real                  eqcqc
        Real                  equad
        Real                  etei
        Real                  rhoquad
        Integer               ncycles
        CQCAtomContainer      *qcAtoms
        CQCParameter          *qcParameters
        CCoordinates3         *coordinates3
        CCoordinates3         *gradients3
        CQCMMInteractionState *qcmmstate
        CReal2DArray          *c2o
        CReal2DArray          *o2c
        CQCOnePDM             *densityp
        CQCOnePDM             *densityq
        CBlockStorage         *fitintegrals
        CCoordinates3         *qccoordinates3
        CDFTGrid              *dftgrid
        CReal2DArray          *lowdintransformation
        CReal2DArray          *orthogonalizer
        CReal2DArray          *overlapeigenvectors
        CSymmetricMatrix      *inversefitmatrix
        CSymmetricMatrix      *oneelectronmatrix
        CSymmetricMatrix      *overlap
        CSymmetricMatrix      *wdensity
        CReal1DArray          *fpotential
        CReal1DArray          *fselfoverlap
        CReal1DArray          *overlapeigenvalues
        CReal1DArray          *wvector

    cdef CQCModelDFTState *QCModelDFTState_Allocate   ( )
    cdef void              QCModelDFTState_Deallocate (       CQCModelDFTState **self )
    cdef void              QCModelDFTState_Finalize   (       CQCModelDFTState  *self, Boolean QKEEPDATA )
    cdef void              QCModelDFTState_Initialize (       CQCModelDFTState  *self, CCoordinates3 *coordinates3, CQCMMInteractionState *qcmmstate, CCoordinates3 *gradients3 )
    cdef CQCModelDFTState *QCModelDFTState_Setup      ( CQCAtomContainer *qcAtoms, CQCParameter *qcParameters, Real alphaCharge, Real betaCharge, QCOnePDMOccupancyType occupancyType, Boolean isSpinRestricted, \
                                                                                                                                                    Integer numberFractionalHOOs     , Integer numberFractionalLUOs     , \
                                                                                                                                                    Integer numberFractionalAlphaHOOs, Integer numberFractionalAlphaLUOs, \
                                                                                                                                                    Integer numberFractionalBetaHOOs , Integer numberFractionalBetaLUOs , \
                                                                                                                                                                                           Real fermiBroadening )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class QCModelDFTState:

    cdef CQCModelDFTState *cObject
    cdef public object     isOwner
