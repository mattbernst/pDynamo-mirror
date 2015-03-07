#-------------------------------------------------------------------------------
# . File      : pMolecule.QCMMInteractionState.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions    cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Real1DArray     cimport CReal1DArray, Real1DArray , Real1DArray_Clone
from pCore.Status          cimport Status
from pCore.SymmetricMatrix cimport CSymmetricMatrix

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "QCMMInteractionState.h":

    ctypedef struct CQCMMInteractionState "QCMMInteractionState":
        CReal1DArray     *qcCharges
        CReal1DArray     *qcmmPotentials
        CSymmetricMatrix *qcqcPotentials

    cdef CQCMMInteractionState *QCMMInteractionState_Allocate               ( Integer n, Boolean includeQCQC, Status *status )
    cdef void                   QCMMInteractionState_AllocateQCQCPotentials ( CQCMMInteractionState  *self, Status *status )
    cdef void                   QCMMInteractionState_Deallocate             ( CQCMMInteractionState **self )
    cdef void                   QCMMInteractionState_Initialize             ( CQCMMInteractionState  *self )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class QCMMInteractionState:

    cdef CQCMMInteractionState *cObject
    cdef public object          isOwner
