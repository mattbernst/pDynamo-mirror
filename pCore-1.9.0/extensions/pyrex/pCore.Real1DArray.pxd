#-------------------------------------------------------------------------------
# . File      : pCore.Real1DArray.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from libc cimport stdlib

from pCore.BPData       cimport BPData, BPData_Allocate, BPData_CopyToView, BPData_Deallocate, BPData_InitializeView, BPData_Set
from pCore.cDefinitions cimport Boolean, Integer, Real
from pCore.Memory       cimport Memory_Allocate_Array
from pCore.Status       cimport Status, Status_Continue

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "Real1DArray.h":

    # . The real array type.
    ctypedef struct CReal1DArray "Real1DArray":
        Boolean  isOwner
        Boolean  isView
        Integer  length
        Integer  offset
        Integer  size
        Integer  stride
        Real    *data

    # . Functions.
    cdef Real          Real1DArray_AbsoluteMaximum      ( CReal1DArray  *self )
    cdef Integer       Real1DArray_AbsoluteMaximumIndex ( CReal1DArray  *self )
    cdef void          Real1DArray_AddScalar            ( CReal1DArray  *self, Real value )
    cdef void          Real1DArray_AddScaledArray       ( CReal1DArray  *self, Real value, CReal1DArray *other, Status *status )
    cdef CReal1DArray *Real1DArray_Allocate             ( Integer length, Status *status )
    cdef CReal1DArray *Real1DArray_Clone                ( CReal1DArray  *self, Status *status )
    cdef void          Real1DArray_CopyTo               ( CReal1DArray  *self, CReal1DArray *other, Status *status )
    cdef void          Real1DArray_Deallocate           ( CReal1DArray **self )
    cdef void          Real1DArray_Divide               ( CReal1DArray  *self, CReal1DArray *other, Status *status )
    cdef Real          Real1DArray_Dot                  ( CReal1DArray  *self, CReal1DArray *other, Status *status )
    cdef void          Real1DArray_Exp                  ( CReal1DArray  *self )
    cdef Real          Real1DArray_GetItem              ( CReal1DArray  *self, Integer i, Status *status )
    cdef Integer       Real1DArray_Length               ( CReal1DArray  *self )
    cdef void          Real1DArray_Ln                   ( CReal1DArray  *self )
    cdef void          Real1DArray_Multiply             ( CReal1DArray  *self, CReal1DArray *other, Status *status )
    cdef Real          Real1DArray_Norm2                ( CReal1DArray  *self )
    cdef void          Real1DArray_Normalize            ( CReal1DArray  *self, Real *tolerance, Status *status )
    cdef void          Real1DArray_Resize               ( CReal1DArray  *self, Integer length, Real *initializer, Status *status )
    cdef void          Real1DArray_Reciprocate          ( CReal1DArray  *self )
    cdef void          Real1DArray_Reverse              ( CReal1DArray  *self )
    cdef Real          Real1DArray_RootMeanSquare       ( CReal1DArray  *self )
    cdef void          Real1DArray_Scale                ( CReal1DArray  *self, Real scale )
    cdef void          Real1DArray_Set                  ( CReal1DArray  *self, Real value )
    cdef void          Real1DArray_SetItem              ( CReal1DArray  *self, Integer i, Real value, Status *status )
    cdef void          Real1DArray_Slice                ( CReal1DArray  *self, Integer start, Integer stop, Integer stride, CReal1DArray *slice, Status *status )
    cdef void          Real1DArray_Sort                 ( CReal1DArray  *self )
#    cdef void          Real1DArray_SortIndex            ( CReal1DArray  *self, Integer1DArray *indices, Status *status )
    cdef void          Real1DArray_SortUnique           ( CReal1DArray  *self )
    cdef Real          Real1DArray_Sum                  ( CReal1DArray  *self )
    cdef void          Real1DArray_ViewOfRaw            ( CReal1DArray  *self, Integer offset, Integer length, Integer stride, Real *data, Integer size, Status *status )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class Real1DArray:

    cdef BPData        *bpData
    cdef CReal1DArray  *cObject
    cdef public object  isOwner
    cdef public object  owner
