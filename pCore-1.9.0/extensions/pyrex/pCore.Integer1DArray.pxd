#-------------------------------------------------------------------------------
# . File      : pCore.Integer1DArray.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions cimport Boolean, Integer
from pCore.Status       cimport Status, Status_Continue

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "Integer1DArray.h":

    # . The integer array type.
    ctypedef struct CInteger1DArray "Integer1DArray":
        Boolean  isOwner
        Boolean  isView
        Integer  length
        Integer  offset
        Integer  size
        Integer  stride
        Integer *data

    # . Functions.
    cdef CInteger1DArray *Integer1DArray_Allocate   ( Integer length, Status *status )
    cdef CInteger1DArray *Integer1DArray_Clone      ( CInteger1DArray  *self, Status *status )
    cdef void             Integer1DArray_CopyTo     ( CInteger1DArray  *self, CInteger1DArray *other, Status *status )
    cdef void             Integer1DArray_Deallocate ( CInteger1DArray **self )
    cdef Integer          Integer1DArray_GetItem    ( CInteger1DArray  *self, Integer i, Status *status )
    cdef Integer          Integer1DArray_Length     ( CInteger1DArray  *self )
    cdef void             Integer1DArray_Resize     ( CInteger1DArray  *self, Integer length, Integer *initializer, Status *status )
    cdef void             Integer1DArray_Reverse    ( CInteger1DArray  *self )
    cdef void             Integer1DArray_Set        ( CInteger1DArray  *self, Integer value )
    cdef void             Integer1DArray_SetItem    ( CInteger1DArray  *self, Integer i, Integer value, Status *status )
    cdef void             Integer1DArray_Slice      ( CInteger1DArray  *self, Integer start, Integer stop, Integer stride, CInteger1DArray *slice, Status *status )
    cdef void             Integer1DArray_Sort       ( CInteger1DArray  *self )
#    cdef void             Integer1DArray_SortIndex  ( CInteger1DArray  *self, Integer1DArray *indices, Status *status )
    cdef void             Integer1DArray_SortUnique ( CInteger1DArray  *self )
    cdef void             Integer1DArray_ViewOfRaw  ( CInteger1DArray  *self, Integer offset, Integer length, Integer stride, Integer *data, Integer size, Status *status )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class Integer1DArray:

    cdef CInteger1DArray  *cObject
    cdef public object     isOwner
    cdef public object     owner
