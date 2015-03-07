#-------------------------------------------------------------------------------
# . File      : pCore.Integer2DArray.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions    cimport Boolean, Integer
from pCore.Integer1DArray  cimport CInteger1DArray, Integer1DArray, Integer1DArray_ViewOfRaw
from pCore.Status          cimport Status, Status_Continue

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "Integer2DArray.h":

    # . The array type.
    ctypedef struct CInteger2DArray "Integer2DArray":
        Boolean  isOwner
        Boolean  isView
        Integer  length
        Integer  length0
        Integer  length1
        Integer  offset
        Integer  size
        Integer  stride0
        Integer  stride1
        Integer *data

    # . Functions.
    cdef void             Integer2DArray_1DSlice    ( CInteger2DArray  *self, Integer start0, Integer stop0, Integer stride0, Integer start1, Integer stop1, Integer stride1, CInteger1DArray *view, Status *status )
    cdef CInteger2DArray *Integer2DArray_Allocate   ( Integer length0, Integer length1, Status *status )
    cdef CInteger2DArray *Integer2DArray_Clone      ( CInteger2DArray  *self, Status *status )
    cdef void             Integer2DArray_CopyTo     ( CInteger2DArray  *self, CInteger2DArray *other, Status *status )
    cdef void             Integer2DArray_Deallocate ( CInteger2DArray **self )
    cdef Integer          Integer2DArray_GetItem    ( CInteger2DArray  *self, Integer i, Integer j, Status *status )
    cdef Integer          Integer2DArray_Length     ( CInteger2DArray  *self, Integer dimension )
    cdef void             Integer2DArray_Resize     ( CInteger2DArray  *self, Integer length0, Integer *initializer, Status *status )
    cdef void             Integer2DArray_Set        ( CInteger2DArray  *self, Integer value )
    cdef void             Integer2DArray_SetItem    ( CInteger2DArray  *self, Integer i, Integer j, Integer value, Status *status )
    cdef void             Integer2DArray_Slice      ( CInteger2DArray  *self, Integer start0, Integer stop0, Integer stride0, Integer start1, Integer stop1, Integer stride1, CInteger2DArray *view, Status *status )
    cdef void             Integer2DArray_ViewOfRaw  ( CInteger2DArray  *self, Integer offset, Integer length0, Integer stride0, Integer length1, Integer stride1, Integer *data, Integer size, Status *status )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class Integer2DArray:

    cdef CInteger2DArray *cObject
    cdef public object    isOwner
    cdef public object    owner
