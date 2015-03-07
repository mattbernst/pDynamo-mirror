#-------------------------------------------------------------------------------
# . File      : pCore.Real2DArray.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Real1DArray  cimport CReal1DArray, Real1DArray, Real1DArray_ViewOfRaw
from pCore.Status       cimport Status, Status_Continue

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "Real2DArray.h":

    # . The array type.
    ctypedef struct CReal2DArray "Real2DArray":
        Boolean  isOwner
        Boolean  isView
        Integer  length
        Integer  length0
        Integer  length1
        Integer  offset
        Integer  size
        Integer  stride0
        Integer  stride1
        Real    *data

    # . Functions.
    cdef void          Real2DArray_1DSlice             ( CReal2DArray  *self, Integer start0, Integer stop0, Integer stride0, Integer start1, Integer stop1, Integer stride1, CReal1DArray *view, Status *status )
    cdef CReal2DArray *Real2DArray_Allocate            ( Integer length0, Integer length1, Status *status )
    cdef CReal2DArray *Real2DArray_Clone               ( CReal2DArray  *self, Status *status )
#    cdef void          Real2DArray_ColumnSlice         ( CReal2DArray  *self, Integer column, Real1DArray *slice, Status *status )
    cdef void          Real2DArray_CopyTo              ( CReal2DArray  *self, CReal2DArray *other, Status *status )
    cdef void          Real2DArray_Deallocate          ( CReal2DArray **self )
    cdef Real          Real2DArray_Determinant         ( CReal2DArray  *self, Status *status )
    cdef void          Real2DArray_Exp                 ( CReal2DArray  *self )
    cdef Real          Real2DArray_GetItem             ( CReal2DArray  *self, Integer i, Integer j, Status *status )
    cdef Integer       Real2DArray_Length              ( CReal2DArray  *self, Integer dimension )
    cdef void          Real2DArray_MatrixMultiply      ( Boolean aTranspose, Boolean bTranspose, Real alpha, CReal2DArray *a, CReal2DArray *b, Real beta, CReal2DArray *c, Status *status )
    cdef void          Real2DArray_ProjectOutOf1DArray ( CReal2DArray  *self, CReal1DArray *vector, Status *status )
    cdef void          Real2DArray_Resize              ( CReal2DArray  *self, Integer length0, Real *initializer, Status *status )
#    cdef Real          Real2DArray_RootMeanSquare      ( CReal1DArray  *self )
#    cdef void          Real2DArray_RowSlice            ( CReal2DArray  *self, Integer row, Real1DArray *slice, Status *status )
    cdef void          Real2DArray_Scale               ( CReal2DArray  *self, Real value )
    cdef void          Real2DArray_Set                 ( CReal2DArray  *self, Real value )
    cdef void          Real2DArray_SetItem             ( CReal2DArray  *self, Integer i, Integer j, Real value, Status *status )
    cdef void          Real2DArray_Slice               ( CReal2DArray  *self, Integer start0, Integer stop0, Integer stride0, Integer start1, Integer stop1, Integer stride1, CReal2DArray *view, Status *status )
    cdef void          Real2DArray_Transpose           ( CReal2DArray  *self, Status *status )
    cdef void          Real2DArray_VectorMultiply      ( Boolean aTranspose, Real alpha, CReal2DArray *a, CReal1DArray *x, Real beta, CReal1DArray *y, Status *status )
    cdef void          Real2DArray_ViewOfRaw           ( CReal2DArray  *self, Integer offset, Integer length0, Integer stride0, Integer length1, Integer stride1, Real *data, Integer size, Status *status )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class Real2DArray:

    cdef CReal2DArray  *cObject
    cdef public object  isOwner
    cdef public object  owner
