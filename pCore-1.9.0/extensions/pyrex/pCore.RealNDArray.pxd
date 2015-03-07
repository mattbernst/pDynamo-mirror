#-------------------------------------------------------------------------------
# . File      : pCore.RealNDArray.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.ArrayView    cimport ArrayView_Allocate, ArrayView_Clone, ArrayView_Slice, ArrayViewItemIterator_Allocate, ArrayViewItemIterator_Deallocate, \
                                ArrayViewItemIterator_Initialize, ArrayViewItemIterator_Next, ArrayViewItemIterator_NextWithIndices, CArrayView, CArrayViewItemIterator
from pCore.cDefinitions cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Memory       cimport Memory_Allocate_Array_Integer, Memory_Deallocate_Integer
from pCore.Real1DArray  cimport CReal1DArray, Real1DArray
from pCore.Real2DArray  cimport CReal2DArray, Real2DArray
from pCore.SliceOld     cimport CMultiSliceX, CSliceX, MultiSliceX_Allocate, MultiSliceX_Deallocate, MultiSliceX_SetRank, SliceX_SetFromScalar, SliceX_SetFromSlice
from pCore.Status       cimport CStatusHandler, Status, Status_Continue, Status_MemoryAllocationFailure, StatusHandler_Continue, StatusHandler_Initialize, StatusHandler_Record

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "RealNDArray.h":

    # . The real array type.
    ctypedef struct CRealNDArray "RealNDArray":
        Boolean     isOwner 
        Boolean     isView  
        Integer     capacity
        CArrayView *view    
        Real       *data    

    # . Functions.
    cdef CRealNDArray *RealNDArray_Allocate    ( Integer rank, Integer *extents, Status *status )
    cdef CRealNDArray *RealNDArray_AllocateRaw ( Status *status )
    cdef CRealNDArray *RealNDArray_Clone       ( CRealNDArray  *self, Boolean doShallow, Status *status )
    cdef void          RealNDArray_CopyTo      ( CRealNDArray  *self, CRealNDArray *other, Status *status )
    cdef void          RealNDArray_Deallocate  ( CRealNDArray **self )
    cdef Integer       RealNDArray_Extent      ( CRealNDArray  *self, Integer dimension, Status *status )
    cdef Real          RealNDArray_GetItem     ( CRealNDArray  *self, Integer *indices, Status *status )
    cdef Boolean       RealNDArray_IsCompact   ( CRealNDArray  *self )
    cdef Boolean       RealNDArray_IsUniform   ( CRealNDArray  *self )
    cdef Integer       RealNDArray_Rank        ( CRealNDArray  *self )
    cdef void          RealNDArray_Reshape1D   ( CRealNDArray  *self, CReal1DArray *view, Status *status )
    cdef void          RealNDArray_Scale       ( CRealNDArray  *self, Real value )
    cdef void          RealNDArray_Set         ( CRealNDArray  *self, Real value )
    cdef void          RealNDArray_SetItem     ( CRealNDArray  *self, Integer *indices, Real value, Status *status )
    cdef Integer       RealNDArray_Size        ( CRealNDArray  *self )
    cdef void          RealNDArray_Slice1D     ( CRealNDArray  *self, CArrayView *view, CReal1DArray *slice, Status *status )
    cdef void          RealNDArray_Slice2D     ( CRealNDArray  *self, CArrayView *view, CReal2DArray *slice, Status *status )
    cdef CRealNDArray *RealNDArray_SliceND     ( CRealNDArray  *self, CArrayView *view, Status *status )
    cdef void          RealNDArray_TailSlice2D ( CRealNDArray  *self, Integer *indices, CReal2DArray *slice, Status *status )
    cdef CRealNDArray *RealNDArray_ViewOfRaw   ( CArrayView *view, Real *data, Integer capacity, Status *status )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class RealNDArray:

    cdef CRealNDArray  *cObject
    cdef Integer       *indices
    cdef public Integer rank
    cdef public object  isOwner
    cdef public object  owner

    cdef _Slice1D ( self, CArrayView *cView )
    cdef _Slice2D ( self, CArrayView *cView )
    cdef _SliceND ( self, CArrayView *cView )
    cdef _ViewOf1DArray ( self, Real1DArray other, CArrayView *cView )
    cdef _ViewOf2DArray ( self, Real2DArray other, CArrayView *cView )
    cdef _ViewOfNDArray ( self, RealNDArray other, CArrayView *cView )
