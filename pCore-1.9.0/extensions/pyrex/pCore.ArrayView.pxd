#-------------------------------------------------------------------------------
# . File      : pCore.ArrayView.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.SliceOld     cimport CMultiSliceX
from pCore.Status       cimport Status, Status_Continue, Status_Success

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "ArrayView.h":

    ctypedef struct CArrayView "ArrayView":
        Integer  offset 
        Integer  rank   
        Integer  size   
        Integer *extents
        Integer *strides

    ctypedef struct CArrayViewItemIterator "ArrayViewItemIterator":
        pass

    cdef CArrayView *ArrayView_Allocate            ( Integer rank, Status *status )
    cdef CArrayView *ArrayView_AllocateRaw         ( Status *status )
    cdef CArrayView *ArrayView_AllocateWithExtents ( Integer rank, Integer *extents, Status *status )
    cdef Boolean     ArrayView_CheckCapacity       ( CArrayView  *self, Integer capacity )
    cdef CArrayView *ArrayView_Clone               ( CArrayView  *self, Status *status )
    cdef void        ArrayView_Deallocate          ( CArrayView **self )
    cdef Integer     ArrayView_Extent              ( CArrayView  *self, Integer dimension, Status *status )
    cdef Boolean     ArrayView_IsCompact           ( CArrayView  *self )
    cdef Boolean     ArrayView_IsConformable       ( CArrayView  *self, CArrayView *other )
    cdef Boolean     ArrayView_IsUniform           ( CArrayView  *self )
    cdef Integer     ArrayView_ItemIndex           ( CArrayView  *self, Integer rank, Integer *indices )
    cdef Integer     ArrayView_Rank                ( CArrayView  *self )
    cdef Boolean     ArrayView_Reshape1D           ( CArrayView  *self, Integer *offset, Integer *extent, Integer *stride, Status *status )
    cdef Integer     ArrayView_Size                ( CArrayView  *self )
    cdef Status      ArrayView_Slice               ( CArrayView  *self, CMultiSliceX *multiSlice, CArrayView **view )

    cdef CArrayViewItemIterator *ArrayViewItemIterator_Allocate        ( CArrayView *target, Status *status )
    cdef void                    ArrayViewItemIterator_Deallocate      ( CArrayViewItemIterator **self )
    cdef void                    ArrayViewItemIterator_Initialize      ( CArrayViewItemIterator  *self )
    cdef Integer                 ArrayViewItemIterator_Next            ( CArrayViewItemIterator  *self )
    cdef Integer                 ArrayViewItemIterator_NextWithIndices ( CArrayViewItemIterator  *self, Integer *indices )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
#cdef class ArrayView:

#    cdef CArrayView    *cObject
#    cdef public object  isOwner
