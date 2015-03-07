#-------------------------------------------------------------------------------
# . File      : pCore.SliceOld.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Status       cimport Status

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "SliceOld.h":

    ctypedef struct CSliceX "SliceX":
        Boolean isScalar
        Integer extent  
        Integer start   
        Integer stop    
        Integer stride  

    ctypedef struct CMultiSliceX "MultiSliceX":
        Integer   rank
        Integer   size
        CSliceX  *items

    cdef Status  SliceX_Allocate      ( CSliceX **self )
    cdef void    SliceX_Deallocate    ( CSliceX **self )
    cdef void    SliceX_Initialize    ( CSliceX  *self )
    cdef Status  SliceX_SetFromScalar ( CSliceX  *self, Integer index, Integer rawExtent )
    cdef Status  SliceX_SetFromSlice  ( CSliceX  *self, Integer start, Integer stop, Integer stride, Integer rawExtent )

    cdef Status  MultiSliceX_Allocate    ( CMultiSliceX **self, Integer rank )
    cdef Status  MultiSliceX_AllocateRaw ( CMultiSliceX **self )
    cdef void    MultiSliceX_Deallocate  ( CMultiSliceX **self )
    cdef void    MultiSliceX_SetRank     ( CMultiSliceX  *self )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
#cdef class MultiSlice:

#    cdef CMultiSlice   *cObject
#    cdef public object  isOwner
