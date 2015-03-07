#-------------------------------------------------------------------------------
# . File      : pCore.BPData.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Buffer protocol functions."""

from CoreObjects import CLibraryError

#===================================================================================================================================
# . Allocation.
#===================================================================================================================================
cdef BPData *BPData_Allocate ( Integer rank ):
    """Allocate BP data."""
    cdef BPData *self = NULL
    self = < BPData * > stdlib.malloc ( sizeof ( BPData ) )
    if self != NULL:
        ndim = max ( rank, 0 )
        self.buf        = NULL
#        self.obj        = NULL
        self.len        = 0
        self.itemsize   = 0
        self.ndim       = ndim
        self.format     = NULL
        self.shape      = NULL
        self.strides    = NULL
        self.suboffsets = NULL
        if ndim > 0:
            itemsize = sizeof ( Py_ssize_t )
            self.shape      = < Py_ssize_t * > Memory_Allocate_Array ( ndim, itemsize )
            self.strides    = < Py_ssize_t * > Memory_Allocate_Array ( ndim, itemsize )
            self.suboffsets = < Py_ssize_t * > Memory_Allocate_Array ( ndim, itemsize )
            isOK = ( self.shape != NULL ) and ( self.strides != NULL ) and ( self.suboffsets != NULL )
        else:
            isOK = True
    else:
        isOK = False
    if not isOK:
        BPData_Deallocate ( self )
        raise CLibraryError ( "BP data object allocation failure." )
    return self

#===================================================================================================================================
# . Copy to a view.
#===================================================================================================================================
cdef void BPData_CopyToView ( BPData *self, Py_buffer *view ):
    """Copy the BP data to a view."""
    view.buf        = self.buf
    view.format     = self.format
    view.internal   = NULL
    view.itemsize   = self.itemsize
    view.len        = self.len
    view.ndim       = self.ndim
#    view.obj        = NULL #self.obj
    view.readonly   = 0
    view.shape      = self.shape
    view.strides    = self.strides
    view.suboffsets = self.suboffsets

#===================================================================================================================================
# . Deallocation.
#===================================================================================================================================
cdef BPData *BPData_Deallocate ( BPData *self ):
    """Deallocate a BP data."""
    if self != NULL:
        stdlib.free ( self.shape      )
        stdlib.free ( self.strides    )
        stdlib.free ( self.suboffsets )
        stdlib.free ( self )
    return NULL

#===================================================================================================================================
# . Initialize a view.
#===================================================================================================================================
cdef void BPData_InitializeView ( Py_buffer *view, char *format, Py_ssize_t itemsize ):
    """Initialize a view."""
    view.buf        = NULL
    view.format     = format
    view.internal   = NULL
    view.itemsize   = itemsize
    view.len        = 0
    view.ndim       = 0
#    view.obj        = NULL #obj
    view.readonly   = 0
    view.shape      = NULL
    view.strides    = NULL
    view.suboffsets = NULL

#===================================================================================================================================
# . Set the object fields.
#===================================================================================================================================
cdef void BPData_Set ( BPData *self, void *buf, char *format, Py_ssize_t itemsize, Py_ssize_t size, Integer *shape, Integer *strides ):
    """Set the BP data."""
    # . ndim is already set.
    cdef Integer i
    self.buf      = buf
    self.format   = format
    self.itemsize = itemsize
    self.len      = size * itemsize
#    self.obj      = NULL
    for i from 0 <= i < self.ndim:
        self.shape     [i] = shape   [i]
        self.strides   [i] = strides [i] * itemsize
        self.suboffsets[i] = -1
