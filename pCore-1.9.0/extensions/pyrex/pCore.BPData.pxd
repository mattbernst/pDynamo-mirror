#-------------------------------------------------------------------------------
# . File      : pCore.BPData.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from libc cimport stdlib

from pCore.cDefinitions cimport Boolean, Integer, Real
from pCore.Memory       cimport Memory_Allocate_Array
from pCore.Status       cimport Status, Status_Continue

#===================================================================================================================================
# . Structures.
#===================================================================================================================================
# . The buffer type.
ctypedef struct BPData:
    void       *buf
#    PyObject   *obj
    Py_ssize_t  len
    Py_ssize_t  itemsize
    int         ndim
    char       *format
    Py_ssize_t *shape
    Py_ssize_t *strides
    Py_ssize_t *suboffsets

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
cdef BPData *BPData_Allocate       ( Integer rank )
cdef void    BPData_CopyToView     ( BPData    *self, Py_buffer *view )
cdef BPData *BPData_Deallocate     ( BPData    *self )
cdef void    BPData_InitializeView ( Py_buffer *view, char *format, Py_ssize_t itemsize )#, object self )
cdef void    BPData_Set            ( BPData    *self, void *buf, char *format, Py_ssize_t itemsize, Py_ssize_t size, Integer *shape, Integer *strides )
