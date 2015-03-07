#-------------------------------------------------------------------------------
# . File      : pCore.Memory.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions cimport Boolean, CFalse, CSize, CTrue, Integer, Real

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "Memory.h":

    # . Functions.
    cdef void    *Memory_Allocate_Array                    ( CSize nelements, CSize nsize )
    cdef Boolean *Memory_Allocate_Array_Boolean            ( CSize nelements )
    cdef Boolean *Memory_Allocate_Array_Boolean_Initialize ( CSize nelements, Boolean value )
    cdef Real    *Memory_Allocate_Array_Real               ( CSize nelements )
    cdef Real    *Memory_Allocate_Array_Real_Initialize    ( CSize nelements, Real    value )
    cdef Integer *Memory_Allocate_Array_Integer            ( CSize nelements )
    cdef Integer *Memory_Allocate_Array_Integer_Initialize ( CSize nelements, Integer value )
    cdef void     Memory_Deallocate_Integer                ( Integer **address )
    cdef void     Memory_Deallocate_Real                   ( Real    **address )
