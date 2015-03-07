#-------------------------------------------------------------------------------
# . File      : pCore.Status.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions cimport Boolean, CFalse, CTrue, Integer, Real

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "Status.h":

    ctypedef enum Status:
        Status_DimensionError          =  0
        Status_DivideByZero            =  1
        Status_IndexOutOfRange         =  2
        Status_InvalidArgument         =  3
        Status_LinearEquationFailure   =  4
        Status_Null                    =  5
        Status_NullVector              =  6
        Status_OutOfMemory             =  7
        Status_OverFlow                =  8
        Status_SingularMatrix          =  9
        Status_Success                 = 10
        Status_ValueError              = 11
        Status_ArrayDimensionMismatch  = 12
        Status_ArrayLengthMismatch     = 13
        Status_Continue                = 14
        Status_InvalidArrayOperation   = 15
        Status_InvalidDimension        = 16
        Status_MemoryAllocationFailure = 17
        Status_NegativeArrayLength     = 18

cdef extern from "StatusHandler.h":

    ctypedef struct CStatusHandler "StatusHandler":
        pass

    cdef Status  StatusHandler_Allocate   ( CStatusHandler **self )
    cdef Boolean StatusHandler_Continue   ( CStatusHandler  *self )
    cdef void    StatusHandler_Deallocate ( CStatusHandler **self )
    cdef void    StatusHandler_Initialize ( CStatusHandler  *self )
    cdef void    StatusHandler_Record     ( CStatusHandler  *self, Status status )
