#-------------------------------------------------------------------------------
# . File      : pCore.CubicSpline.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions cimport Boolean, Integer, Real
from pCore.Real1DArray  cimport CReal1DArray, Real1DArray_Allocate, Real1DArray_Deallocate
from pCore.Status       cimport Status, Status_Null

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "CubicSpline.h":

    # . The cubic spline type.
    ctypedef struct CCubicSpline "CubicSpline":
        Integer length
        CReal1DArray *h
        CReal1DArray *x
        CReal1DArray *y

    # . Functions.
    cdef Status CubicSpline_Allocate             ( CCubicSpline **self, Integer length, Boolean QX, Boolean QY )
    cdef Status CubicSpline_Clone                ( CCubicSpline  *self, CCubicSpline **clone )
    cdef void   CubicSpline_Deallocate           ( CCubicSpline **self )
    cdef Status CubicSpline_MakeFromReal1DArrays ( CCubicSpline **self, CReal1DArray **x, CReal1DArray **y, Integer lowerderivative, Real lowervalue, Integer upperderivative, Real uppervalue )
    cdef Status CubicSpline_Evaluate             ( CCubicSpline  *self, Real x, Real *f, Real *g, Real *h )
    cdef Status CubicSpline_FindExtrema          ( CCubicSpline  *self, CReal1DArray **maxima, CReal1DArray **minima )
    cdef Real   CubicSpline_Integrate            ( CCubicSpline  *self, Real a, Real b, Status *status )
    cdef Real   CubicSpline_IntegrateFull        ( CCubicSpline  *self, Status *status )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class CubicSpline:

    cdef CCubicSpline  *cObject
    cdef public object  isOwner
