#-------------------------------------------------------------------------------
# . File      : pCore.BicubicSpline.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions cimport Boolean, Integer, Real
from pCore.Real1DArray  cimport CReal1DArray
from pCore.Real2DArray  cimport CReal2DArray
from pCore.RealNDArray  cimport CRealNDArray
from pCore.Status       cimport Status

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "CubicSpline.h":

    # . Spline boundary conditions.
    ctypedef enum BicubicSplineType:
        BicubicSplineType_Clamped  = 0
        BicubicSplineType_Natural  = 1
        BicubicSplineType_NotAKnot = 2
        BicubicSplineType_Periodic = 3

    # . The bicubic spline type.
    ctypedef struct CBicubicSpline "BicubicSpline":
        BicubicSplineType type
        Integer           lengthx
        Integer           lengthy
        CReal1DArray      *x
        CReal1DArray      *y
        CReal2DArray      *f
        CRealNDArray      *coefficients

    # . Functions.
    cdef CBicubicSpline *BicubicSpline_Allocate            ( Integer lengthx, Integer lengthy, Boolean doX, Boolean doY, Boolean doF, Status *status )
    cdef CBicubicSpline *BicubicSpline_Clone               ( CBicubicSpline  *self, Status *status )
    cdef void            BicubicSpline_Deallocate          ( CBicubicSpline **self )
    cdef void            BicubicSpline_Evaluate            ( CBicubicSpline  *self, Real x, Real y, Real *f, Real *g1, Real *g2, Status *status )
    cdef CBicubicSpline *BicubicSpline_MakeFromReal2DArray ( CReal1DArray **x, CReal1DArray **y, CReal2DArray **f, BicubicSplineType type, Status *status )
