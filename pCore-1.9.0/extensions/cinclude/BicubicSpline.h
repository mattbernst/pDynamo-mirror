/*------------------------------------------------------------------------------
! . File      : BicubicSpline.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _BICUBICSPLINE
# define _BICUBICSPLINE

# include "Boolean.h"
# include "Integer.h"
# include "Real.h"
# include "Real1DArray.h"
# include "Real2DArray.h"
# include "RealNDArray.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Spline boundary conditions. */
typedef enum {
    BicubicSplineType_Clamped  = 0,
    BicubicSplineType_Natural  = 1,
    BicubicSplineType_NotAKnot = 2,
    BicubicSplineType_Periodic = 3
} BicubicSplineType ;

/* . The bicubic spline type. */
typedef struct {
    BicubicSplineType type    ;
    Integer           lengthx ; /* . Number of points >= 2. Number of intervals is length - 1. */
    Integer           lengthy ;
    Real1DArray      *x       ; /* . X-values. */
    Real1DArray      *y       ; /* . Y-values. */
    Real2DArray      *f       ; /* . Function values. */
    RealNDArray      *coefficients ; /* . Coefficients required for evaluating the spline. */
} BicubicSpline ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern BicubicSpline *BicubicSpline_Allocate            ( const Integer lengthx, const Integer lengthy, const Boolean doX, const Boolean doY, const Boolean doF, Status *status ) ;
extern BicubicSpline *BicubicSpline_Clone               ( const BicubicSpline  *self, Status *status ) ;
extern void           BicubicSpline_Deallocate          (       BicubicSpline **self ) ;
extern void           BicubicSpline_Evaluate            ( const BicubicSpline  *self, const Real x, const Real y, Real *f, Real *g1, Real *g2, Status *status ) ;
extern BicubicSpline *BicubicSpline_MakeFromReal2DArray ( Real1DArray **x, Real1DArray **y, Real2DArray **f, const BicubicSplineType type, Status *status ) ;

# endif
