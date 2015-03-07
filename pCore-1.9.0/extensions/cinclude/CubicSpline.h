/*------------------------------------------------------------------------------
! . File      : CubicSpline.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _CUBICSPLINE
# define _CUBICSPLINE

# include "Boolean.h"
# include "Integer.h"
# include "Real.h"
# include "Real1DArray.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The cubic spline type. */
typedef struct {
    Integer length ; /* . Number of points >= 2. Number of intervals is length - 1. */
    Real1DArray *h ;   /* . Second derivatives. */
    Real1DArray *x ;   /* . X-values. */
    Real1DArray *y ;   /* . Y-values. */
} CubicSpline ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define CubicSpline_FastEvaluateFG( self, l, u, d, s, t, f, g ) \
    { \
        auto Real hl, hu, yl, yu ; \
        hl = Real1DArray_Item ( self->h, l ) * d / 6.0e+00 ; \
        hu = Real1DArray_Item ( self->h, u ) * d / 6.0e+00 ; \
        yl = Real1DArray_Item ( self->y, l ) ; \
        yu = Real1DArray_Item ( self->y, u ) ; \
        f  = t * yl + s * yu + d * ( t * ( t * t - 1.0e+00 ) * hl + s * ( s * s - 1.0e+00 ) * hu ) ; \
        g  = ( yu - yl ) / d + ( - ( 3.0e+00 * t * t - 1.0e+00 ) * hl + ( 3.0e+00 * s * s - 1.0e+00 ) * hu ) ; \
    }

# define CubicSpline_FastEvaluateLUDST( self, x0, l, u, d, s, t ) \
    { \
        auto Integer i ; \
        auto Real1DArray *abscissa = self->x ; \
        l = 0 ; \
        u = self->length - 1 ; \
        while ( ( u - l ) > 1 ) \
        { \
            i = ( u + l ) >> 1 ; \
            if ( Real1DArray_Item ( abscissa, i ) > x0 ) u = i ; \
            else                                         l = i ; \
        } \
        d = ( Real1DArray_Item ( abscissa, u ) - Real1DArray_Item ( abscissa, l ) ) ; \
        s = ( x0 - Real1DArray_Item ( abscissa, l ) ) / d ; \
        t = ( Real1DArray_Item ( abscissa, u ) - x0 ) / d ; \
    }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern Status CubicSpline_Allocate             (       CubicSpline **self, const Integer length, const Boolean QX, const Boolean QY ) ;
extern Status CubicSpline_Clone                ( const CubicSpline  *self, CubicSpline **clone ) ;
extern void   CubicSpline_Deallocate           (       CubicSpline **self ) ;
extern Status CubicSpline_Evaluate             ( const CubicSpline  *self, const Real x, Real *f, Real *g, Real *h ) ;
extern void   CubicSpline_EvaluateLUDST        ( const CubicSpline  *self, const Real x, Integer *l, Integer *u, Real *d, Real *s, Real *t ) ;
extern Status CubicSpline_FindExtrema          ( const CubicSpline  *self, Real1DArray **maxima, Real1DArray **minima ) ;
extern Real   CubicSpline_Integrate            ( const CubicSpline  *self, const Real a, const Real b, Status *status ) ;
extern Real   CubicSpline_IntegrateFull        ( const CubicSpline  *self, Status *status ) ;
extern Status CubicSpline_MakeFromReal1DArrays (       CubicSpline **self, Real1DArray **x, Real1DArray **y, const Integer lowerderivative, const Real lowervalue, const Integer upperderivative, const Real uppervalue ) ;

# endif
