/*------------------------------------------------------------------------------
! . File      : BicubicSpline.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . Bicubic splines.
! . Abscissae should be in strictly ascending order.
!=================================================================================================================================*/

# include <stdio.h>
# include <stdlib.h>

# include "BicubicSpline.h"
# include "math.h"
# include "Memory.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void    BicubicSpline_EvaluateCoefficientTable     ( BicubicSpline *self, const Real2DArray *p, const Real2DArray *q, const Real2DArray *r ) ;
static void    BicubicSpline_Get1DDerivatives             ( const Real1DArray *x, const Real1DArray *u, Real1DArray *d, BicubicSplineType type, Real1DArray *Ad, Real1DArray *Asd, Real1DArray *qdy, Real1DArray *lll ) ;
static Integer BicubicSpline_Locate                       ( const Real1DArray *abscissa, const Real x, const Boolean isPeriodic, Status *status ) ;
static void    BicubicSpline_Setup                        ( BicubicSpline *self, Status *status ) ;
static void    BicubicSpline_TridiagonalLDLtSolve         ( const Integer n, Real1DArray *d, Real1DArray *l, Real1DArray *b ) ;
static void    BicubicSpline_TridiagonalLDLtSolvePeriodic ( const Integer n, Real1DArray *d, Real1DArray *lsd, Real1DArray *lll, Real1DArray *b ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
BicubicSpline *BicubicSpline_Allocate ( const Integer lengthx, const Integer lengthy, const Boolean doX, const Boolean doY, const Boolean doF, Status *status )
{
    BicubicSpline *self = NULL ;
    if ( ( lengthx >= 2 ) && ( lengthy >= 2 ) )
    {
        MEMORY_ALLOCATE ( self, BicubicSpline ) ;
        if ( self != NULL )
        {
            auto Integer lengths[4] ;
            auto Status  localstatus   ;
            Status_Set ( &localstatus, Status_Continue ) ;
            self->type         = BicubicSplineType_Natural ;
            self->lengthx      = lengthx ;
            self->lengthy      = lengthy ;
            self->x            = NULL    ;
            self->y            = NULL    ;
            self->f            = NULL    ;
            self->coefficients = NULL    ;
            /* . Array allocation. */
            lengths[0] = lengthx - 1 ; lengths[1] = lengthy - 1 ; lengths[2] = 4 ; lengths[3] = 4 ;
            self->coefficients = RealNDArray_Allocate ( 4, lengths, &localstatus ) ;
            if ( doX ) self->x = Real1DArray_Allocate ( lengthx, &localstatus ) ;
            if ( doY ) self->y = Real1DArray_Allocate ( lengthy, &localstatus ) ;
            if ( doF ) self->f = Real2DArray_Allocate ( lengthx, lengthy, &localstatus ) ;
            if ( ! Status_OK ( &localstatus ) ) BicubicSpline_Deallocate ( &self ) ;
        }
        if ( self == NULL ) Status_Set ( status, Status_OutOfMemory ) ;
    }
    else Status_Set ( status, Status_DimensionError ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
BicubicSpline *BicubicSpline_Clone ( const BicubicSpline *self, Status *status )
{
    BicubicSpline *clone = NULL ;
    if ( self != NULL )
    {
        clone = BicubicSpline_Allocate ( self->lengthx, self->lengthy, True, True, True, status ) ;
        if ( clone != NULL )
        {
            clone->type = self->type ;
            Real1DArray_CopyTo ( self->x, clone->x, status ) ;
            Real1DArray_CopyTo ( self->y, clone->y, status ) ;
            Real2DArray_CopyTo ( self->f, clone->f, status ) ;
            RealNDArray_CopyTo ( self->coefficients, clone->coefficients, status ) ;
        }
    }
    return clone ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void BicubicSpline_Deallocate ( BicubicSpline **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) )
    {
        Real1DArray_Deallocate ( &((*self)->x) ) ;
        Real1DArray_Deallocate ( &((*self)->y) ) ;
        Real2DArray_Deallocate ( &((*self)->f) ) ;
        RealNDArray_Deallocate ( &((*self)->coefficients) ) ;
        MEMORY_DEALLOCATE ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Evaluation (function and first derivatives only).
!---------------------------------------------------------------------------------------------------------------------------------*/
void BicubicSpline_Evaluate ( const BicubicSpline *self, const Real x, const Real y, Real *f, Real *g1, Real *g2, Status *status )
{
    if ( f  != NULL ) (*f)  = 0.0e+00 ;
    if ( g1 != NULL ) (*g1) = 0.0e+00 ;
    if ( g2 != NULL ) (*g2) = 0.0e+00 ;
    if ( self != NULL )
    {
        auto Integer ix, iy ;

        /* . Locate the indices. */
        ix = BicubicSpline_Locate ( self->x, x, ( self->type == BicubicSplineType_Periodic ), status ) ;
        iy = BicubicSpline_Locate ( self->y, y, ( self->type == BicubicSplineType_Periodic ), status ) ;

        /* . Compute the function and derivative values. */
        if ( ( ix >= 0 ) && ( iy >= 0 ) )
        {
            auto Integer     i, indices[2] ;
            auto Real        dudx, dudy, dx, dy, u ;
            auto Real2DArray C ;

            /* . Factors. */
            indices[0] = ix ; indices[1] = iy ;
            RealNDArray_TailSlice2D ( self->coefficients, indices, &C, status ) ;
            dx   = x - Real1DArray_Item ( self->x, ix ) ;
            dy   = y - Real1DArray_Item ( self->y, iy ) ;
            u    = 0.0e+00 ;
            dudx = 0.0e+00 ;
            dudy = 0.0e+00 ;
            for ( i = 3 ; i >= 0 ; i-- )
            {
                u    = Real2DArray_Item ( &C, i, 0 ) + dy * ( Real2DArray_Item ( &C, i, 1 ) + dy * ( Real2DArray_Item ( &C, i, 2 ) + dy * Real2DArray_Item ( &C, i, 3 ) ) ) + u * dx ;
                dudx = Real2DArray_Item ( &C, 1, i ) + dx * ( 2.0e+00 * Real2DArray_Item ( &C, 2, i ) + 3.0e+00 * dx * Real2DArray_Item ( &C, 3, i ) ) + dudx * dy ;
                dudy = Real2DArray_Item ( &C, i, 1 ) + dy * ( 2.0e+00 * Real2DArray_Item ( &C, i, 2 ) + 3.0e+00 * dy * Real2DArray_Item ( &C, i, 3 ) ) + dudy * dx ;
            }
            if ( f  != NULL ) (*f)  = u ;
            if ( g1 != NULL ) (*g1) = dudx ;
            if ( g2 != NULL ) (*g2) = dudy ;
        }
        else Status_Set ( status, Status_IndexOutOfRange ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Evaluate the coefficient table.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void BicubicSpline_EvaluateCoefficientTable ( BicubicSpline *self, const Real2DArray *p, const Real2DArray *q, const Real2DArray *r )
{
    auto Integer     i, j, indices[2] ;
    auto Real        a, b, c, d, dx, dy ;
    auto Real2DArray C, *u ;
    u = self->f ;
    for ( j = 0 ; j < self->lengthy - 1 ; j++ )
    {
        dy = 1.0e+00 / ( Real1DArray_Item ( self->y, j+1 ) - Real1DArray_Item ( self->y, j ) ) ;
        for ( i = 0 ; i < self->lengthx - 1 ; i++ )
        {
            dx = 1.0e+00 / ( Real1DArray_Item ( self->x, i+1 ) - Real1DArray_Item ( self->x, i ) ) ;
            indices[0] = i ; indices[1] = j ;
            RealNDArray_TailSlice2D ( self->coefficients, indices, &C, NULL ) ;

            Real2DArray_Item ( &C, 0, 0 ) = Real2DArray_Item ( u, i, j ) ;
            Real2DArray_Item ( &C, 1, 0 ) = Real2DArray_Item ( p, i, j ) ;
            Real2DArray_Item ( &C, 0, 1 ) = Real2DArray_Item ( q, i, j ) ;
            Real2DArray_Item ( &C, 1, 1 ) = Real2DArray_Item ( r, i, j ) ;

            a = ( Real2DArray_Item ( u, i+1, j ) - Real2DArray_Item ( u, i, j ) ) * dx ;
            Real2DArray_Item ( &C, 2, 0 ) = ( 3.0e+00 * a - 2.0e+00 * Real2DArray_Item ( p, i, j ) - Real2DArray_Item ( p, i+1, j ) ) * dx ;
            Real2DArray_Item ( &C, 3, 0 ) = ( Real2DArray_Item ( p, i+1, j ) + Real2DArray_Item ( p, i, j ) - 2.0e+00 * a ) * ( dx * dx ) ;

            a = ( Real2DArray_Item ( u, i, j+1 ) - Real2DArray_Item ( u, i, j ) ) * dy ;
            Real2DArray_Item ( &C, 0, 2 ) = ( 3.0e+00 * a - 2.0e+00 * Real2DArray_Item ( q, i, j ) - Real2DArray_Item ( q, i, j+1 ) ) * dy ;
            Real2DArray_Item ( &C, 0, 3 ) = ( Real2DArray_Item ( q, i, j+1 ) + Real2DArray_Item ( q, i, j ) - 2.0e+00 * a ) * ( dy * dy ) ;

            a = ( Real2DArray_Item ( q, i+1, j ) - Real2DArray_Item ( q, i, j ) ) * dx ;
            Real2DArray_Item ( &C, 2, 1 ) = ( 3.0e+00 * a - Real2DArray_Item ( r, i+1, j ) - 2.0e+00 * Real2DArray_Item ( r, i, j ) ) * dx ;
            Real2DArray_Item ( &C, 3, 1 ) = ( Real2DArray_Item ( r, i+1, j ) + Real2DArray_Item ( r, i, j ) - 2.0e+00 * a ) * ( dx * dx ) ;

            a = ( Real2DArray_Item ( p, i, j+1 ) - Real2DArray_Item ( p, i, j ) ) * dy ;
            Real2DArray_Item ( &C, 1, 2 ) = ( 3.0e+00 * a - Real2DArray_Item ( r, i, j+1 ) - 2.0e+00 * Real2DArray_Item ( r, i, j ) ) * dy ;
            Real2DArray_Item ( &C, 1, 3 ) = ( Real2DArray_Item ( r, i, j+1 ) + Real2DArray_Item ( r, i, j ) - 2.0e+00 * a ) * ( dy * dy ) ;

            a = ( Real2DArray_Item ( u, i+1, j+1 ) + Real2DArray_Item ( u, i, j ) - Real2DArray_Item ( u, i+1, j ) - Real2DArray_Item ( u, i, j+1 ) ) * ( dx * dx * dy * dy ) -
                ( Real2DArray_Item ( p, i, j+1 ) - Real2DArray_Item ( p, i, j ) ) * ( dx * dy * dy ) - ( Real2DArray_Item ( q, i+1, j ) - Real2DArray_Item ( q, i, j ) ) * ( dx * dx * dy ) + Real2DArray_Item ( r, i, j ) * ( dx * dy ) ;
            b = ( Real2DArray_Item ( p, i+1, j+1 ) + Real2DArray_Item ( p, i, j ) - Real2DArray_Item ( p, i+1, j ) - Real2DArray_Item ( p, i, j+1 ) ) * ( dx * dy * dy ) - ( Real2DArray_Item ( r, i+1, j ) - Real2DArray_Item ( r, i, j ) ) * ( dx * dy ) ;
            c = ( Real2DArray_Item ( q, i+1, j+1 ) + Real2DArray_Item ( q, i, j ) - Real2DArray_Item ( q, i+1, j ) - Real2DArray_Item ( q, i, j+1 ) ) * ( dx * dx * dy ) - ( Real2DArray_Item ( r, i, j+1 ) - Real2DArray_Item ( r, i, j ) ) * ( dx * dy ) ;
            d = ( Real2DArray_Item ( r, i+1, j+1 ) + Real2DArray_Item ( r, i, j ) - Real2DArray_Item ( r, i+1, j ) - Real2DArray_Item ( r, i, j+1 ) ) * ( dx * dy ) ;
            Real2DArray_Item ( &C, 2, 2 ) =    9.0e+00 * a - 3.0e+00 * b - 3.0e+00 * c + d ;
            Real2DArray_Item ( &C, 2, 3 ) = ( -6.0e+00 * a + 2.0e+00 * b + 3.0e+00 * c - d ) * dy ;
            Real2DArray_Item ( &C, 3, 2 ) = ( -6.0e+00 * a + 3.0e+00 * b + 2.0e+00 * c - d ) * dx ;
            Real2DArray_Item ( &C, 3, 3 ) = (  4.0e+00 * a - 2.0e+00 * b - 2.0e+00 * c + d ) * dx * dy ;
/*
printf ( "\nCoefficients for (%d,%d):\n", i, j ) ;
Real2DArray_Print ( &C ) ;
*/
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get the 1-D derivatives.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void BicubicSpline_Get1DDerivatives ( const Real1DArray *x, const Real1DArray *y, Real1DArray *d, BicubicSplineType type, Real1DArray *Ad, Real1DArray *Asd, Real1DArray *qdy, Real1DArray *lll )
{
    Integer     i, n ;
    Real        r ;
    Real1DArray a, b, c ;

    /* . Initialization. */
    n = x->length ;

    /* . Setup. */
    for ( i = 0 ; i < n-1 ; i++ )
    {
        Real1DArray_Item ( Asd, i ) = 1.0e+00 / ( Real1DArray_Item ( x, i+1 ) - Real1DArray_Item ( x, i ) ) ;
        Real1DArray_Item ( qdy, i ) = ( Real1DArray_Item ( y, i+1 ) - Real1DArray_Item ( y, i ) ) * pow ( Real1DArray_Item ( Asd, i ), 2 ) ;
    }
    for ( i = 1 ; i < n-1 ; i++ )
    {
        Real1DArray_Item ( Ad, i ) = 2.0e+00 * ( Real1DArray_Item ( Asd, i-1 ) + Real1DArray_Item ( Asd, i ) ) ;
        Real1DArray_Item (  d, i ) = 3.0e+00 * ( Real1DArray_Item ( qdy, i-1 ) + Real1DArray_Item ( qdy, i ) ) ;
    }

    /* . Branch on the type. */
    switch ( type )
    {
       case BicubicSplineType_Clamped:
           Real1DArray_Item ( d, 1   ) -= Real1DArray_Item ( d, 0   ) * Real1DArray_Item ( Asd, 0   ) ;
           Real1DArray_Item ( d, n-2 ) -= Real1DArray_Item ( d, n-1 ) * Real1DArray_Item ( Asd, n-2 ) ;
           Real1DArray_Slice ( Ad ,  1, n-1, 1, &a, NULL ) ;
           Real1DArray_Slice ( Asd,  1, n-1, 1, &b, NULL ) ;
           Real1DArray_Slice ( d  ,  1, n-1, 1, &c, NULL ) ;
           BicubicSpline_TridiagonalLDLtSolve ( n-2, &a, &b, &c ) ;
           break ;
       case BicubicSplineType_Natural:
           Real1DArray_Item ( Ad, 0   ) = 2.0e+00 * Real1DArray_Item ( Asd, 0   ) ;
           Real1DArray_Item (  d, 0   ) = 3.0e+00 * Real1DArray_Item ( qdy, 0   ) ;
           Real1DArray_Item ( Ad, n-1 ) = 2.0e+00 * Real1DArray_Item ( Asd, n-2 ) ;
           Real1DArray_Item (  d, n-1 ) = 3.0e+00 * Real1DArray_Item ( qdy, n-2 ) ;
           BicubicSpline_TridiagonalLDLtSolve ( n, Ad, Asd, d ) ;
           break ;
       case BicubicSplineType_NotAKnot:
           r = Real1DArray_Item ( Asd, 1   ) / Real1DArray_Item ( Asd, 0   ) ;
           Real1DArray_Item ( Ad, 0   ) = Real1DArray_Item ( Asd, 0   ) / ( 1.0e+00 + r ) ;
           Real1DArray_Item ( d,  0   ) = ( ( 3.0e+00 * r + 2.0e+00 ) * Real1DArray_Item ( qdy, 0   ) + r * Real1DArray_Item ( qdy, 1   ) ) / pow ( ( 1.0e+00 + r ), 2 ) ;
           r = Real1DArray_Item ( Asd, n-3 ) / Real1DArray_Item ( Asd, n-2 ) ;
           Real1DArray_Item ( Ad, n-1 ) = Real1DArray_Item ( Asd, n-2 ) / ( 1.0e+00 + r ) ;
           Real1DArray_Item ( d,  n-1 ) = ( ( 3.0e+00 * r + 2.0e+00 ) * Real1DArray_Item ( qdy, n-2 ) + r * Real1DArray_Item ( qdy, n-3 ) ) / pow ( ( 1.0e+00 + r ), 2 ) ;
           BicubicSpline_TridiagonalLDLtSolve ( n, Ad, Asd, d ) ;
           break ;
       case BicubicSplineType_Periodic:
           Real1DArray_Set ( lll, 0.0e+00 ) ;
           Real1DArray_Item ( Ad,  0   ) = 2.0e+00 * ( Real1DArray_Item ( Asd, 0 ) + Real1DArray_Item ( Asd, n-2 ) ) ;
           Real1DArray_Item ( d,   0   ) = 3.0e+00 * ( Real1DArray_Item ( qdy, 0 ) + Real1DArray_Item ( qdy, n-2 ) ) ;
           Real1DArray_Item ( lll, 0   ) = Real1DArray_Item ( Asd, n-2 ) ;
           Real1DArray_Item ( lll, n-3 ) = Real1DArray_Item ( Asd, n-3 ) ;
           BicubicSpline_TridiagonalLDLtSolvePeriodic ( n-1, Ad, Asd, lll, d ) ;
           Real1DArray_Item ( d, n-1 ) = Real1DArray_Item ( d, 0 ) ;
           break ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Locate the index of a point such that point lies between i and i+1.
! . Points outside the range return -1.
! . This should never happen for periodic abscissa.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Integer BicubicSpline_Locate ( const Real1DArray *abscissa, const Real x, const Boolean isPeriodic, Status *status )
{
    Integer index = -1 ;
    if ( abscissa != NULL )
    {
        auto Integer n ;
        auto Real    a, al, au ;

        /* . Initialization. */
        n  = abscissa->length ;
        a  = x ;
        al = Real1DArray_Item ( abscissa, 0     ) ;
        au = Real1DArray_Item ( abscissa, n - 1 ) ;

        /* . Adjust x for a periodic abscissa. */
        if ( ( isPeriodic ) && ( ( a < al ) || ( a > au ) ) )
        {
            auto Real dx, rF, rI ;
            dx = au - al ;
            rF = modf ( ( a - al ) / dx, &rI );
            if ( rF >= 0.0e+00 ) a = al + rF * dx ;
            else                 a = au + rF * dx ;
            if      ( a < al ) a = al ;
            else if ( a > au ) a = au ;
        }

        /* . Locate the index. */
        if ( ( a >= al ) && ( a <= au ) )
        {
            auto Integer i, l, u ;
            l = 0     ;
            u = n - 1 ;
            while ( ( u - l ) > 1 )
            {
                i = ( u + l ) >> 1 ;
                if ( Real1DArray_Item ( abscissa, i ) >= a ) u = i ;
                else                                         l = i ;
            }
            index = l ;
        }
    }
    return index ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make a spline given (x,y,f) and a type.
!---------------------------------------------------------------------------------------------------------------------------------*/
BicubicSpline *BicubicSpline_MakeFromReal2DArray ( Real1DArray **x, Real1DArray **y, Real2DArray **f, const BicubicSplineType type, Status *status )
{
    BicubicSpline *self = NULL ;
    Boolean        isOK ;
    isOK = ( x != NULL ) && ( (*x) != NULL ) && ( y != NULL ) && ( (*y) != NULL ) && ( f != NULL ) && ( (*f) != NULL ) ;
    if ( isOK )
    {
        auto Integer lengthx, lengthy ;

        /* . Checks. */
        lengthx = (*x)->length ;
        lengthy = (*y)->length ;
        isOK = ( lengthx > 1 ) && ( lengthy > 1 ) && ( (*f)->length0 == lengthx ) && ( (*f)->length1 = lengthy ) ;
        if ( isOK )
        {
            auto Integer i ;
            for ( i = 1 ; i < lengthx ; i++ ) { if ( Real1DArray_Item ( (*x), i ) <=  Real1DArray_Item ( (*x), i-1 ) ) { isOK = False ; break ; } }
            for ( i = 1 ; i < lengthy ; i++ ) { if ( Real1DArray_Item ( (*y), i ) <=  Real1DArray_Item ( (*y), i-1 ) ) { isOK = False ; break ; } }
        }

        /* . OK so far. */
        if ( isOK )
        {
            /* . Allocate space. */
            self = BicubicSpline_Allocate ( lengthx, lengthy, False, False, False, status ) ;
            if ( self != NULL )
            {
                auto Status localstatus ;
                Status_Set ( &localstatus, Status_Continue ) ;
                /* . Finish construction. */
                self->type = type ;
                self->x = (*x) ; (*x) = NULL ;
                self->y = (*y) ; (*y) = NULL ;
                self->f = (*f) ; (*f) = NULL ;
                BicubicSpline_Setup ( self, &localstatus ) ;
                if ( ! Status_OK ( &localstatus ) )
                {
                    BicubicSpline_Deallocate ( &self ) ;
                    Status_Set ( status, localstatus ) ;
                }
            }
            else Status_Set ( status, Status_OutOfMemory ) ;
        }
    }
    if ( ! isOK ) Status_Set ( status, Status_InvalidArgument ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Setup the bicubic spline.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void BicubicSpline_Setup ( BicubicSpline *self, Status *status )
{
    Integer      n ;
    Real1DArray *Ad, *Asd, *d, *ll = NULL, *qdu, t, u ;
    Real2DArray *p, *q, *r ;
    Status       localstatus ;

    /* . Initialization . */
    n = Maximum ( self->lengthx, self->lengthy ) ;
    Status_Set ( &localstatus, Status_Continue ) ;

    /* . Allocate space. */
    Ad  = Real1DArray_Allocate ( n, &localstatus ) ;
    Asd = Real1DArray_Allocate ( n, &localstatus ) ;
    d   = Real1DArray_Allocate ( self->lengthy, &localstatus ) ;
    if ( self->type == BicubicSplineType_Periodic ) ll  = Real1DArray_Allocate ( n, &localstatus ) ;
    qdu = Real1DArray_Allocate ( n, &localstatus ) ;
    p   = Real2DArray_Allocate ( self->lengthx, self->lengthy, &localstatus ) ;
    q   = Real2DArray_Allocate ( self->lengthx, self->lengthy, &localstatus ) ;
    r   = Real2DArray_Allocate ( self->lengthx, self->lengthy, &localstatus ) ;
    if ( Status_OK ( &localstatus ) )
    {
        auto Integer i ;

        /* . du/dx. */
        for ( i = 0 ; i < self->lengthy ; i++ )
        {
            Real2DArray_ColumnSlice ( p      , i, &t, &localstatus ) ;
            Real2DArray_ColumnSlice ( self->f, i, &u, &localstatus ) ;
            BicubicSpline_Get1DDerivatives ( self->x, &u, &t, self->type, Ad, Asd, qdu, ll ) ;
        }

        /* . du/dy. */
        for ( i = 0 ; i < self->lengthx ; i++ )
        {
            Real2DArray_RowSlice ( q      , i, &t, &localstatus ) ;
            Real2DArray_RowSlice ( self->f, i, &u, &localstatus ) ;
            BicubicSpline_Get1DDerivatives ( self->y, &u, d, self->type, Ad, Asd, qdu, ll ) ;
            Real1DArray_CopyTo ( d, &t, &localstatus ) ;
        }

        /* . d2u/dxdy. */
        Real2DArray_ColumnSlice ( q, 0, &u, &localstatus ) ;
        Real2DArray_ColumnSlice ( r, 0, &t, &localstatus ) ;
        BicubicSpline_Get1DDerivatives ( self->x, &u, &t, self->type, Ad, Asd, qdu, ll ) ;
        Real2DArray_ColumnSlice ( q, self->lengthy-1, &u, &localstatus ) ;
        Real2DArray_ColumnSlice ( r, self->lengthy-1, &t, &localstatus ) ;
        BicubicSpline_Get1DDerivatives ( self->x, &u, &t, self->type, Ad, Asd, qdu, ll ) ;
        for ( i = 0 ; i < self->lengthx ; i++ )
        {
            Real2DArray_RowSlice ( p, i, &u, &localstatus ) ;
            Real1DArray_Item ( d, 0               ) = Real2DArray_Item ( r, i, 0               ) ;
            Real1DArray_Item ( d, self->lengthy-1 ) = Real2DArray_Item ( r, i, self->lengthy-1 ) ;
            BicubicSpline_Get1DDerivatives ( self->y, &u, d, BicubicSplineType_Clamped, Ad, Asd, qdu, ll ) ;
            Real1DArray_Slice   ( d,            1, self->lengthy-1, 1, &u, &localstatus ) ;
            Real2DArray_1DSlice ( r, i, i+1, 1, 1, self->lengthy-1, 1, &t, &localstatus ) ;
            Real1DArray_CopyTo ( &u, &t, &localstatus ) ;
        }

        /* . Coefficient table. */
        BicubicSpline_EvaluateCoefficientTable ( self, p, q, r ) ;
/*
printf ( "\nAbscissae for %u:\n", self->type ) ;
Real1DArray_Print ( self->x ) ;
Real1DArray_Print ( self->y ) ;
printf ( "\nFunction Values:\n" ) ;
Real2DArray_Print ( self->f ) ;
*/
    }
    else Status_Set ( status, localstatus ) ;

    /* . Finish up. */
    Real1DArray_Deallocate ( &Ad  ) ;
    Real1DArray_Deallocate ( &Asd ) ;
    Real1DArray_Deallocate ( &d   ) ;
    Real1DArray_Deallocate ( &ll  ) ;
    Real1DArray_Deallocate ( &qdu ) ;
    Real2DArray_Deallocate ( &p   ) ;
    Real2DArray_Deallocate ( &q   ) ;
    Real2DArray_Deallocate ( &r   ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Solution of a tridiagonal system.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void BicubicSpline_TridiagonalLDLtSolve ( const Integer n, Real1DArray *d, Real1DArray *l, Real1DArray *b )
{
    auto Integer i ;
    auto Real temp ;
    for ( i = 1 ; i < n ; i++ )
    {
        temp = Real1DArray_Item ( l, i-1 ) ;
        Real1DArray_Item ( l, i-1 ) /= Real1DArray_Item ( d, i-1 ) ;
        Real1DArray_Item ( d, i   ) -= temp * Real1DArray_Item ( l, i-1 ) ;
        Real1DArray_Item ( b, i   ) -= Real1DArray_Item ( l, i-1 ) * Real1DArray_Item ( b, i-1 ) ;
    }
    Real1DArray_Item ( b, n-1 ) /= Real1DArray_Item ( d, n-1 ) ;
    for ( i = n-2 ; i >= 0 ; i-- ) Real1DArray_Item ( b, i ) = ( Real1DArray_Item ( b, i ) / Real1DArray_Item ( d, i ) - Real1DArray_Item ( l, i ) * Real1DArray_Item ( b, i+1 ) ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Solution of a periodic tridiagonal system.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void BicubicSpline_TridiagonalLDLtSolvePeriodic ( const Integer n, Real1DArray *d, Real1DArray *lsd, Real1DArray *lll, Real1DArray *b )
{
    auto Integer i ;
    auto Real temp1, temp2 ;
    for ( i = 0 ; i < n-2 ; i++ )
    {
       temp1 = Real1DArray_Item ( lsd, i ) ;
       temp2 = Real1DArray_Item ( lll, i ) ;
       Real1DArray_Item ( lsd, i   ) /= Real1DArray_Item ( d, i ) ;
       Real1DArray_Item ( lll, i   ) /= Real1DArray_Item ( d, i ) ;
       Real1DArray_Item (   d, i+1 ) -= Real1DArray_Item ( lsd, i ) * temp1 ;
       Real1DArray_Item ( lll, i+1 ) -= Real1DArray_Item ( lll, i ) * temp1 ;
       Real1DArray_Item (   d, n-1 ) -= Real1DArray_Item ( lll, i ) * temp2 ;
    }
    temp2 = Real1DArray_Item ( lll, n-2 ) ;
    Real1DArray_Item ( lll, n-2 ) /= Real1DArray_Item (   d, n-2 ) ;
    Real1DArray_Item (   d, n-1 ) -= Real1DArray_Item ( lll, n-2 ) * temp2 ;
    for ( i = 1 ; i < n-1 ; i++ ) Real1DArray_Item ( b, i   ) -= Real1DArray_Item ( lsd, i-1 ) * Real1DArray_Item ( b, i-1 ) ;
    for ( i = 0 ; i < n-1 ; i++ ) Real1DArray_Item ( b, n-1 ) -= Real1DArray_Item ( lll, i   ) * Real1DArray_Item ( b, i   ) ;
    for ( i = 0 ; i < n   ; i++ ) Real1DArray_Item ( b, i   ) /= Real1DArray_Item (   d, i   ) ;
    Real1DArray_Item ( b, n-2 ) -= Real1DArray_Item ( lll, n-2 ) * Real1DArray_Item ( b, n-1 ) ;
    for ( i = n-3 ; i >= 0 ; i-- ) Real1DArray_Item ( b, i ) -= ( Real1DArray_Item ( lsd, i ) * Real1DArray_Item ( b, i+1 ) + Real1DArray_Item ( lll, i ) * Real1DArray_Item ( b, n-1 ) ) ;
}

