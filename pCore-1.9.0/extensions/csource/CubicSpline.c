/*------------------------------------------------------------------------------
! . File      : CubicSpline.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . Cubic splines.
! . Follows the symmetrical formulation in Numerical Recipes.
!=================================================================================================================================*/

# include <stdlib.h>

# include "CubicSpline.h"
# include "f2clapack.h"
# include "math.h"
# include "Memory.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
Status CubicSpline_Allocate ( CubicSpline **self, const Integer length, const Boolean QX, const Boolean QY )
{
    Status status = Status_Null ;
    if ( self != NULL )
    {
        (*self) = NULL ;
        if ( length >= 2 )
        {
            MEMORY_ALLOCATE ( (*self), CubicSpline ) ;
            if ( (*self) != NULL )
            {
                (*self)->length = length ;
                (*self)->h = NULL ;
                (*self)->x = NULL ;
                (*self)->y = NULL ;
                /* . Array allocation. */
                status = Status_Success ;
                (*self)->h = Real1DArray_Allocate ( length, &status ) ;
                if ( QX && ( status == Status_Success ) ) (*self)->x = Real1DArray_Allocate ( length, &status ) ;
                if ( QY && ( status == Status_Success ) ) (*self)->y = Real1DArray_Allocate ( length, &status ) ;
                if ( status != Status_Success ) CubicSpline_Deallocate ( self ) ;
            }
            if ( (*self) == NULL ) status = Status_OutOfMemory ;
            else                   status = Status_Success     ;
        }
        else status = Status_DimensionError ;
    }
    return status ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
Status CubicSpline_Clone ( const CubicSpline *self, CubicSpline **clone )
{
    Status status = Status_Null ;
    if ( ( self != NULL ) && ( clone != NULL ) )
    {
        status = CubicSpline_Allocate ( clone, self->length, True, True ) ;
        if ( status == Status_Success )
        {
            Real1DArray_CopyTo ( self->h, (*clone)->h, NULL ) ;
            Real1DArray_CopyTo ( self->x, (*clone)->x, NULL ) ;
            Real1DArray_CopyTo ( self->y, (*clone)->y, NULL ) ;
        }
    }
    return status ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void CubicSpline_Deallocate ( CubicSpline **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) )
    {
        Real1DArray_Deallocate ( &((*self)->h) ) ;
        Real1DArray_Deallocate ( &((*self)->x) ) ;
        Real1DArray_Deallocate ( &((*self)->y) ) ;
        MEMORY_DEALLOCATE ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Evaluation.
!---------------------------------------------------------------------------------------------------------------------------------*/
Status CubicSpline_Evaluate ( const CubicSpline *self, const Real x, Real *f, Real *g, Real *h )
{
    Status status = Status_Null ;
    if ( f != NULL ) (*f) = 0.0e+00 ;
    if ( g != NULL ) (*g) = 0.0e+00 ;
    if ( h != NULL ) (*h) = 0.0e+00 ;
    if ( self != NULL )
    {
        auto Integer n = self->length ;
        auto Real1DArray *abscissa = self->x ;

        /* . Check ranges.*/
        if ( ( x < Real1DArray_Item ( abscissa, 0 ) ) || ( x > Real1DArray_Item ( abscissa, n-1 ) ) ) status = Status_InvalidArgument ;
        /* . OK. */
        else
        {
            auto Integer i, l, u ;
            auto Real    hl, hu, d, s, t, yl, yu ;

            /* . Locate the index. */
            l = 0     ;
            u = n - 1 ;
            while ( ( u - l ) > 1 )
            {
                i = ( u + l ) >> 1 ;
                if ( Real1DArray_Item ( abscissa, i ) > x ) u = i ;
                else                                        l = i ;
            }

            /* . Factors. */
            d  = ( Real1DArray_Item ( abscissa, u ) - Real1DArray_Item ( abscissa, l ) ) ;
            s  = ( x - Real1DArray_Item ( abscissa, l ) ) / d ;
            t  = ( Real1DArray_Item ( abscissa, u ) - x ) / d ;
            hl = Real1DArray_Item ( self->h, l ) * d / 6.0e+00 ;
            hu = Real1DArray_Item ( self->h, u ) * d / 6.0e+00 ;
            yl = Real1DArray_Item ( self->y, l ) ;
            yu = Real1DArray_Item ( self->y, u ) ;

            /* . Calculate the function and its derivatives. */
            if ( f != NULL ) (*f) = t * yl + s * yu + d * ( t * ( t * t - 1.0e+00 ) * hl + s * ( s * s - 1.0e+00 ) * hu ) ;
            if ( g != NULL ) (*g) = ( yu - yl ) / d + ( - ( 3.0e+00 * t * t - 1.0e+00 ) * hl + ( 3.0e+00 * s * s - 1.0e+00 ) * hu ) ;
            if ( h != NULL ) (*h) = 6.0e+00 * ( t * hl + s * hu ) / d ;
        }
    }
    return status ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Evaluation of quantities required to evaluate the spline - there is no checking!
!---------------------------------------------------------------------------------------------------------------------------------*/
void CubicSpline_EvaluateLUDST ( const CubicSpline *self, const Real x, Integer *l, Integer *u, Real *d, Real *s, Real *t )
{
    auto Integer      i, l0, u0 ;
    auto Real1DArray *abscissa = self->x ;

    /* . Locate the index. */
    l0 = 0 ;
    u0 = self->length - 1 ;
    while ( ( u0 - l0 ) > 1 )
    {
        i = ( u0 + l0 ) >> 1 ;
        if ( Real1DArray_Item ( abscissa, i ) > x ) u0 = i ;
        else                                        l0 = i ;
    }

    /* . Save data. */
    (*l) = l0 ;
    (*u) = u0 ;
    (*d) = ( Real1DArray_Item ( abscissa, u0 ) - Real1DArray_Item ( abscissa, l0 ) ) ;
    (*s) = ( x - Real1DArray_Item ( abscissa, l0 ) ) / (*d) ;
    (*t) = ( Real1DArray_Item ( abscissa, u0 ) - x ) / (*d) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Find the positions of any extrema (maxima and minima only).
!---------------------------------------------------------------------------------------------------------------------------------*/
Status CubicSpline_FindExtrema ( const CubicSpline *self, Real1DArray **maxima, Real1DArray **minima )
{
    Status status = Status_Null ;
    /* . Check for non-NULL arguments. */
    if ( ( ( maxima != NULL ) && ( (*maxima) != NULL ) ) || ( ( minima != NULL ) && ( (*minima) != NULL ) ) ) status = Status_InvalidArgument ;
    /* . Normal checks. */
    else if ( ( self != NULL ) && ( ( maxima != NULL ) || ( minima != NULL ) ) )
    {
        auto Integer n = self->length ;
        auto Real1DArray *tmaxima = NULL, *tminima = NULL ;

        /* . Allocate temporary space. */
        status = Status_Success ;
        if (   maxima != NULL )                                   tmaxima = Real1DArray_Allocate ( n + 1, &status ) ;
        if ( ( minima != NULL ) && ( status == Status_Success ) ) tminima = Real1DArray_Allocate ( n + 1, &status ) ;

        if ( status == Status_Success )
        {
            auto Integer i, itry, nmaxima = 0, nminima = 0, n = self->length, ntry = 0 ;
            auto Real    a, b, c, d, factor, h, hl, hu, s, sign, x, xl, xu, yl, yu ;
            auto Real1DArray *abscissa = self->x ;

            /* . Each interval will have (in principle) a maximum of one maximum and one minimum. */
            for ( i = 0 ; i < n - 1 ; i++ )
            {
                /* . The interval range. */
                xl = Real1DArray_Item ( abscissa, i     ) ;
                xu = Real1DArray_Item ( abscissa, i + 1 ) ;
                d  = xu - xl ;

                /* . Other factors. */
                hl = Real1DArray_Item ( self->h, i     ) * d / 6.0e+00 ;
                hu = Real1DArray_Item ( self->h, i + 1 ) * d / 6.0e+00 ;
                yl = Real1DArray_Item ( self->y, i     ) ;
                yu = Real1DArray_Item ( self->y, i + 1 ) ;

                /* . Factors for a quadratic equation in terms of s. */
                a = 3.0e+00 * ( hu - hl ) ;
                b = 6.0e+00 * hl ;
                c = ( yu - yl ) / d - 2.0e+00 * hl - hu ;
                factor = b * b - 4.0e+00 * a * c ;

                /* . Solve the first derivative quadratic. */
                if ( factor >= 0.0e+00 )
                {
                    if ( factor == 0.0e+00 ) ntry = 1 ;
                    else { factor = sqrt ( factor ) ; ntry = 2 ; }
                    for ( itry = 0, sign = -1.0e+00 ; itry < ntry ; itry++ )
                    {
                        s = ( - b + sign * factor ) / ( 2.0e+00 * a ) ;
                        /* . Only check for lower boundary extrema except for the last interval when both are checked. */
                        if ( ( s >= 0.0e+00 ) && ( ( ( s < 1.0e+00 ) && ( i < n - 2 ) ) || ( ( s <= 1.0e+00 ) && ( i == n - 2 ) ) ) )
                        {
                            x = d * s + xl ;
                            /* . Apply the second derivative test. */
                            /* . The more general extremum test is unnecessary for a cubic function because if the second derivative is zero
                            ! .  the point is an inflection point (when the third derivative is non-zero) or the function is a constant
                            ! . (as all derivatives are zero). */
                            CubicSpline_Evaluate ( self, x, NULL, NULL, &h ) ;
                                 if ( ( h > 0.0e+00 ) && ( tminima != NULL ) ) { Real1DArray_Item ( tminima, nminima ) = x ; nminima ++ ; }
                            else if ( ( h < 0.0e+00 ) && ( tmaxima != NULL ) ) { Real1DArray_Item ( tmaxima, nmaxima ) = x ; nmaxima ++ ; }
                        }
                        sign *= -1.0e00 ;
                    }
                }
            }

            /* . Save the extrema. */
            if ( nmaxima > 0 )
            {
                Real1DArray_Resize ( tmaxima, nmaxima, NULL, &status ) ;
                (*maxima) = tmaxima ;
            }
            if ( ( nminima > 0 ) && ( status = Status_Success ) )
            {
                Real1DArray_Resize ( tminima, nminima, NULL, &status ) ;
                (*minima) = tminima ;
            }
        }
    }
    return status ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the integral of the spline in the range [a,b].
!---------------------------------------------------------------------------------------------------------------------------------*/
Real CubicSpline_Integrate ( const CubicSpline *self, const Real a, const Real b, Status *status )
{
    Real integral = 0.0e+00 ;
    if ( ( self != NULL ) && ( a != b ) )
    {
        /* . Check the range of the integral.*/
        if ( ( a < Real1DArray_Item ( self->x, 0 ) ) || ( a > Real1DArray_Item ( self->x, self->length-1 ) ) ||
             ( b < Real1DArray_Item ( self->x, 0 ) ) || ( b > Real1DArray_Item ( self->x, self->length-1 ) ) ||
             ( a > b ) ) Status_Set ( status, Status_InvalidArgument ) ;
        /* . OK. */
        else
        {
            Integer l, lA, lB, u ;
            Real    d, local, sA, sB, sF2, sF4, tA, tB, tF2, tF4 ;

            /* . Get data for a and b. */
            CubicSpline_EvaluateLUDST ( self, a, &lA, &u, &d, &sA, &tA ) ;
            CubicSpline_EvaluateLUDST ( self, b, &lB, &u, &d, &sB, &tB ) ;

            /* . Loop over intervals. */
            for ( l = lA ; l <= lB ; l++ )
            {
                u   = l + 1 ;
                d   = Real1DArray_Item ( self->x, u ) - Real1DArray_Item ( self->x, l ) ;
                sF4 = Real1DArray_Item ( self->h, u ) * d * d / 6.0e+00 ;
                sF2 = 2.0e+00 * ( Real1DArray_Item ( self->y, u ) - sF4 ) ;
                tF4 = Real1DArray_Item ( self->h, l ) * d * d / 6.0e+00 ;
                tF2 = 2.0e+00 * ( Real1DArray_Item ( self->y, l ) - tF4 ) ;
                if ( l == lA ) local  = tA * tA * ( tF2 + tA * tA * tF4 ) - sA * sA * ( sF2 + sA * sA * sF4 ) ;
                else           local  = tF2 + tF4 ;
                if ( l == lB ) local += sB * sB * ( sF2 + sB * sB * sF4 ) - tB * tB * ( tF2 + tB * tB * tF4 ) ;
                else           local += sF2 + sF4 ;
                integral += 0.25e+00 * d * local ;
            }
        }
    }
    return integral ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the integral of the full spline.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real CubicSpline_IntegrateFull ( const CubicSpline *self, Status *status )
{
    Real integral = 0.0e+00 ;
    if ( self != NULL )
    {
        integral = CubicSpline_Integrate ( self, Real1DArray_Item ( self->x, 0 ), Real1DArray_Item ( self->x, self->length-1 ), status ) ;
    }
    return integral ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make a spline given (x,y) pairs and boundary conditions.
! . The X and Y arrays are taken by the spline.
! . X may be absent, in which case integer values are assumed.
! . Periodic splines requires solving a cyclic tridiagonal system (to be added if necessary).
! . Checks should probably be added to ensure the data in X is OK - i.e. ordered and intervals are non-zero.
!---------------------------------------------------------------------------------------------------------------------------------*/
Status CubicSpline_MakeFromReal1DArrays ( CubicSpline **self, Real1DArray **x, Real1DArray **y, const Integer lowerderivative, const Real lowervalue, const Integer upperderivative, const Real uppervalue )
{
    Status status = Status_Null ;
    if ( self != NULL )
    {
        (*self) = NULL ;
        if ( ( y != NULL ) && ( (*y) != NULL ) )
        {
            auto Boolean QNEEDX ;
            auto Integer n      ;

            /* . Set some quantities. */
            n      = (*y)->length     ;
            QNEEDX = ( x == NULL ) || ( (*x) == NULL ) ;

            /* . Basic checks. */
            if ( ! QNEEDX && ( (*x)->length != n ) ) status = Status_DimensionError ;
            else if ( ( lowerderivative < 1 ) || ( lowerderivative > 2 ) || ( upperderivative < 1 ) || ( upperderivative > 2 ) ) status = Status_InvalidArgument ;
            /* . Everything looks OK so far. */
            else
            {
                auto Real1DArray *diagonal = NULL, *subdiagonal = NULL, *superdiagonal = NULL ;

                /* . Allocate space. */
                status = CubicSpline_Allocate ( self, n, QNEEDX, False ) ;
                if ( status == Status_Success ) diagonal      = Real1DArray_Allocate ( n    , &status ) ;
                if ( status == Status_Success ) subdiagonal   = Real1DArray_Allocate ( n - 1, &status ) ;
                if ( status == Status_Success ) superdiagonal = Real1DArray_Allocate ( n - 1, &status ) ;
                if ( status == Status_Success )
                {
                    auto Integer i ;
                    auto Real    dl, du ;
                    auto Real1DArray *rhs ;

                    /* . Alias h to rhs and initialize. */
                    rhs = (*self)->h ;
                    Real1DArray_Set ( rhs, 0.0e+00 ) ;

                    /* . Ensure the remaining data is set. */
                    if ( QNEEDX ) for ( i = 0 ; i < n ; i++ ) Real1DArray_Item ( (*self)->x, i ) = ( Real ) i ;
                    else { (*self)->x = (*x) ; (*x) = NULL ; }
                    (*self)->y = (*y) ; (*y) = NULL ;

                    /* . Reverse the array if x is decreasing. */
                    if ( Real1DArray_Item ( (*self)->x, 0 ) > Real1DArray_Item ( (*self)->x, n-1 ) )
                    {
                        Real1DArray_Reverse ( (*self)->x ) ;
                        Real1DArray_Reverse ( (*self)->y ) ;
                    }

                    /* . Set up the tridiagonal system. */
                    /* . Lower boundary. */
                    if ( lowerderivative == 1 )
                    {
                        dl = ( Real1DArray_Item ( (*self)->x, 1 ) - Real1DArray_Item ( (*self)->x, 0 ) ) ;
                        Real1DArray_Item ( diagonal,      0 ) = dl / 3.0e+00 ;
                        Real1DArray_Item ( superdiagonal, 0 ) = dl / 6.0e+00 ;
                        Real1DArray_Item ( rhs,           0 ) = ( Real1DArray_Item ( (*self)->y, 1 ) - Real1DArray_Item ( (*self)->y, 0 ) ) / dl - lowervalue ;
                    }
                    else if ( lowerderivative == 2 )
                    {
                        Real1DArray_Item ( diagonal,      0 ) = 1.0e+00    ;
                        Real1DArray_Item ( superdiagonal, 0 ) = 0.0e+00    ;
                        Real1DArray_Item ( rhs,           0 ) = lowervalue ;
                    }
                    /* . Interior conditions. */
                    for ( i = 1 ; i < n - 1 ; i++ )
                    {
                        dl = ( Real1DArray_Item ( (*self)->x, i   ) - Real1DArray_Item ( (*self)->x, i-1 ) ) ;
                        du = ( Real1DArray_Item ( (*self)->x, i+1 ) - Real1DArray_Item ( (*self)->x, i   ) ) ;
                        Real1DArray_Item (   subdiagonal, i-1 ) =   dl        / 6.0e+00 ;
                        Real1DArray_Item (      diagonal, i   ) = ( dl + du ) / 3.0e+00 ;
                        Real1DArray_Item ( superdiagonal, i   ) =        du   / 6.0e+00 ;
                        Real1DArray_Item ( rhs, i ) = ( Real1DArray_Item ( (*self)->y, i+1 ) - Real1DArray_Item ( (*self)->y, i ) ) / du +
                                                      ( Real1DArray_Item ( (*self)->y, i-1 ) - Real1DArray_Item ( (*self)->y, i ) ) / dl ;
                    }
                    /* . Upper boundary. */
                    if ( upperderivative == 1 )
                    {
                        du = ( Real1DArray_Item ( (*self)->x, n-1 ) - Real1DArray_Item ( (*self)->x, n-2 ) ) ;
                        Real1DArray_Item ( subdiagonal, n-2 ) = du / 6.0e+00 ;
                        Real1DArray_Item ( diagonal,    n-1 ) = du / 3.0e+00 ;
                        Real1DArray_Item ( rhs,         n-1 ) = uppervalue - ( Real1DArray_Item ( (*self)->y, n-1 ) - Real1DArray_Item ( (*self)->y, n-2 ) ) / du ;
                    }
                    else if ( upperderivative == 2 )
                    {
                        Real1DArray_Item ( subdiagonal, n-2 ) = 0.0e+00    ;
                        Real1DArray_Item ( diagonal,    n-1 ) = 1.0e+00    ;
                        Real1DArray_Item ( rhs,         n-1 ) = uppervalue ;
                    }

                    /* . Solve the tridiagonal system. */
                    {
                        auto integer info   = 0 ;
                        auto integer length = n ;
                        auto integer nrhs   = 1 ;
                        dgtsv_ ( &length, &nrhs, Real1DArray_Data ( subdiagonal ), Real1DArray_Data ( diagonal ), Real1DArray_Data ( superdiagonal ), Real1DArray_Data ( rhs ), &length, &info ) ;
                        /* . Set the status. */
                        if ( info == 0 ) status = Status_Success               ;
                        else             status = Status_LinearEquationFailure ;
                    }
                }

                /* . Deallocate space. */
                Real1DArray_Deallocate ( &diagonal      ) ;
                Real1DArray_Deallocate ( &subdiagonal   ) ;
                Real1DArray_Deallocate ( &superdiagonal ) ;
            }
            if ( status != Status_Success ) CubicSpline_Deallocate ( self ) ;
        }
    }
    return status ;
}
