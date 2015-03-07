/*------------------------------------------------------------------------------
! . File      : Correlation.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . Functions for calculating correlation functions.
!=================================================================================================================================*/

/* . For cross-correlation the function calculated is X(t)Y(0) so need a second call to get X(0)Y(t). */

# include <math.h>
# include <stdio.h>

# include "Correlation.h"
# include "Macros.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Parameters.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define DEFAULTSMALL 1.0e-10

/*----------------------------------------------------------------------------------------------------------------------------------
! . Private procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void Correlation_DotProductDirect      ( const Real2DArray *x, const Real2DArray *y, Real1DArray *c, Status *status ) ;
static void Correlation_DotProductNormalize   ( const Real2DArray *x, const Real2DArray *y, const Real small, Real1DArray *c, Status *status ) ;
static void Correlation_DotProductRemoveMeans ( Real2DArray *x, Status *status ) ;

static void Correlation_SimpleDirect          ( const Real1DArray *x, const Real1DArray *y, Real1DArray *c ) ;
static void Correlation_SimpleNormalize       ( const Real1DArray *x, const Real1DArray *y, const Real small, Real1DArray *c ) ;
static void Correlation_SimpleRemoveMean      ( Real1DArray *x ) ;

/*==================================================================================================================================
! . Public procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Dot-product auto or cross-correlation function.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real1DArray *Correlation_MakeDotProduct ( Real2DArray *x, Real2DArray *y, Real1DArray *c, const Boolean useFFT, const Boolean normalize, const Boolean removeMean, const Integer tCorrelation, const Real *tolerance, Status *status )
{
    Real1DArray *f = NULL ;
    if ( x != NULL )
    {
        /* . Get the correlation array. */
        /* . The presence of an array overrides all other options. */
        if ( c == NULL ) f = Real1DArray_Allocate ( tCorrelation, status ) ;
        else             f = c ;
        if ( f != NULL )
        {
            /* . Remove means. */
            if ( removeMean )
            {
                Correlation_DotProductRemoveMeans ( x, status ) ;
                if ( y != x ) Correlation_DotProductRemoveMeans ( y, status ) ;
            }

            /* . Calculate the correlation function. */
            Correlation_DotProductDirect ( x, y, f, status ) ;
/*
            if ( useFFT ) Correlation_DotProductFFT    ( x, y, f, status ) ;
            else          Correlation_DotProductDirect ( x, y, f, status ) ;
*/

            /* . Normalize. */
            if ( normalize )
            {
                auto Real small ;
                if ( tolerance == NULL ) small = DEFAULTSMALL ;
                else                     small = fabs ( *tolerance ) ;
                Correlation_DotProductNormalize ( x, y, small, f, status ) ;
            }
        }
    }
    return f ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Simple auto or cross-correlation function.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real1DArray *Correlation_MakeSimple ( Real1DArray *x, Real1DArray *y, Real1DArray *c, const Boolean useFFT, const Boolean normalize, const Boolean removeMean, const Integer tCorrelation, const Real *tolerance, Status *status )
{
    Real1DArray *f = NULL ;
    if ( x != NULL )
    {
        /* . Get the correlation array. */
        /* . The presence of an array overrides all other options. */
        if ( c == NULL ) f = Real1DArray_Allocate ( tCorrelation, status ) ;
        else             f = c ;
        if ( f != NULL )
        {
            /* . Remove means. */
            if ( removeMean )
            {
                Correlation_SimpleRemoveMean ( x ) ;
                if ( y != x ) Correlation_SimpleRemoveMean ( y ) ;
            }

            /* . Calculate the correlation function. */
            Correlation_SimpleDirect ( x, y, f ) ;
/*
            if ( useFFT ) Correlation_SimpleFFT    ( x, y, f, status ) ;
            else          Correlation_SimpleDirect ( x, y, f ) ;
*/

            /* . Normalize. */
            if ( normalize )
            {
                auto Real small ;
                if ( tolerance == NULL ) small = DEFAULTSMALL ;
                else                     small = fabs ( *tolerance ) ;
                Correlation_SimpleNormalize ( x, y, small, f ) ;
            }
        }
    }
    return f ;
}

/*==================================================================================================================================
! . Private procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Dot-product auto or cross-correlation function via a direct method.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void Correlation_DotProductDirect ( const Real2DArray *x, const Real2DArray *y, Real1DArray *c, Status *status )
{
    if ( ( c != NULL ) && ( x != NULL ) )
    {
        auto       Integer      i, t, tCorrelation, tMaximum, tRun ;
        auto       Real         sum   ;
        auto       Real1DArray  arow, brow ;
        auto const Real2DArray *a, *b ;

        /* . Assign array aliases. */
        a = x ;
        if ( y == NULL ) b = x ;
        else             b = y ;

        /* . Get problem size. */
        tRun         = Minimum ( Real2DArray_Length ( a, 0 ), Real2DArray_Length ( b, 0 ) ) ;
        tCorrelation = Minimum ( Real1DArray_Length ( c ) - 1, tRun - 1 ) ;

        /* . Initialization. */
        Real1DArray_Set ( c, 0.0e+00 ) ;

        /* . Loop over the elements in the correlation function. */
        for ( t = 0 ; t <= tCorrelation ; t++ )
        {
            tMaximum = tRun - t ;
            sum      = 0.0e+00  ;
            for ( i = 0 ; i < tMaximum ; i++ )
            {
                Real2DArray_RowSlice ( a, t+i, &arow, status ) ;
                Real2DArray_RowSlice ( b,   i, &brow, status ) ;
                sum += Real1DArray_Dot ( &arow, &brow, status ) ;
            }
            Real1DArray_Item ( c, t ) = sum / ( ( Real ) tMaximum ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Normalize a dot-product correlation function.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void Correlation_DotProductNormalize ( const Real2DArray *x, const Real2DArray *y, const Real small, Real1DArray *c, Status *status )
{
    if ( ( x != NULL ) && ( c != NULL ) )
    {
        auto Real scale ;
        /* . Auto-correlation. */
        if ( y == NULL ) scale = Real1DArray_Item ( c, 0 ) ;
        /* . Cross-correlation. */
        else
        {
            auto Integer     i ;
            auto Real        x2, y2 ;
            auto Real1DArray column ;
            x2 = 0.0e+00 ;
            y2 = 0.0e+00 ;
            for ( i = 0 ; i < Real2DArray_Length ( x, 1 ) ; i++ )
            {
                Real2DArray_ColumnSlice ( x, i, &column, status ) ;
                x2 += Real1DArray_Dot ( &column, &column, status ) ;
            }
            for ( i = 0 ; i < Real2DArray_Length ( y, 1 ); i++ )
            {
                Real2DArray_ColumnSlice ( y, i, &column, status ) ;
                y2 += Real1DArray_Dot ( &column, &column, status ) ;
            }

            scale = sqrt ( fabs ( x2 * y2 ) ) ;
        }
        /* . Scaling. */
        if ( fabs ( scale ) > small ) Real1DArray_Scale ( c, 1.0 / scale ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Remove the means from a 2-D array of data.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void Correlation_DotProductRemoveMeans ( Real2DArray *x, Status *status )
{
    auto Integer n ;
    n = Real2DArray_Length ( x, 0 ) ;
    if ( ( x != NULL ) && ( n > 0 ) )
    {
        auto Integer     i      ;
        auto Real        mean   ;
        auto Real1DArray column ;
        for ( i = 0 ; i < Real2DArray_Length ( x, 1 ) ; i++ )
        {
            Real2DArray_ColumnSlice ( x, i, &column, status ) ;
            mean = Real1DArray_Sum ( &column ) / ( ( Real ) n ) ;
            Real1DArray_AddScalar ( &column, -mean ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Simple auto or cross-correlation function via a direct method.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void Correlation_SimpleDirect ( const Real1DArray *x, const Real1DArray *y, Real1DArray *c )
{
    if ( ( c != NULL ) && ( x != NULL ) )
    {
        auto       Integer      i, t, tCorrelation, tMaximum, tRun ;
        auto       Real         sum   ;
        auto const Real1DArray *a, *b ;

        /* . Assign array aliases. */
        a = x ;
        if ( y == NULL ) b = x ;
        else             b = y ;

        /* . Get problem size. */
        tRun         = Minimum ( Real1DArray_Length ( a ), Real1DArray_Length ( b ) ) ;
        tCorrelation = Minimum ( Real1DArray_Length ( c ) - 1, tRun - 1 ) ;

        /* . Initialization. */
        Real1DArray_Set ( c, 0.0e+00 ) ;

        /* . Loop over the elements in the correlation function. */
        for ( t = 0 ; t <= tCorrelation ; t++ )
        {
            tMaximum = tRun - t ;
            sum      = 0.0e+00  ;
            for ( i = 0 ; i < tMaximum ; i++ ) sum += ( Real1DArray_Item ( a, t+i ) * Real1DArray_Item ( b, i ) ) ;
            Real1DArray_Item ( c, t ) = sum / ( ( Real ) tMaximum ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Normalize a simple correlation function.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void Correlation_SimpleNormalize ( const Real1DArray *x, const Real1DArray *y, const Real small, Real1DArray *c )
{
    if ( ( x != NULL ) && ( c != NULL ) )
    {
        auto Real scale ;
        /* . Auto-correlation. */
        if ( y == NULL ) scale = Real1DArray_Item ( c, 0 ) ;
        /* . Cross-correlation. */
        else             scale = sqrt ( fabs ( Real1DArray_Dot ( x, x, NULL ) * Real1DArray_Dot ( y, y, NULL ) ) ) ;
        /* . Scaling. */
        if ( fabs ( scale ) > small ) Real1DArray_Scale ( c, 1.0 / scale ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Remove the mean from an array of data.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void Correlation_SimpleRemoveMean ( Real1DArray *x )
{
    auto Integer n ;
    n = Real1DArray_Length ( x ) ;
    if ( ( x != NULL ) && ( n > 0 ) )
    {
        auto Real mean ;
        mean = Real1DArray_Sum ( x ) / ( ( Real ) n ) ;
        Real1DArray_AddScalar ( x, -mean ) ;
    }
}
