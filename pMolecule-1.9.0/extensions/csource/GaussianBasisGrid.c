/*------------------------------------------------------------------------------
! . File      : GaussianBasisGrid.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . Grid procedures for Gaussian basis sets.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>

# include "GaussianBasisGrid.h"

/* . Note here that the order of derivatives does not have to correspond to the order of the Cartesian basis functions. */

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void GaussianBasis_GridValue   ( const GaussianBasis *iBasis, const Real *rI, const Real *rG, const Integer stride, Real *f ) ;
static void GaussianBasis_GridValueD  ( const GaussianBasis *iBasis, const Real *rI, const Real *rG, const Integer stride, Real *f, Real *fX, Real *fY, Real *fZ ) ;
static void GaussianBasis_GridValueD2 ( const GaussianBasis *iBasis, const Real *rI, const Real *rG, const Integer stride, Real *f, Real *fX, Real *fY, Real *fZ ,
                                                                                               Real *fXX, Real *fXY, Real *fXZ, Real *fYY, Real *fYZ,  Real *fZZ ) ;
static void GaussianBasis_GridValueD3 ( const GaussianBasis *iBasis, const Real *rI, const Real *rG, const Integer stride ,
                                                                                    Real *f, Real *fX, Real *fY, Real *fZ ,
                                                        Real *fXX, Real *fXY, Real *fXZ, Real *fYY, Real *fYZ,  Real *fZZ ,
                                                               Real *fXXX, Real *fXXY, Real *fXXZ, Real *fXYY, Real *fXYZ ,
                                                               Real *fXZZ, Real *fYYY, Real *fYYZ, Real *fYZZ, Real *fZZZ ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the values of the basis functions and their derivatives at a block of grid points.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Probably better to do for multiple points at once rather than separately. */
/* . Storage starts with the current value of numberOfFunctions in data. */
void GaussianBasis_GridPointValues ( const GaussianBasis         *iBasis ,
                                     const Real                  *rI     ,
                                     const Coordinates3          *rG     ,
                                           GridFunctionDataBlock *data   )
{
    if ( ( iBasis != NULL ) && ( rI != NULL ) && ( rG != NULL ) && ( data != NULL ) )
    {
        auto Integer fStart, g ;
        fStart = data->numberOfFunctions ;
        switch ( data->order )
        {
            case 0:
                for ( g = 0 ; g < rG->length0 ; g++ )
                {
                    GaussianBasis_GridValue ( iBasis ,  rI ,
                                              Coordinates3_RowPointer ( rG, g ) ,
                                              data->numberOfPoints ,
                                              Real2DArray_ItemPointer ( data->f, fStart, g ) ) ;
                }
                break ;
            case 1:
                for ( g = 0 ; g < rG->length0 ; g++ )
                {
                    GaussianBasis_GridValueD ( iBasis , rI ,                               
                                               Coordinates3_RowPointer ( rG, g ) ,   
                                               data->numberOfPoints ,   
                                               Real2DArray_ItemPointer ( data->f , fStart, g ) ,   
                                               Real2DArray_ItemPointer ( data->fX, fStart, g ) ,   
                                               Real2DArray_ItemPointer ( data->fY, fStart, g ) ,   
                                               Real2DArray_ItemPointer ( data->fZ, fStart, g ) ) ; 
                }
                break ;
            case 2:
                for ( g = 0 ; g < rG->length0 ; g++ )
                {
                    GaussianBasis_GridValueD2 ( iBasis , rI ,
                                                Coordinates3_RowPointer ( rG, g ) ,
                                                data->numberOfPoints ,
                                                Real2DArray_ItemPointer ( data->f  , fStart, g ) ,
                                                Real2DArray_ItemPointer ( data->fX , fStart, g ) ,
                                                Real2DArray_ItemPointer ( data->fY , fStart, g ) ,
                                                Real2DArray_ItemPointer ( data->fZ , fStart, g ) ,
                                                Real2DArray_ItemPointer ( data->fXX, fStart, g ) ,
                                                Real2DArray_ItemPointer ( data->fXY, fStart, g ) ,
                                                Real2DArray_ItemPointer ( data->fXZ, fStart, g ) ,
                                                Real2DArray_ItemPointer ( data->fYY, fStart, g ) ,
                                                Real2DArray_ItemPointer ( data->fYZ, fStart, g ) ,
                                                Real2DArray_ItemPointer ( data->fZZ, fStart, g ) ) ;
                }
                break ;
            case 3:
                for ( g = 0 ; g < rG->length0 ; g++ )
                {
                    GaussianBasis_GridValueD3 ( iBasis , rI ,
                                                Coordinates3_RowPointer ( rG, g ) ,
                                                data->numberOfPoints ,
                                                Real2DArray_ItemPointer ( data->f   , fStart, g ) ,
                                                Real2DArray_ItemPointer ( data->fX  , fStart, g ) ,
                                                Real2DArray_ItemPointer ( data->fY  , fStart, g ) ,
                                                Real2DArray_ItemPointer ( data->fZ  , fStart, g ) ,
                                                Real2DArray_ItemPointer ( data->fXX , fStart, g ) ,
                                                Real2DArray_ItemPointer ( data->fXY , fStart, g ) ,
                                                Real2DArray_ItemPointer ( data->fXZ , fStart, g ) ,
                                                Real2DArray_ItemPointer ( data->fYY , fStart, g ) ,
                                                Real2DArray_ItemPointer ( data->fYZ , fStart, g ) ,
                                                Real2DArray_ItemPointer ( data->fZZ , fStart, g ) ,
                                                Real2DArray_ItemPointer ( data->fXXX, fStart, g ) ,
                                                Real2DArray_ItemPointer ( data->fXXY, fStart, g ) ,
                                                Real2DArray_ItemPointer ( data->fXXZ, fStart, g ) ,
                                                Real2DArray_ItemPointer ( data->fXYY, fStart, g ) ,
                                                Real2DArray_ItemPointer ( data->fXYZ, fStart, g ) ,
                                                Real2DArray_ItemPointer ( data->fXZZ, fStart, g ) ,
                                                Real2DArray_ItemPointer ( data->fYYY, fStart, g ) ,
                                                Real2DArray_ItemPointer ( data->fYYZ, fStart, g ) ,
                                                Real2DArray_ItemPointer ( data->fYZZ, fStart, g ) ,
                                                Real2DArray_ItemPointer ( data->fZZZ, fStart, g ) ) ;
                }
        }
        /* . Increment the function number. */
        data->numberOfFunctions += iBasis->nbasisw ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the values of the basis functions at a given point.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void GaussianBasis_GridValue ( const GaussianBasis *iBasis, const Real *rI, const Real *rG, const Integer stride, Real *f )
{
    if ( ( iBasis != NULL ) && ( rI != NULL ) && ( rG != NULL ) && ( f != NULL ) )
    {
        Integer i, icbfind, ip, iShell, iStart, ix, iy, iz, ncfunci ;
        Real    dx, dy, dz, eip, r2 ;
        Real    e0[MAXCBF], g0[MAXCBF], x0[MAXAMP1], y0[MAXAMP1], z0[MAXAMP1] ;

        /* . Get the distance squared between atom centers. */
        dx = rG[0] - rI[0] ;
        dy = rG[1] - rI[1] ;
        dz = rG[2] - rI[2] ;
        r2 = dx * dx + dy * dy + dz * dz ;

        /* . Check for negligible contributions. */
        /* if ( r2 >= iBasis->maximum_range ) return ; */

        /* . Form the angular functions. */
        x0[0] = 1.0e+00 ; y0[0] = 1.0e+00 ; z0[0] = 1.0e+00 ;
        for ( i = 1 ; i <= iBasis->maximum_angularmomentum ; i++ )
        {
            x0[i] = dx * x0[i-1] ;
            y0[i] = dy * y0[i-1] ;
            z0[i] = dz * z0[i-1] ;
        }

        /* . Outer loop over shells. */
        for ( iShell = 0 ; iShell < iBasis->nshells ; iShell++ )
        {
            /* . Get information about the shell. */
            icbfind = iBasis->shells[iShell].type->cbfindex ;
            ncfunci = iBasis->shells[iShell].type->ncbf     ;

            /* . Check for negligible contributions. */
            /* if ( r2 >= iBasis->shells[iShell].maximum_range ) return ; */

            /* . Form the exponential factors. */
            for ( i = 0 ; i < ncfunci ; i++ ) { e0[i] = 0.0e+00 ; }
            for ( ip = 0; ip < iBasis->shells[iShell].nprimitives ; ip ++ )
            {
                eip = exp ( - iBasis->shells[iShell].primitives[ip].exponent * r2 ) ;
                for ( i = 0 ; i < ncfunci ; i++ ) { e0[i] += iBasis->shells[iShell].primitives[ip].ccbf[i] * eip ; }
            }

            /* . Form the Cartesian function values. */
            for ( i = 0 ; i < ncfunci ; i++ )
            {
                ix = CBFPOWX[i+icbfind] ;
                iy = CBFPOWY[i+icbfind] ;
                iz = CBFPOWZ[i+icbfind] ;
                g0[i] = x0[ix] * y0[iy] * z0[iz] * e0[i] ;
            }

            /* . Transform to spherical harmonics. */
/*            if ( iBasis->QTOSPHERICAL ) Integral_Real1DArray_Transform ( iangmom, ncfunci, nfunci, g0, 1, 0 ) ; */

            /* . Put the values in the proper place. */
            iStart = iBasis->shells[iShell].nstartw * stride ;
            for ( i = 0 ; i < iBasis->shells[iShell].nbasisw ; i++ ) { f[iStart+i*stride] = g0[i] ; }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the values of the basis functions at a given point and their first derivatives.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void GaussianBasis_GridValueD ( const GaussianBasis *iBasis, const Real *rI, const Real *rG, const Integer stride, Real *f, Real *fX, Real *fY, Real *fZ )
{
    if ( ( iBasis != NULL ) && ( rI != NULL ) && ( rG != NULL ) && ( f != NULL ) && ( fX != NULL ) && ( fY != NULL ) && ( fZ != NULL ) )
    {
        Integer i, icbfind, ii, ip, iShell, iStart, ix, iy, iz, ncfunci ;
        Real    dx, dy, dz, eip, r2 ;
        Real    e0[MAXCBF], e1[MAXCBF], g01[4*MAXCBF], x0[MAXAMP2], y0[MAXAMP2], z0[MAXAMP2] ,
                                                       x1[MAXAMP2], y1[MAXAMP2], z1[MAXAMP2] ;

        /* . Get the distance squared between atom centers. */
        dx = rG[0] - rI[0] ;
        dy = rG[1] - rI[1] ;
        dz = rG[2] - rI[2] ;
        r2 = dx * dx + dy * dy + dz * dz ;

        /* . Check for negligible contributions. */
        /* if ( r2 >= iBasis->maximum_range ) return ; */

        /* . Form the angular functions. */
        x0[0] = 1.0e+00 ; y0[0] = 1.0e+00 ; z0[0] = 1.0e+00 ;
        x1[0] = 0.0e+00 ; y1[0] = 0.0e+00 ; z1[0] = 0.0e+00 ;
        for ( i = 1 ; i <= iBasis->maximum_angularmomentum + 1; i++ )
        {
            x0[i] = dx * x0[i-1] ;
            y0[i] = dy * y0[i-1] ;
            z0[i] = dz * z0[i-1] ;
            x1[i] = ( Real ) ( i ) * x0[i-1] ;
            y1[i] = ( Real ) ( i ) * y0[i-1] ;
            z1[i] = ( Real ) ( i ) * z0[i-1] ;
        }

        /* . Outer loop over shells. */
        for ( iShell = 0 ; iShell < iBasis->nshells ; iShell++ )
        {
            /* . Get information about the shell. */
            icbfind = iBasis->shells[iShell].type->cbfindex ;
            ncfunci = iBasis->shells[iShell].type->ncbf     ;

            /* . Check for negligible contributions. */
            /* if ( r2 >= iBasis->shells[iShell].maximum_range ) return ; */

            /* . Form the exponential factors. */
            for ( i = 0 ; i < ncfunci ; i++ ) { e0[i] = 0.0e+00 ; e1[i] = 0.0e+00 ; }
            for ( ip = 0; ip < iBasis->shells[iShell].nprimitives ; ip ++ )
            {
                eip = exp ( - iBasis->shells[iShell].primitives[ip].exponent * r2 ) ;
                for ( i = 0 ; i < ncfunci ; i++ )
                {
                    e0[i] += iBasis->shells[iShell].primitives[ip].ccbf[i]  * eip ;
                    e1[i] += iBasis->shells[iShell].primitives[ip].ccbf[i]  *
                             iBasis->shells[iShell].primitives[ip].exponent * eip ;
                }
            }
            for ( i = 0 ; i < ncfunci ; i++ ) { e1[i] *= -2.0e+00 ; }

            /* . Form the Cartesian function values. */
            for ( i = 0 ; i < ncfunci ; i++ )
            {
                ix = CBFPOWX[i+icbfind] ;
                iy = CBFPOWY[i+icbfind] ;
                iz = CBFPOWZ[i+icbfind] ;
                g01[i]          = x0[ix] * y0[iy] * z0[iz] * e0[i] ;
                g01[i+  MAXCBF] = ( x1[ix] * e0[i] + x0[ix+1] * e1[i] ) * y0[iy] * z0[iz] ;
                g01[i+2*MAXCBF] = ( y1[iy] * e0[i] + y0[iy+1] * e1[i] ) * x0[ix] * z0[iz] ;
                g01[i+3*MAXCBF] = ( z1[iz] * e0[i] + z0[iz+1] * e1[i] ) * x0[ix] * y0[iy] ;
            }

            /* . Transform to spherical harmonics. */
/*            if ( iBasis->QTOSPHERICAL ) Integral_Real1DArray_Transform ( iangmom, ncfunci, nfunci, g01, 4, MAXCBF ) ; */

            /* . Put the values in the proper place. */
            iStart = iBasis->shells[iShell].nstartw * stride ;
            for ( i = 0 ; i < iBasis->shells[iShell].nbasisw ; i++ )
            {
                ii = iStart + i*stride ;
                f [ii] = g01[i] ;
                fX[ii] = g01[i+  MAXCBF] ;
                fY[ii] = g01[i+2*MAXCBF] ;
                fZ[ii] = g01[i+3*MAXCBF] ;
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the values of the basis functions at a given point and their first and second derivatives.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void GaussianBasis_GridValueD2 ( const GaussianBasis *iBasis, const Real *rI, const Real *rG, const Integer stride, Real *f, Real *fX, Real *fY, Real *fZ,
                                                                                              Real *fXX, Real *fXY, Real *fXZ, Real *fYY, Real *fYZ,  Real *fZZ )
{
    if ( ( iBasis != NULL ) && ( rI != NULL ) && ( rG != NULL ) && ( f != NULL ) && ( fX != NULL ) && ( fY != NULL ) && ( fZ != NULL ) ) /* . More checks here? */
    {
        Integer i, icbfind, ii, ip, iShell, iStart, ix, iy, iz, ncfunci ;
        Real    dx, dy, dz, e, eip, ee, r2 ;
        Real    e0[MAXCBF], e1[MAXCBF], e2[MAXCBF], g012[10*MAXCBF], x0[MAXAMP3], y0[MAXAMP3], z0[MAXAMP3] ,
                                                                     x1[MAXAMP3], y1[MAXAMP3], z1[MAXAMP3] ,
                                                                     x2[MAXAMP3], y2[MAXAMP3], z2[MAXAMP3] ;

        /* . Get the distance squared between atom centers. */
        dx = rG[0] - rI[0] ;
        dy = rG[1] - rI[1] ;
        dz = rG[2] - rI[2] ;
        r2 = dx * dx + dy * dy + dz * dz ;

        /* . Check for negligible contributions. */
        /* if ( r2 >= iBasis->maximum_range ) return ; */

        /* . Form the angular functions. */
        x0[0] = 1.0e+00 ; y0[0] = 1.0e+00 ; z0[0] = 1.0e+00 ;
        x1[0] = 0.0e+00 ; y1[0] = 0.0e+00 ; z1[0] = 0.0e+00 ;
        x2[0] = 0.0e+00 ; y2[0] = 0.0e+00 ; z2[0] = 0.0e+00 ;
        for ( i = 1 ; i <= iBasis->maximum_angularmomentum + 2; i++ )
        {
            x0[i] = dx * x0[i-1] ;
            y0[i] = dy * y0[i-1] ;
            z0[i] = dz * z0[i-1] ;
            x1[i] = ( Real ) ( i ) * x0[i-1] ;
            y1[i] = ( Real ) ( i ) * y0[i-1] ;
            z1[i] = ( Real ) ( i ) * z0[i-1] ;
            x2[i] = ( Real ) ( i ) * x1[i-1] ;
            y2[i] = ( Real ) ( i ) * y1[i-1] ;
            z2[i] = ( Real ) ( i ) * z1[i-1] ;
        }

        /* . Outer loop over shells. */
        for ( iShell = 0 ; iShell < iBasis->nshells ; iShell++ )
        {
            /* . Get information about the shell. */
            icbfind = iBasis->shells[iShell].type->cbfindex ;
            ncfunci = iBasis->shells[iShell].type->ncbf     ;

            /* . Check for negligible contributions. */
            /* if ( r2 >= iBasis->shells[iShell].maximum_range ) return ; */

            /* . Form the exponential factors. */
            for ( i = 0 ; i < ncfunci ; i++ ) { e0[i] = 0.0e+00 ; e1[i] = 0.0e+00 ; e2[i] = 0.0e+00 ; }
            for ( ip = 0; ip < iBasis->shells[iShell].nprimitives ; ip ++ )
            {
                e   = iBasis->shells[iShell].primitives[ip].exponent ;
                ee  = e * e ;
                eip = exp ( - e * r2 ) ;
                for ( i = 0 ; i < ncfunci ; i++ )
                {
                    e0[i] += iBasis->shells[iShell].primitives[ip].ccbf[i]       * eip ;
                    e1[i] += iBasis->shells[iShell].primitives[ip].ccbf[i]  * e  * eip ;
                    e2[i] += iBasis->shells[iShell].primitives[ip].ccbf[i]  * ee * eip ;
                }
            }
            for ( i = 0 ; i < ncfunci ; i++ ) { e1[i] *= -2.0e+00 ; e2[i] *= 4.0e+00 ; }

            /* . Form the Cartesian function values. */
            for ( i = 0 ; i < ncfunci ; i++ )
            {
                ix = CBFPOWX[i+icbfind] ;
                iy = CBFPOWY[i+icbfind] ;
                iz = CBFPOWZ[i+icbfind] ;
                g012[i]          = x0[ix] * y0[iy] * z0[iz] * e0[i] ;
                g012[i+  MAXCBF] = ( x1[ix] * e0[i] + x0[ix+1] * e1[i] ) * y0[iy] * z0[iz] ;
                g012[i+2*MAXCBF] = ( y1[iy] * e0[i] + y0[iy+1] * e1[i] ) * x0[ix] * z0[iz] ;
                g012[i+3*MAXCBF] = ( z1[iz] * e0[i] + z0[iz+1] * e1[i] ) * x0[ix] * y0[iy] ;
                g012[i+4*MAXCBF] = ( x2[ix] * e0[i] + ( dx * x1[ix] + x1[ix+1] ) * e1[i] + x0[ix+2] * e2[i] ) * y0[iy] * z0[iz] ;
                g012[i+5*MAXCBF] = ( x1[ix]*y1[iy]*e0[i] + ( x1[ix]*y0[iy+1] + x0[ix+1]*y1[iy] )*e1[i] + x0[ix+1]*y0[iy+1]*e2[i] ) * z0[iz] ;
                g012[i+6*MAXCBF] = ( x1[ix]*z1[iz]*e0[i] + ( x1[ix]*z0[iz+1] + x0[ix+1]*z1[iz] )*e1[i] + x0[ix+1]*z0[iz+1]*e2[i] ) * y0[iy] ;
                g012[i+7*MAXCBF] = ( y2[iy] * e0[i] + ( dy * y1[iy] + y1[iy+1] ) * e1[i] + y0[iy+2] * e2[i] ) * x0[ix] * z0[iz] ;
                g012[i+8*MAXCBF] = ( y1[iy]*z1[iz]*e0[i] + ( y1[iy]*z0[iz+1] + y0[iy+1]*z1[iz] )*e1[i] + y0[iy+1]*z0[iz+1]*e2[i] ) * x0[ix] ;
                g012[i+9*MAXCBF] = ( z2[iz] * e0[i] + ( dz * z1[iz] + z1[iz+1] ) * e1[i] + z0[iz+2] * e2[i] ) * x0[ix] * y0[iy] ;
            }

            /* . Transform to spherical harmonics. */
/*            if ( iBasis->QTOSPHERICAL ) Integral_Real1DArray_Transform ( iangmom, ncfunci, nfunci, g012, 10, MAXCBF ) ; */

            /* . Put the values in the proper place. */
            iStart = iBasis->shells[iShell].nstartw * stride ;
            for ( i = 0 ; i < iBasis->shells[iShell].nbasisw ; i++ )
            {
                ii = iStart + i*stride ;
                f  [ii] = g012[i] ;
                fX [ii] = g012[i+  MAXCBF] ;
                fY [ii] = g012[i+2*MAXCBF] ;
                fZ [ii] = g012[i+3*MAXCBF] ;
                fXX[ii] = g012[i+4*MAXCBF] ;
                fXY[ii] = g012[i+5*MAXCBF] ;
                fXZ[ii] = g012[i+6*MAXCBF] ;
                fYY[ii] = g012[i+7*MAXCBF] ;
                fYZ[ii] = g012[i+8*MAXCBF] ;
                fZZ[ii] = g012[i+9*MAXCBF] ;
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the values of the basis functions at a given point and their first, second and third derivatives.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void GaussianBasis_GridValueD3 ( const GaussianBasis *iBasis, const Real *rI, const Real *rG, const Integer stride ,
                                                                                    Real *f, Real *fX, Real *fY, Real *fZ ,
                                                        Real *fXX, Real *fXY, Real *fXZ, Real *fYY, Real *fYZ,  Real *fZZ ,
                                                               Real *fXXX, Real *fXXY, Real *fXXZ, Real *fXYY, Real *fXYZ ,
                                                               Real *fXZZ, Real *fYYY, Real *fYYZ, Real *fYZZ, Real *fZZZ )
{
    if ( ( iBasis != NULL ) && ( rI != NULL ) && ( rG != NULL ) && ( f != NULL ) && ( fX != NULL ) && ( fY != NULL ) && ( fZ != NULL ) ) /* . More checks here? */
    {
        Integer i, icbfind, ii, ip, iShell, iStart, ix, iy, iz, ncfunci ;
        Real    dx, dy, dz, e, eip, ee, eee, r2 ;
        Real    e0[MAXCBF], e1[MAXCBF], e2[MAXCBF], e3[MAXCBF], g013[20*MAXCBF], x0[MAXAMP4], y0[MAXAMP4], z0[MAXAMP4] ,
                                                                                 x1[MAXAMP4], y1[MAXAMP4], z1[MAXAMP4] ,
                                                                                 x2[MAXAMP4], y2[MAXAMP4], z2[MAXAMP4] ,
                                                                                 x3[MAXAMP4], y3[MAXAMP4], z3[MAXAMP4] ;

        /* . Get the distance squared between atom centers. */
        dx = rG[0] - rI[0] ;
        dy = rG[1] - rI[1] ;
        dz = rG[2] - rI[2] ;
        r2 = dx * dx + dy * dy + dz * dz ;

        /* . Check for negligible contributions. */
        /* if ( r2 >= iBasis->maximum_range ) return ; */

        /* . Form the angular functions. */
        x0[0] = 1.0e+00 ; y0[0] = 1.0e+00 ; z0[0] = 1.0e+00 ;
        x1[0] = 0.0e+00 ; y1[0] = 0.0e+00 ; z1[0] = 0.0e+00 ;
        x2[0] = 0.0e+00 ; y2[0] = 0.0e+00 ; z2[0] = 0.0e+00 ;
        x3[0] = 0.0e+00 ; y3[0] = 0.0e+00 ; z3[0] = 0.0e+00 ;
        for ( i = 1 ; i <= iBasis->maximum_angularmomentum + 3 ; i++ )
        {
            x0[i] = dx * x0[i-1] ;
            y0[i] = dy * y0[i-1] ;
            z0[i] = dz * z0[i-1] ;
            x1[i] = ( Real ) ( i ) * x0[i-1] ;
            y1[i] = ( Real ) ( i ) * y0[i-1] ;
            z1[i] = ( Real ) ( i ) * z0[i-1] ;
            x2[i] = ( Real ) ( i ) * x1[i-1] ;
            y2[i] = ( Real ) ( i ) * y1[i-1] ;
            z2[i] = ( Real ) ( i ) * z1[i-1] ;
            x3[i] = ( Real ) ( i ) * x2[i-1] ;
            y3[i] = ( Real ) ( i ) * y2[i-1] ;
            z3[i] = ( Real ) ( i ) * z2[i-1] ;
        }

        /* . Outer loop over shells. */
        for ( iShell = 0 ; iShell < iBasis->nshells ; iShell++ )
        {
            /* . Get information about the shell. */
            icbfind = iBasis->shells[iShell].type->cbfindex ;
            ncfunci = iBasis->shells[iShell].type->ncbf     ;

            /* . Check for negligible contributions. */
            /* if ( r2 >= iBasis->shells[iShell].maximum_range ) return ; */

            /* . Form the exponential factors. */
            for ( i = 0 ; i < ncfunci ; i++ ) { e0[i] = 0.0e+00 ; e1[i] = 0.0e+00 ; e2[i] = 0.0e+00 ; e3[i] = 0.0e+00 ; }
            for ( ip = 0; ip < iBasis->shells[iShell].nprimitives ; ip ++ )
            {
                e   = iBasis->shells[iShell].primitives[ip].exponent ;
                ee  = e  * e ;
                eee = ee * e ;
                eip = exp ( - e * r2 ) ;
                for ( i = 0 ; i < ncfunci ; i++ )
                {
                    e0[i] += iBasis->shells[iShell].primitives[ip].ccbf[i]        * eip ;
                    e1[i] += iBasis->shells[iShell].primitives[ip].ccbf[i]  * e   * eip ;
                    e2[i] += iBasis->shells[iShell].primitives[ip].ccbf[i]  * ee  * eip ;
                    e3[i] += iBasis->shells[iShell].primitives[ip].ccbf[i]  * eee * eip ;
                }
            }
            for ( i = 0 ; i < ncfunci ; i++ ) { e1[i] *= -2.0e+00 ; e2[i] *= 4.0e+00 ; e3[i] *= -8.0e+00 ; }

            /* . Form the Cartesian function values. */
            for ( i = 0 ; i < ncfunci ; i++ )
            {
                ix = CBFPOWX[i+icbfind] ;
                iy = CBFPOWY[i+icbfind] ;
                iz = CBFPOWZ[i+icbfind] ;
                g013[i]          = x0[ix] * y0[iy] * z0[iz] * e0[i] ;
                g013[i+   MAXCBF] = ( x1[ix] * e0[i] + x0[ix+1] * e1[i] ) * y0[iy] * z0[iz] ;
                g013[i+ 2*MAXCBF] = ( y1[iy] * e0[i] + y0[iy+1] * e1[i] ) * x0[ix] * z0[iz] ;
                g013[i+ 3*MAXCBF] = ( z1[iz] * e0[i] + z0[iz+1] * e1[i] ) * x0[ix] * y0[iy] ;
                g013[i+ 4*MAXCBF] = ( x2[ix] * e0[i] + ( dx * x1[ix] + x1[ix+1] ) * e1[i] + x0[ix+2] * e2[i] ) * y0[iy] * z0[iz] ;
                g013[i+ 5*MAXCBF] = ( x1[ix]*y1[iy]*e0[i] + ( x1[ix]*y0[iy+1] + x0[ix+1]*y1[iy] )*e1[i] + x0[ix+1]*y0[iy+1]*e2[i] ) * z0[iz] ;
                g013[i+ 6*MAXCBF] = ( x1[ix]*z1[iz]*e0[i] + ( x1[ix]*z0[iz+1] + x0[ix+1]*z1[iz] )*e1[i] + x0[ix+1]*z0[iz+1]*e2[i] ) * y0[iy] ;
                g013[i+ 7*MAXCBF] = ( y2[iy] * e0[i] + ( dy * y1[iy] + y1[iy+1] ) * e1[i] + y0[iy+2] * e2[i] ) * x0[ix] * z0[iz] ;
                g013[i+ 8*MAXCBF] = ( y1[iy]*z1[iz]*e0[i] + ( y1[iy]*z0[iz+1] + y0[iy+1]*z1[iz] )*e1[i] + y0[iy+1]*z0[iz+1]*e2[i] ) * x0[ix] ;
                g013[i+ 9*MAXCBF] = ( z2[iz] * e0[i] + ( dz * z1[iz] + z1[iz+1] ) * e1[i] + z0[iz+2] * e2[i] ) * x0[ix] * y0[iy] ;
                g013[i+10*MAXCBF] = ( x3[ix] * e0[i] + ( x1[ix] + 2.0e+00*dx*x2[ix] + x2[ix+1] )*e1[i] + ( dx*dx*x1[ix] + dx*x1[ix+1] + x1[ix+2] )*e2[i] + x0[ix+3] * e3[i] ) * y0[iy] * z0[iz] ;
                g013[i+11*MAXCBF] = ( x2[ix]*y1[iy]*e0[i] + ( ( dx*x1[ix] + x1[ix+1] )*y1[iy]   + x2[ix]*y0[iy+1] )*e1[i] +
                                                            ( ( dx*x1[ix] + x1[ix+1] )*y0[iy+1] + x0[ix+2]*y1[iy] )*e2[i] + x0[ix+2]*y0[iy+1]*e3[i] ) * z0[iz] ;
                g013[i+12*MAXCBF] = ( x2[ix]*z1[iz]*e0[i] + ( ( dx*x1[ix] + x1[ix+1] )*z1[iz]   + x2[ix]*z0[iz+1] )*e1[i] +
                                                            ( ( dx*x1[ix] + x1[ix+1] )*z0[iz+1] + x0[ix+2]*z1[iz] )*e2[i] + x0[ix+2]*z0[iz+1]*e3[i] ) * y0[iy] ;
                g013[i+13*MAXCBF] = ( y2[iy]*x1[ix]*e0[i] + ( ( dy*y1[iy] + y1[iy+1] )*x1[ix]   + y2[iy]*x0[ix+1] )*e1[i] +
                                                            ( ( dy*y1[iy] + y1[iy+1] )*x0[ix+1] + y0[iy+2]*x1[ix] )*e2[i] + y0[iy+2]*x0[ix+1]*e3[i] ) * z0[iz] ;
                g013[i+14*MAXCBF] = ( x1[ix]*y1[iy]*z1[iz]*e0[i] + ( x1[ix]*y1[iy  ]*z0[iz+1] + x1[ix  ]*y0[iy+1]*z1[iz  ] + x0[ix+1]*y1[iy  ]*z1[iz] )*e1[i] +
                                                                   ( x1[ix]*y0[iy+1]*z0[iz+1] + x0[ix+1]*y1[iy  ]*z0[iz+1] + x0[ix+1]*y0[iy+1]*z1[iz] )*e2[i] + x0[ix+1]*y0[ix+1]*z0[iz+1]*e3[i] ) ;
                g013[i+15*MAXCBF] = ( z2[iz]*x1[ix]*e0[i] + ( ( dz*z1[iz] + z1[iz+1] )*x1[ix]   + z2[iz]*x0[ix+1] )*e1[i] +
                                                            ( ( dz*z1[iz] + z1[iz+1] )*x0[ix+1] + z0[iz+2]*x1[ix] )*e2[i] + z0[iz+2]*x0[ix+1]*e3[i] ) * y0[iy] ;
                g013[i+16*MAXCBF] = ( y3[iy] * e0[i] + ( y1[iy] + 2.0e+00*dy*y2[iy] + y2[iy+1] )*e1[i] + ( dy*dy*y1[iy] + dy*y1[iy+1] + y1[iy+2] )*e2[i] + y0[iy+3] * e3[i] ) * x0[ix] * z0[iz] ;
                g013[i+17*MAXCBF] = ( y2[iy]*z1[iz]*e0[i] + ( ( dy*y1[iy] + y1[iy+1] )*z1[iz]   + y2[iy]*z0[iz+1] )*e1[i] +
                                                            ( ( dy*y1[iy] + y1[iy+1] )*z0[iz+1] + y0[iy+2]*z1[iz] )*e2[i] + y0[iy+2]*z0[iz+1]*e3[i] ) * x0[ix] ;
                g013[i+18*MAXCBF] = ( z2[iz]*y1[iy]*e0[i] + ( ( dz*z1[iz] + z1[iz+1] )*y1[iy]   + z2[iz]*y0[iy+1] )*e1[i] +
                                                            ( ( dz*z1[iz] + z1[iz+1] )*y0[iy+1] + z0[iz+2]*y1[iy] )*e2[i] + z0[iz+2]*y0[iy+1]*e3[i] ) * x0[ix] ;
                g013[i+19*MAXCBF] = ( z3[iz] * e0[i] + ( z1[iz] + 2.0e+00*dz*z2[iz] + z2[iz+1] )*e1[i] + ( dz*dz*z1[iz] + dz*z1[iz+1] + z1[iz+2] )*e2[i] + z0[iz+3] * e3[i] ) * x0[ix] * y0[iy] ;
            }

            /* . Transform to spherical harmonics. */
/*            if ( iBasis->QTOSPHERICAL ) Integral_Real1DArray_Transform ( iangmom, ncfunci, nfunci, g013, 20, MAXCBF ) ; */

            /* . Put the values in the proper place. */
            iStart = iBasis->shells[iShell].nstartw * stride ;
            for ( i = 0 ; i < iBasis->shells[iShell].nbasisw ; i++ )
            {
                ii = iStart + i*stride ;
                f   [ii] = g013[i] ;
                fX  [ii] = g013[i+   MAXCBF] ;
                fY  [ii] = g013[i+ 2*MAXCBF] ;
                fZ  [ii] = g013[i+ 3*MAXCBF] ;
                fXX [ii] = g013[i+ 4*MAXCBF] ;
                fXY [ii] = g013[i+ 5*MAXCBF] ;
                fXZ [ii] = g013[i+ 6*MAXCBF] ;
                fYY [ii] = g013[i+ 7*MAXCBF] ;
                fYZ [ii] = g013[i+ 8*MAXCBF] ;
                fZZ [ii] = g013[i+ 9*MAXCBF] ;
                fXXX[ii] = g013[i+10*MAXCBF] ;
                fXXY[ii] = g013[i+11*MAXCBF] ;
                fXXZ[ii] = g013[i+12*MAXCBF] ;
                fXYY[ii] = g013[i+13*MAXCBF] ;
                fXYZ[ii] = g013[i+14*MAXCBF] ;
                fXZZ[ii] = g013[i+15*MAXCBF] ;
                fYYY[ii] = g013[i+16*MAXCBF] ;
                fYYZ[ii] = g013[i+17*MAXCBF] ;
                fYZZ[ii] = g013[i+18*MAXCBF] ;
                fZZZ[ii] = g013[i+19*MAXCBF] ;
            }
        }
    }
}
