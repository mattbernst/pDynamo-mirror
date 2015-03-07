/*------------------------------------------------------------------------------
! . File      : SymmetryParameterGradients.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==============================================================================
! . Notes:
!
!   Both r/M and f/M representations are catered for but care should be taken
!   to ensure that they are used consistently.
!
!=============================================================================*/

# include <math.h>

# include "Memory.h"
# include "SymmetryParameterGradients.h"
# include "Units.h"

/*------------------------------------------------------------------------------
! . Allocation.
!-----------------------------------------------------------------------------*/
SymmetryParameterGradients *SymmetryParameterGradients_Allocate ( void )
{
    SymmetryParameterGradients *self = NULL ;
    self           = ( SymmetryParameterGradients * ) Memory_Allocate ( sizeof ( SymmetryParameterGradients ) ) ;
    self->dEda     = 0.0e+00 ;
    self->dEdb     = 0.0e+00 ;
    self->dEdc     = 0.0e+00 ;
    self->dEdalpha = 0.0e+00 ;
    self->dEdbeta  = 0.0e+00 ;
    self->dEdgamma = 0.0e+00 ;
    self->dEdM     = Matrix33_Allocate ( ) ;
    Matrix33_Set ( self->dEdM, 0.0e+00 ) ;
    return self ;
}

/*------------------------------------------------------------------------------
! . Deallocation.
!-----------------------------------------------------------------------------*/
void SymmetryParameterGradients_Deallocate ( SymmetryParameterGradients **self )
{
    if ( (*self) != NULL )
    {
        Matrix33_Deallocate ( &((*self)->dEdM) ) ;
        Memory_Deallocate   (   (*self) ) ;
    }
}

/*------------------------------------------------------------------------------
! . Convert dEdM to dEda, dEdb, dEdc, dEdalpha, dEdbeta, dEdgamma.
! . This depends upon how M is defined in symmetryParameters.
! . The procedure works for both the r/M and f/M representations.
!-----------------------------------------------------------------------------*/
void SymmetryParameterGradients_CrystalDerivatives ( SymmetryParameterGradients *self, const SymmetryParameters *symmetryParameters )
{
    if ( ( self != NULL ) && ( symmetryParameters != NULL ) )
    {
        auto Real alpha, beta, cosalpha, cosbeta, cosgamma, fact12, fact22, gamma, singamma ;
        /* . Some factors. */
        alpha    = symmetryParameters->alpha * UNITS_ANGLE_DEGREES_TO_RADIANS ;
        beta     = symmetryParameters->beta  * UNITS_ANGLE_DEGREES_TO_RADIANS ;
        gamma    = symmetryParameters->gamma * UNITS_ANGLE_DEGREES_TO_RADIANS ;
        cosalpha = cos ( alpha ) ;
        cosbeta  = cos ( beta  ) ;
        cosgamma = cos ( gamma ) ;
        singamma = sin ( gamma ) ;
        fact12   = ( cosalpha - cosbeta * cosgamma ) ;
        fact22   = sqrt ( 1.0e+00 - cosalpha * cosalpha - cosbeta * cosbeta - cosgamma * cosgamma + 2.0e+00 * cosalpha * cosbeta * cosgamma ) ;

        /* . The derivatives - a, b, c. */
        self->dEda =            Matrix33_Item ( self->dEdM, 0, 0 ) ;
        self->dEdb = cosgamma * Matrix33_Item ( self->dEdM, 0, 1 ) + singamma * Matrix33_Item ( self->dEdM, 1, 1 ) ;
        self->dEdc = cosbeta  * Matrix33_Item ( self->dEdM, 0, 2 ) + fact12   * Matrix33_Item ( self->dEdM, 1, 2 ) / singamma +
                                                                                fact22   * Matrix33_Item ( self->dEdM, 2, 2 ) / singamma ;

        /* . The derivatives - alpha, beta, gamma. */
        self->dEdalpha = - ( Matrix33_Item ( self->dEdM, 1, 2 ) - fact12 * Matrix33_Item ( self->dEdM, 2, 2 ) / fact22 ) *
                           ( symmetryParameters->c * sin ( alpha ) * UNITS_ANGLE_DEGREES_TO_RADIANS ) / singamma ;
        self->dEdbeta  = - ( singamma * Matrix33_Item ( self->dEdM, 0, 2 ) - cosgamma * Matrix33_Item ( self->dEdM, 1, 2 ) +
                           ( cosalpha * cosgamma - cosbeta ) * Matrix33_Item ( self->dEdM, 2, 2 ) / fact22 ) *
                           ( symmetryParameters->c * sin ( beta ) * UNITS_ANGLE_DEGREES_TO_RADIANS ) / singamma ;
        self->dEdgamma =   ( symmetryParameters->b * ( - singamma * Matrix33_Item ( self->dEdM, 0, 1 ) + cosgamma * Matrix33_Item ( self->dEdM, 1, 1 ) ) +
                             symmetryParameters->c * ( ( cosbeta - fact12 * cosgamma / ( singamma * singamma ) ) * Matrix33_Item ( self->dEdM, 1, 2 ) +
                                                     ( ( cosgamma - cosalpha * cosbeta ) / fact22 - fact22 * cosgamma / ( singamma * singamma ) ) *
                                                     Matrix33_Item ( self->dEdM, 2, 2 ) ) ) * UNITS_ANGLE_DEGREES_TO_RADIANS ;
    }
}

/*------------------------------------------------------------------------------
! . Convert r/M to f/M derivatives (in-place).
! . Will need selection.
!-----------------------------------------------------------------------------*/
void SymmetryParameterGradients_FractionalDerivatives ( SymmetryParameterGradients *self, const SymmetryParameters *symmetryParameters,
                                                                             const Coordinates3 *coordinates3, Coordinates3 *gradients3 )
{
    if ( ( self != NULL ) && ( symmetryParameters != NULL ) && ( coordinates3 != NULL ) && ( gradients3 != NULL ) && ( coordinates3->length0 == gradients3->length0 ) )
    {
        auto Real fx, fy, fz, gx, gy, gz, im00, im01, im02, im10, im11, im12, im20, im21, im22, m00, m01, m02, m10, m11, m12, m20, m21, m22, x, y, z ;
        auto Integer    i ;

        /* . Get local variables for M and inverseM. */
        m00  = Matrix33_Item ( symmetryParameters->M, 0, 0 ) ;
        m01  = Matrix33_Item ( symmetryParameters->M, 0, 1 ) ;
        m02  = Matrix33_Item ( symmetryParameters->M, 0, 2 ) ;
        m10  = Matrix33_Item ( symmetryParameters->M, 1, 0 ) ;
        m11  = Matrix33_Item ( symmetryParameters->M, 1, 1 ) ;
        m12  = Matrix33_Item ( symmetryParameters->M, 1, 2 ) ;
        m20  = Matrix33_Item ( symmetryParameters->M, 2, 0 ) ;
        m21  = Matrix33_Item ( symmetryParameters->M, 2, 1 ) ;
        m22  = Matrix33_Item ( symmetryParameters->M, 2, 2 ) ;
        im00 = Matrix33_Item ( symmetryParameters->inverseM, 0, 0 ) ;
        im01 = Matrix33_Item ( symmetryParameters->inverseM, 0, 1 ) ;
        im02 = Matrix33_Item ( symmetryParameters->inverseM, 0, 2 ) ;
        im10 = Matrix33_Item ( symmetryParameters->inverseM, 1, 0 ) ;
        im11 = Matrix33_Item ( symmetryParameters->inverseM, 1, 1 ) ;
        im12 = Matrix33_Item ( symmetryParameters->inverseM, 1, 2 ) ;
        im20 = Matrix33_Item ( symmetryParameters->inverseM, 2, 0 ) ;
        im21 = Matrix33_Item ( symmetryParameters->inverseM, 2, 1 ) ;
        im22 = Matrix33_Item ( symmetryParameters->inverseM, 2, 2 ) ;

        /* . Loop over the coordinates/gradients. */
        for ( i = 0 ; i < coordinates3->length0 ; i++ )
        {
            Coordinates3_GetRow ( coordinates3, i,  x,  y,  z ) ;
            Coordinates3_GetRow   ( gradients3,   i, gx, gy, gz ) ;

            /* . Get f. */
            fx = im00 * x + im01 * y + im02 * z ;
            fy = im10 * x + im11 * y + im12 * z ;
            fz = im20 * x + im21 * y + im22 * z ;

            /* . dEdM. */
            Matrix33_Item ( self->dEdM, 0, 0 ) += fx * gx ;
            Matrix33_Item ( self->dEdM, 0, 1 ) += fy * gx ;
            Matrix33_Item ( self->dEdM, 0, 2 ) += fz * gx ;
            Matrix33_Item ( self->dEdM, 1, 0 ) += fx * gy ;
            Matrix33_Item ( self->dEdM, 1, 1 ) += fy * gy ;
            Matrix33_Item ( self->dEdM, 1, 2 ) += fz * gy ;
            Matrix33_Item ( self->dEdM, 2, 0 ) += fx * gz ;
            Matrix33_Item ( self->dEdM, 2, 1 ) += fy * gz ;
            Matrix33_Item ( self->dEdM, 2, 2 ) += fz * gz ;

            /* . dEdf. */
            Coordinates3_Item ( gradients3, i, 0 ) = m00 * gx + m10 * gy + m20 * gz ;
            Coordinates3_Item ( gradients3, i, 1 ) = m01 * gx + m11 * gy + m21 * gz ;
            Coordinates3_Item ( gradients3, i, 2 ) = m02 * gx + m12 * gy + m22 * gz ;
        }
    }
}

/*------------------------------------------------------------------------------
! . Calculate the derivatives due to image terms (r/M formalism).
! . |transformation3| is the fractional transformation (without M).
! . Will need selection.
!-----------------------------------------------------------------------------*/
void SymmetryParameterGradients_ImageDerivatives ( SymmetryParameterGradients *self, const SymmetryParameters *symmetryParameters,
                                                                                                const Transformation3 *transformation3,
                                                                 const Coordinates3 *coordinates3, const Coordinates3 *gradients3  )
{
    if ( ( self != NULL ) && ( symmetryParameters != NULL ) && ( transformation3 != NULL ) && ( coordinates3 != NULL ) && ( gradients3 != NULL ) && ( coordinates3->length0 == gradients3->length0 ) )
    {
        auto Real         dx, dy, dz, gx, gy, gz, r00, r01, r02, r10, r11, r12, r20, r21, r22, sum, t, x, y, z ;
        auto Integer            a, b, i ;
        auto Matrix33 *di, *ms, *si ;

        /* . Allocate space. */
        di = Matrix33_Allocate ( ) ;
        ms = Matrix33_Allocate ( ) ;
        si = Matrix33_Allocate ( ) ;

        /* . Intermediate quantities (MS and SM^-1). */
        Matrix33_CopyTo ( transformation3->rotation, ms, NULL ) ; Matrix33_PreMultiplyBy  ( ms, symmetryParameters->M        ) ;
        Matrix33_CopyTo ( transformation3->rotation, si, NULL ) ; Matrix33_PostMultiplyBy ( si, symmetryParameters->inverseM ) ;

        /* . Loop over the elements of M. */
        for ( a = 0 ; a < 3 ; a++ )
        {
            for ( b = 0 ; b < 3 ; b++ )
            {
                /* . Rotational contribution. */
                /* . Inverse derivative. */
                Matrix33_InverseDerivative ( symmetryParameters->M, a, b, di ) ;
                Matrix33_PreMultiplyBy ( di, ms ) ;

                /* . Non-inverse derivative. */
                Matrix33_GetRow       ( si, b, dx, dy, dz ) ;
                Matrix33_IncrementRow ( di, a, dx, dy, dz ) ;

                /* . Local variables. */
                r00  = Matrix33_Item ( di, 0, 0 ) ;
                r01  = Matrix33_Item ( di, 0, 1 ) ;
                r02  = Matrix33_Item ( di, 0, 2 ) ;
                r10  = Matrix33_Item ( di, 1, 0 ) ;
                r11  = Matrix33_Item ( di, 1, 1 ) ;
                r12  = Matrix33_Item ( di, 1, 2 ) ;
                r20  = Matrix33_Item ( di, 2, 0 ) ;
                r21  = Matrix33_Item ( di, 2, 1 ) ;
                r22  = Matrix33_Item ( di, 2, 2 ) ;

                /* . Translational contribution. */
                t = Vector3_Item ( transformation3->translation, b ) ;

                /* . Loop over the coordinates/gradients. */
                for ( i = 0, sum = 0.0e+00 ; i < coordinates3->length0 ; i++ )
                {
                    Coordinates3_GetRow ( coordinates3, i,  x,  y,  z ) ;
                    Coordinates3_GetRow ( gradients3,   i, gx, gy, gz ) ;

                    /* . Rotational contribution. */
                    dx = r00 * x + r01 * y + r02 * z ;
                    dy = r10 * x + r11 * y + r12 * z ;
                    dz = r20 * x + r21 * y + r22 * z ;

                    /* . Translational contribution. */
                    switch ( a )
                    {
                        case 0: dx += t ; break ;
                        case 1: dy += t ; break ;
                        case 2: dz += t ; break ;
                    }

                    /* . dEdM. */
                    sum += dx * gx + dy * gy + dz * gz ;
                }
                /* . Set the component of dEdM. */
                Matrix33_Item ( self->dEdM, a, b ) += sum ;
            }
        }

        /* . Deallocate space. */
        Matrix33_Deallocate ( &di ) ;
        Matrix33_Deallocate ( &ms ) ;
        Matrix33_Deallocate ( &si ) ;

    }
}

