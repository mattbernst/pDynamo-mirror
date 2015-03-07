/*------------------------------------------------------------------------------
! . File      : GaussianBasisCore.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==============================================================================
! . Core procedures for Gaussian basis sets.
!=============================================================================*/

# include <math.h>
# include <stdarg.h>
# include <stdlib.h>

# include "GaussianBasisCore.h"
# include "GaussianBasisSubsidiary.h"
# include "Memory.h"
# include "RysQuadrature.h"
# include "SymmetricMatrix.h"

/* . This needs to be rationalized. */

/*# define CHECKNORMALIZATION*/
# define DIMENSIONWARNING
# ifdef CHECKNORMALIZATION
static void CheckNormalization ( GaussianBasis *self ) ;
# endif

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cartesian basis function and shell parameters.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The number of Cartesian basis functions for a shell. Equivalent to MAXCBF. */
const Integer CBFFUNCT[MAXAMP1] = { 1, 3, 6, 10, 15 } ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the full transformation for the basis (Cartesian to spherical and orthogonalizing).
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasis_CalculateFullTransformation ( GaussianBasis *self )
{
    if ( self != NULL )
    {
        auto Real               r[3] = { 0.0e+00, 0.0e+00, 0.0e+00 } ;
        auto Integer                   i, ishell, j ;
        auto Real2DArray          *c2s = NULL, *integrals = NULL, *inverse = NULL, *s2c = NULL, *transformation = NULL ;
        auto SymmetricMatrix *sintegrals = NULL, *temp = NULL ;

        /* . Deallocate the previous transformations. */
        Real2DArray_Deallocate ( &(self->c2o) ) ;
        Real2DArray_Deallocate ( &(self->o2c) ) ;

        /* . Calculate the normalization integrals. */
        if ( self->type == GaussianBasisType_Coulomb ) integrals = GaussianBasis_2Coulomb ( self, r, self, r ) ;
        else                                           integrals = GaussianBasis_2Overlap ( self, r, self, r ) ;

        /* . Create a symmetric matrix with the integrals. */
        sintegrals = SymmetricMatrix_Allocate ( self->nbasisw ) ;
        SymmetricMatrix_Set ( sintegrals, 0.0e+00 ) ;
        for ( i = 0 ; i < self->nbasisw ; i++ )
        {
            for ( j = 0 ; j < self->nbasisw ; j++ ) SymmetricMatrix_Set_Component ( sintegrals, i, j, Real2DArray_Item ( integrals, i, j ) ) ;
        }

        /* . If the basis is spherical harmonic create the appropriate transformations and then transform the integrals. */
        if ( self->QSPHERICAL )
        {
            auto Real2DArray *mc, *ms ;

            /* . Create the forwards and backwards transformations. */
            c2s = Real2DArray_Allocate ( self->nbasisw, self->nbasis, NULL ) ; Real2DArray_Set ( c2s, 0.0e+00 ) ;
            s2c = Real2DArray_Allocate ( self->nbasisw, self->nbasis, NULL ) ; Real2DArray_Set ( s2c, 0.0e+00 ) ;
            for ( ishell = 0 ; ishell < self->nshells ; ishell++ )
            {
                mc = self->shells[ishell].c2s ;
                ms = self->shells[ishell].s2c ;
/*
printf ( "XXXX: %d %d\n", ishell, self->nshells ) ;
printf ( "TRANSFORMATION MC\n" ) ;
Matrix_Print ( mc ) ;
printf ( "TRANSFORMATION MS\n" ) ;
Matrix_Print ( ms ) ;
*/
                if ( mc != NULL )
                {
                    for ( i = 0 ; i < mc->length0 ; i++ )
                    {
                        for ( j = 0 ; j < mc->length1 ; j++ ) Real2DArray_Item ( c2s, i + self->shells[ishell].nstartw, j + self->shells[ishell].nstart ) =  Real2DArray_Item ( mc, i, j ) ;
                    }
                }
                else
                {
                    for ( i = 0 ; i < self->shells[ishell].nbasisw ; i++ ) Real2DArray_Item ( c2s, i + self->shells[ishell].nstartw, i + self->shells[ishell].nstart ) = 1.0e+00 ;
                }
                if ( ms != NULL )
                {
                    for ( i = 0 ; i < ms->length0 ; i++ )
                    {
                        for ( j = 0 ; j < ms->length1 ; j++ ) Real2DArray_Item ( s2c, i + self->shells[ishell].nstartw, j + self->shells[ishell].nstart ) =  Real2DArray_Item ( ms, i, j ) ;
                    }
                }
                else
                {
                    for ( i = 0 ; i < self->shells[ishell].nbasisw ; i++ ) Real2DArray_Item ( s2c, i + self->shells[ishell].nstartw, i + self->shells[ishell].nstart ) = 1.0e+00 ;
                }
            }

            /* . Transform the integrals. */
/*
printf ( "INTEGRALS BEFORE\n" ) ;
SymmetricMatrix_Print ( sintegrals ) ;
printf ( "TRANSFORMATION\n" ) ;
Matrix_Print ( c2s ) ;
*/
            SymmetricMatrix_Transform_In_Place ( sintegrals, c2s ) ;
/*
printf ( "INTEGRALS AFTER\n" ) ;
SymmetricMatrix_Print ( sintegrals ) ;
*/
        }

        /* . Determine the orthogonalizing transformation. */
        temp = SymmetricMatrix_Clone ( sintegrals ) ;
        transformation = SymmetricMatrix_OrthogonalizingTransformation ( temp, NULL, NULL, NULL ) ;

/*
printf ( "OVERLAP INTEGRALS\n" ) ;
SymmetricMatrix_Print ( sintegrals ) ;
printf ( "TRANSFORMATION\n" ) ;
Matrix_Print ( transformation ) ;
*/

        /* . Create the inverse transformation - S * X. */
        inverse = Real2DArray_Allocate ( sintegrals->dimension, transformation->length1, NULL ) ;
        SymmetricMatrix_PostMatrixMultiply ( sintegrals, transformation, False, inverse ) ;

/*
printf ( "\n\nDIMENSIONS = %d %d %d %d %d\n\n", self->nbasisw, transformation->length0, self->nbasis, sintegrals->dimension, transformation->length1 ) ;
{
auto Real maxv1, maxv2 ;
auto Integer i, j ;
SymmetricMatrix_Transform_In_Place ( sintegrals, transformation ) ;
maxv1 = 0.0e+00 ;
maxv2 = 0.0e+00 ;
for ( i = 0 ; i < sintegrals->dimension ; i++ )
{
    for ( j = 0 ; j < i ; j++ ) maxv1 = Maximum ( maxv1, fabs ( SymmetricMatrix_Get_Component ( sintegrals, i, j ) ) ) ;
    maxv2 = Maximum ( maxv2, fabs ( SymmetricMatrix_Get_Component ( sintegrals, i, i ) - 1.0e+00 ) ) ;
}
printf ( "\n\nMaximum Deviation Off-Diagonal = %25.15f %10d\n", maxv1, self->type ) ;
printf ( "\n\nMaximum Deviation On -Diagonal = %25.15f %10d\n", maxv2, self->type ) ;
if ( ( maxv1 > 1.0e-10 ) || ( maxv2 > 1.0e-10 ) )
{
    printf ( "\nIntegral Matrix:\n" ) ;
    for ( i = 0 ; i < sintegrals->dimension ; i++ )
    {
        for ( j = 0 ; j <= i ; j++ ) printf ( "%15.10f", SymmetricMatrix_Get_Component ( sintegrals, i, j ) ) ;
        printf ( "\n" ) ;
    }
    printf ( "\n" ) ;
}
}
transformation = Real2DArray_Allocate ( sintegrals->dimension, sintegrals->dimension ) ;
Matrix_Set ( transformation, 0.0e+00 ) ;
for ( i = 0 ; i < transformation->length0 ; i++ ) Real2DArray_Item ( transformation, i, i ) = 1.0e+00 ;
*/

# ifdef DIMENSIONWARNING
if ( transformation->length0 != transformation->length1 ) printf ( "\n\nWARNING - DIMENSION CHANGE = %d TO %d\n\n", transformation->length0, transformation->length1 ) ;
# endif

        /* . Set the transformations. */
        if ( self->QSPHERICAL )
        {
            self->c2o = Real2DArray_Allocate ( self->nbasisw, transformation->length1, NULL ) ; Real2DArray_Set ( self->c2o, 0.0e+00 ) ;
            self->o2c = Real2DArray_Allocate ( self->nbasisw, transformation->length1, NULL ) ; Real2DArray_Set ( self->o2c, 0.0e+00 ) ;
            Real2DArray_MatrixMultiply ( False, False, 1.0e+00, c2s, transformation, 0.0e+00, self->c2o, NULL ) ;
            Real2DArray_MatrixMultiply ( False, False, 1.0e+00, s2c, inverse,        0.0e+00, self->o2c, NULL ) ;
        }
        else
        {
            self->c2o = Real2DArray_Clone ( transformation, NULL ) ;
            self->o2c = Real2DArray_Clone ( inverse       , NULL ) ;
        }

/*
{
auto Real maxv1 = 0.0e+00, maxv2 = 0.0e+00 ;
auto Integer i, j ;
auto Real2DArray *temp ;
temp = Real2DArray_Allocate ( self->c2o->length1, self->o2c->length1 ) ;
Matrix_MatrixMultiply ( True, False, 1.0e+00, self->c2o, self->o2c, 0.0e+00, temp ) ;
for ( i = 0 ; i < temp->length0 ; i++ )
{
    for ( j = 0 ; j < temp->length1 ; j++ )
    {
        if ( j == i ) maxv1 = Maximum ( maxv1, fabs ( Real2DArray_Item ( temp, i, j ) - 1.0e+00 ) ) ;
        else          maxv2 = Maximum ( maxv2, fabs ( Real2DArray_Item ( temp, i, j ) ) ) ;
    }
}
printf ( "\n\n Maximum deviations X = %25.15f %25.15f\n\n", maxv1, maxv2 ) ;
Matrix_Deallocate ( &temp ) ;
}
*/

        /* . Clean up. */
        Real2DArray_Deallocate ( &c2s            ) ;
        Real2DArray_Deallocate ( &integrals      ) ;
        Real2DArray_Deallocate ( &inverse        ) ;
        Real2DArray_Deallocate ( &s2c            ) ;
        Real2DArray_Deallocate ( &transformation ) ;
        SymmetricMatrix_Deallocate ( &sintegrals ) ;
        SymmetricMatrix_Deallocate ( &temp       ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Normalize a basis.
! . The normalization is performed for Cartesian functions.
! . Calls to GaussianBasis_ScaleExponents and to GaussianBasis_UnnormalizePrimitives must occur before a call to this function.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define INTEGRALTOLERANCE 1.0e-20
void GaussianBasis_Normalize ( GaussianBasis *self )
{
    if ( self != NULL )
    {
        auto Real       integral, norm, r[3] = { 0.0e+00, 0.0e+00, 0.0e+00 } ;
        auto Integer          am, i, ishell, m, n, p ;
        auto Real2DArray *integrals = NULL ;

        /* . Fill the primitive ccbf for the basis. */
        for ( ishell = 0 ; ishell < self->nshells ; ishell++ )
        {
            for ( p = 0 ; p < self->shells[ishell].nprimitives ; p++ )
            {
                for ( am = self->shells[ishell].type->angularmomentum_low, n = 0 ; am <= self->shells[ishell].type->angularmomentum_high ; am++ )
                {
                    for ( i = 0 ; i < CBFFUNCT[am] ; i++, n++ ) self->shells[ishell].primitives[p].ccbf[n] = self->shells[ishell].primitives[p].coefficients[am-self->shells[ishell].type->angularmomentum_low] ;
                }
            }
        }

        /* . Calculate the normalization integrals. */
        if ( self->type == GaussianBasisType_Coulomb ) integrals = GaussianBasis_2Coulomb ( self, r, self, r ) ;
        else                                           integrals = GaussianBasis_2Overlap ( self, r, self, r ) ;

        /* . Normalize the ccbf taking into account the values of the diagonal integrals. */
        for ( ishell = n = 0 ; ishell < self->nshells ; ishell++ )
        {
            for ( am = self->shells[ishell].type->angularmomentum_low, m = 0 ; am <= self->shells[ishell].type->angularmomentum_high ; am++ )
            {
                for ( i = 0 ; i < CBFFUNCT[am] ; i++, m++, n++ )
                {
                    integral = Real2DArray_Item ( integrals, n, n ) ;
                    if ( integral > INTEGRALTOLERANCE )
                    {
                        norm = 1.0e+00 / sqrt ( integral ) ;
                        for ( p = 0 ; p < self->shells[ishell].nprimitives ; p++ ) self->shells[ishell].primitives[p].ccbf[m] *= norm ;
                    }
                }
            }
        }

        /* . Free the integrals. */
        Real2DArray_Deallocate ( &integrals ) ;

        /* . Calculate the full transformation for the basis. */
        GaussianBasis_CalculateFullTransformation ( self ) ;

# ifdef CHECKNORMALIZATION
CheckNormalization ( self ) ;
# endif

    }
}
# undef INTEGRALTOLERANCE

/*----------------------------------------------------------------------------------------------------------------------------------
! . Scale the exponents of the shells of a Gaussian basis set.
! . Scaling is done from exponent0.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasis_ScaleExponents ( GaussianBasis *self, ... )
{
    if ( self != NULL )
    {
        if ( self->nshells > 0 )
        {
            auto Real  zeta, zeta2 ;
            auto Integer     i, p        ;
            auto va_list argp        ;
            /* . Loop over the shells of the basis. */
            va_start ( argp, self ) ;
            for ( i = 0 ; i < self->nshells ; i++ )
            {
                zeta  = va_arg ( argp, Real ) ;
                zeta2 = zeta * zeta ;
                for ( p = 0 ; p < self->shells[i].nprimitives ; p++ ) self->shells[i].primitives[p].exponent = self->shells[i].primitives[p].exponent0 * zeta2 ;
            }
            va_end ( argp ) ;
        }
        /* . Renormalize the basis. */
        GaussianBasis_Normalize ( self ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Unnormalize the primitives of a basis.
! . Only done if QNORMALIZEDPRIMITIVES is True and done from coefficients0
! . using exponents0.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasis_UnnormalizePrimitives ( GaussianBasis *self )
{
    if ( ( self != NULL ) && ( self->QNORMALIZEDPRIMITIVES ) )
    {
        auto Integer     am, i, ishell, p ;
        auto Real  ex, fac ;
        /* . Loop over the coefficients of each angular momentum for each shell. */
        for ( ishell = 0 ; ishell < self->nshells ; ishell++ )
        {
            for ( am = self->shells[ishell].type->angularmomentum_low ; am <= self->shells[ishell].type->angularmomentum_high ; am++ )
            {
                for ( p = 0 ; p < self->shells[ishell].nprimitives ; p++ )
                {
                    ex  = 2.0e+00 * self->shells[ishell].primitives[p].exponent0 ;
                    fac = PI32 / ( ex * sqrt ( ex ) ) ;
                    for ( i = 1 ; i <= am ; i++ ) fac *= ( Real ) ( 2*i - 1 ) / ( 2.0e+00 * ex ) ;
                    self->shells[ishell].primitives[p].coefficients [am-self->shells[ishell].type->angularmomentum_low] =
                    self->shells[ishell].primitives[p].coefficients0[am-self->shells[ishell].type->angularmomentum_low] / sqrt ( fac ) ;
                }
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate Coulomb integrals between two bases.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real2DArray *GaussianBasis_2Coulomb ( const GaussianBasis *ibasis, const Real *ri, const GaussianBasis *jbasis, const Real *rj )
{
    auto Real2DArray *integrals = NULL ;
    if ( ( ibasis != NULL ) && ( jbasis != NULL ) )
    {
        auto Real             aandb, ab, ai, aj, bp01, b00, b10, c1x, c1y, c1z, c3x, c3y, c3z, dfi, dfij, fac, fac2, f00, rho, rij2, ti, tij, u2, xc00, xcp00, xij, yc00, ycp00, yij, zc00, zcp00, zij ;
        auto Real             g[MAXCBF*MAXCBF], xint[MAXAMP1*MAXAMP1*MAXRYS], yint[MAXAMP1*MAXAMP1*MAXRYS], zint[MAXAMP1*MAXAMP1*MAXRYS] ;
        auto Integer                i, iammax, icbfind, ip, ishell, ix, iy, iz, j, jammax, jcbfind, jdim, jdimm, jp, jshell, jxix, jyiy, jziz, m, n, ncfunci, ncfuncj, nroots ;
        auto RysQuadrature roots ;

        /* . Allocate space. */
        integrals = Real2DArray_Allocate ( ibasis->nbasisw, jbasis->nbasisw, NULL ) ;
        Real2DArray_Set ( integrals, 0.0e+00 ) ;

        /* . Calculate some distance factors. */
	xij  = ri[0] - rj[0] ;
	yij  = ri[1] - rj[1] ;
	zij  = ri[2] - rj[2] ;
	rij2 = xij * xij + yij * yij + zij * zij ;

        /* . Outer loop over shells. */
        for ( ishell = 0 ; ishell < ibasis->nshells ; ishell++ )
        {

            /* . Get information about the shell. */
            iammax  = ibasis->shells[ishell].type->angularmomentum_high ;
            icbfind = ibasis->shells[ishell].type->cbfindex ;
            ncfunci = ibasis->shells[ishell].type->ncbf     ;

            /* . Inner loop over shells. */
            for ( jshell = 0 ; jshell < jbasis->nshells ; jshell++ )
            {

                /* . Get information about the shell. */
                jammax  = jbasis->shells[jshell].type->angularmomentum_high ;
                jdim    = jammax + 1 ;
                jdimm   = ( iammax + 1 ) * ( jammax + 1 ) ;
                jcbfind = jbasis->shells[jshell].type->cbfindex ;
                ncfuncj = jbasis->shells[jshell].type->ncbf     ;

                /* . Get the number of roots. */
                nroots = ( iammax + jammax ) / 2 + 1 ;

                /* . Initialize the integral blocks. */
                for ( i = 0 ; i < ( ncfunci * ncfuncj ) ; i++ ) g[i] = 0.0e+00 ;

                /* . Outer loop over primitives. */
                for ( ip = 0 ; ip < ibasis->shells[ishell].nprimitives ; ip++ )
                {
                    /* . Get some information for the primitive. */
	            ai  = ibasis->shells[ishell].primitives[ip].exponent ;
                    dfi = PI252 / ai ;

                    /* . Inner loop over primitives. */
                    for ( jp = 0 ; jp < jbasis->shells[jshell].nprimitives ; jp++ )
                    {
                        /* . Get some information for the primitive. */
	                aj = jbasis->shells[jshell].primitives[jp].exponent ;

                        /* . Calculate some factors. */
                        ab    = ai * aj ;
                        aandb = ai + aj ;
                        rho   = ab / aandb ;
                        dfij  = dfi / ( aj * sqrt ( aandb ) ) ;

                        /* . Calculate some displacements. */
                        c1x  = ai * ( ri[0] - rj[0] ) ;
                        c1y  = ai * ( ri[1] - rj[1] ) ;
                        c1z  = ai * ( ri[2] - rj[2] ) ;
                        c3x  = aj * ( rj[0] - ri[0] ) ;
                        c3y  = aj * ( rj[1] - ri[1] ) ;
                        c3z  = aj * ( rj[2] - ri[2] ) ;

                        /* . Calculate the rys polynomial roots. */
                        RysQuadrature_Roots ( &roots, nroots, ( rho * rij2 ) ) ;

                        /* . Loop over the roots and construct the subsidiary integrals. */
                        for ( m = 0 ; m < nroots ; m++ )
                        {
                            u2    = roots.roots[m] * rho ;
                            f00   = roots.weights[m] ;
                            fac   = 1.0e+00 / ( ab + u2 * aandb ) ;
                            fac2  = 0.5e+00 * fac ;
                            bp01  = ( ai + u2 ) * fac2 ;
                            b00   =        u2   * fac2 ;
                            b10   = ( aj + u2 ) * fac2 ;
                            xcp00 = u2 * c1x * fac ;
                            xc00  = u2 * c3x * fac ;
                            ycp00 = u2 * c1y * fac ;
                            yc00  = u2 * c3y * fac ;
                            zcp00 = u2 * c1z * fac ;
                            zc00  = u2 * c3z * fac ;
                            Subsidiary_Integral_Nuclear2C ( iammax, jammax, b00, b10, bp01, f00, xc00, xcp00, yc00, ycp00, zc00, zcp00, jdim, &xint[m*jdimm], &yint[m*jdimm], &zint[m*jdimm] ) ;
                        }

                        /* . Assemble the integrals. */
                        for ( i = 0, n = 0 ; i < ncfunci ; i++ )
                        {
   	                    ix = CBFPOWX[i+icbfind] * jdim ;
	                    iy = CBFPOWY[i+icbfind] * jdim ;
	                    iz = CBFPOWZ[i+icbfind] * jdim ;
                            ti = dfij * ibasis->shells[ishell].primitives[ip].ccbf[i] ;
                            for ( j = 0 ; j < ncfuncj ; j++ )
                            {
	                        jxix = CBFPOWX[j+jcbfind] + ix ;
	                        jyiy = CBFPOWY[j+jcbfind] + iy ;
	                        jziz = CBFPOWZ[j+jcbfind] + iz ;
                                for ( m = 0, fac = 0.0e+00 ; m < nroots ; m++ ) fac += xint[jxix+m*jdimm] * yint[jyiy+m*jdimm] * zint[jziz+m*jdimm] ;
                                tij = ti * jbasis->shells[jshell].primitives[jp].ccbf[j] ;
                                g[n] += tij * fac ;
                                n++ ;
                            }
                        }
                    }
                }

                /* . Transform the integrals. */
	        /* if ( qcAtoms->QTOSPHERICAL ) Integral_Block_Transform_M ( g, iangmom, ncfunci, nfunci, jangmom, ncfuncj, nfuncj ) ; */

                /* . Put the integrals in the proper place. */
                for ( i = 0, n = 0 ; i < ibasis->shells[ishell].nbasisw ; i++ )
                {
                    for ( j = 0 ; j < jbasis->shells[jshell].nbasisw ; j++, n++ ) Real2DArray_Item ( integrals, i+ibasis->shells[ishell].nstartw, j+jbasis->shells[jshell].nstartw ) = g[n] ;
                }
            }
        }
    }
    return integrals ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the overlap between two bases.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real2DArray *GaussianBasis_2Overlap ( const GaussianBasis *ibasis, const Real *ri, const GaussianBasis *jbasis, const Real *rj )
{
    auto Real2DArray *overlap = NULL ;
    if ( ( ibasis != NULL ) && ( jbasis != NULL ) )
    {
        auto Real  aa, aainv, ai, aj, arri, expfac, fac, rij2, xij2, yij2, zij2 ;
        auto Real  ar[3], ari[3], s[MAXCBF*MAXCBF], xo[MAXAMP1*MAXAMP1], yo[MAXAMP1*MAXAMP1], zo[MAXAMP1*MAXAMP1] ;
        auto Integer     i, iammax, icbfind, ip, ishell, ix, iy, iz, j, jammax, jcbfind, jdim, jp, jshell, jxix, jyiy, jziz, n, ncfunci, ncfuncj ;

        /* . Allocate space. */
        overlap = Real2DArray_Allocate ( ibasis->nbasisw, jbasis->nbasisw, NULL ) ;
        Real2DArray_Set ( overlap, 0.0e+00 ) ;

        /* . Get the distance squared between atom centers. */
        xij2 = ri[0] - rj[0] ;
        yij2 = ri[1] - rj[1] ;
        zij2 = ri[2] - rj[2] ;
        rij2 = xij2 * xij2 + yij2 * yij2 + zij2 * zij2 ;

        /* . Outer loop over shells. */
        for ( ishell = 0 ; ishell < ibasis->nshells ; ishell++ )
        {

           /* . Get information about the shell. */
           iammax  = ibasis->shells[ishell].type->angularmomentum_high ;
           icbfind = ibasis->shells[ishell].type->cbfindex ;
           ncfunci = ibasis->shells[ishell].type->ncbf     ;

           /* . Inner loop over shells. */
           for ( jshell = 0 ; jshell < jbasis->nshells ; jshell++ )
           {

              /* . Get information about the shell. */
              jammax  = jbasis->shells[jshell].type->angularmomentum_high ;
              jdim    = jammax + 1 ;
              jcbfind = jbasis->shells[jshell].type->cbfindex ;
              ncfuncj = jbasis->shells[jshell].type->ncbf     ;

              /* . Initialize the integral block. */
              for ( i = 0 ; i < ( ncfunci * ncfuncj ) ; i++ ) s[i] = 0.0e+00 ;

              /* . Outer loop over primitives. */
              for ( ip = 0 ; ip < ibasis->shells[ishell].nprimitives ; ip++ )
              {
                 /* . Get some information for the primitive. */
	         ai   = ibasis->shells[ishell].primitives[ip].exponent ;
	         arri = ai * rij2 ;
	         for ( i = 0 ; i < 3 ; i++ ) ari[i] = ai * ri[i] ;

                 /* . Inner loop over primitives. */
                 for ( jp = 0 ; jp < jbasis->shells[jshell].nprimitives ; jp++ )
                 {
                    /* . Get some information for the primitive. */
	            aj    = jbasis->shells[jshell].primitives[jp].exponent ;
	            aa    = ai + aj ;
	            aainv = 1.0e+00 / aa ;
	            fac   = aj * arri * aainv ;
	            if ( fac > PRIMITIVE_OVERLAP_TOLERANCE ) continue ;
                    expfac = exp ( - fac ) ;
	            for ( i = 0 ; i < 3 ; i++ ) ar[i] = ( ari[i] + aj * rj[i] ) * aainv ;

                    /* . Calculate the overlap integrals. */
                    Subsidiary_Integral_Overlap2 ( xo, yo, zo, aa, ar, ri, rj, iammax, jammax ) ;

                    /* . Add in the contributions to the full integrals. */
                    n = 0 ;
                    for ( i = 0 ; i < ncfunci ; i++ )
                    {
   	               ix = CBFPOWX[i+icbfind] * jdim ;
	               iy = CBFPOWY[i+icbfind] * jdim ;
	               iz = CBFPOWZ[i+icbfind] * jdim ;
                       for ( j = 0 ; j < ncfuncj ; j++ )
                       {
	                  jxix = CBFPOWX[j+jcbfind] + ix ;
	                  jyiy = CBFPOWY[j+jcbfind] + iy ;
	                  jziz = CBFPOWZ[j+jcbfind] + iz ;
    		          s[n] += expfac * ibasis->shells[ishell].primitives[ip].ccbf[i] *
                                           jbasis->shells[jshell].primitives[jp].ccbf[j] * xo[jxix] * yo[jyiy] * zo[jziz] ;
                          n++ ;
                       }
                    }
                 }
              }

              /* . Transform to spherical harmonics. */
  /*            if ( QTOSPHERICAL ) Integral_Block_Transform_M ( s, iangmom, ncfunci, nfunci, jangmom, ncfuncj, nfuncj ) ; */

              /* . Put the integrals in the correct place. */
              for ( i = 0, n = 0 ; i < ibasis->shells[ishell].nbasisw ; i++ )
              {
                  for ( j = 0 ; j < jbasis->shells[jshell].nbasisw ; j++, n++ ) Real2DArray_Item ( overlap, i+ibasis->shells[ishell].nstartw, j+jbasis->shells[jshell].nstartw ) = s[n] ;
              }
           }
        }
    }
    return overlap ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the derivative overlap between two bases.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasis_2OverlapD ( const GaussianBasis *ibasis, const Real *ri, const GaussianBasis *jbasis, const Real *rj, Real2DArray **sijx, Real2DArray **sijy, Real2DArray **sijz )
{
   Real  aa, aainv, ai, aj, arri, denfac, expfac, fac, rij2, xij2, yij2, zij2 ;
   Real  ar[3], ari[3], sx[MAXCBF*MAXCBF],       sy[MAXCBF*MAXCBF],       sz[MAXCBF*MAXCBF],
                          xd[MAXAMP1*MAXAMP1],     yd[MAXAMP1*MAXAMP1],     zd[MAXAMP1*MAXAMP1],
                          xo[MAXAMP1*(MAXAMP1+1)], yo[MAXAMP1*(MAXAMP1+1)], zo[MAXAMP1*(MAXAMP1+1)] ;
   Integer     i, iammax, icbfind, ip, ishell, ix, iy, iz, j, jammax, jcbfind, jdim, jp, jshell, jxix, jyiy, jziz, n, ncfunci, ncfuncj ;
   if ( ( ibasis == NULL ) || ( jbasis == NULL ) )
   {
      (*sijx) = NULL ; (*sijy) = NULL ; (*sijz) = NULL ;
   }
   else
   {
      /* . Allocate space. */
      (*sijx) = Real2DArray_Allocate ( ibasis->nbasisw, jbasis->nbasisw, NULL ) ; Real2DArray_Set ( (*sijx), 0.0e+00 ) ;
      (*sijy) = Real2DArray_Allocate ( ibasis->nbasisw, jbasis->nbasisw, NULL ) ; Real2DArray_Set ( (*sijy), 0.0e+00 ) ;
      (*sijz) = Real2DArray_Allocate ( ibasis->nbasisw, jbasis->nbasisw, NULL ) ; Real2DArray_Set ( (*sijz), 0.0e+00 ) ;

      /* . Get the distance squared between atom centers. */
      xij2 = ri[0] - rj[0] ;
      yij2 = ri[1] - rj[1] ;
      zij2 = ri[2] - rj[2] ;
      rij2 = xij2 * xij2 + yij2 * yij2 + zij2 * zij2 ;

      /* . Outer loop over shells. */
      for ( ishell = 0 ; ishell < ibasis->nshells ; ishell++ )
      {

         /* . Get information about the shell. */
         iammax  = ibasis->shells[ishell].type->angularmomentum_high ;
         icbfind = ibasis->shells[ishell].type->cbfindex ;
         ncfunci = ibasis->shells[ishell].type->ncbf     ;

         /* . Inner loop over shells. */
         for ( jshell = 0 ; jshell < jbasis->nshells ; jshell++ )
         {

            /* . Get information about the shell. */
            jammax  = jbasis->shells[jshell].type->angularmomentum_high ;
            jdim    = jammax + 1 ;
            jcbfind = jbasis->shells[jshell].type->cbfindex ;
            ncfuncj = jbasis->shells[jshell].type->ncbf     ;

            /* . Initialize the integral block. */
            for ( i = 0 ; i < ( ncfunci * ncfuncj ) ; i++ ) { sx[i] = 0.0e+00 ; sy[i] = 0.0e+00 ; sz[i] = 0.0e+00 ; }

            /* . Outer loop over primitives. */
            for ( ip = 0 ; ip < ibasis->shells[ishell].nprimitives ; ip++ )
            {
               /* . Get some information for the primitive. */
	       ai   = ibasis->shells[ishell].primitives[ip].exponent ;
	       arri = ai * rij2 ;
	       for ( i = 0 ; i < 3 ; i++ ) ari[i] = ai * ri[i] ;

               /* . Inner loop over primitives. */
               for ( jp = 0 ; jp < jbasis->shells[jshell].nprimitives ; jp++ )
               {
                  /* . Get some information for the primitive. */
	          aj    = jbasis->shells[jshell].primitives[jp].exponent ;
	          aa    = ai + aj ;
	          aainv = 1.0e+00 / aa ;
	          fac   = aj * arri * aainv ;
	          if ( fac > PRIMITIVE_OVERLAP_TOLERANCE ) continue ;
                  expfac = exp ( - fac ) ;
	          for ( i = 0 ; i < 3 ; i++ ) ar[i] = ( ari[i] + aj * rj[i] ) * aainv ;

                  /* . Calculate the overlap integrals and derivatives. */
                  Subsidiary_Integral_Overlap2    ( xo, yo, zo, aa, ar, ri, rj, iammax+1, jammax ) ;
                  Subsidiary_Integral_Derivative2 ( xo, yo, zo, ai, iammax, jammax, jdim, xd, yd, zd ) ;

                  /* . Add in the contributions to the full integrals. */
                  for ( i = 0, n = 0 ; i < ncfunci ; i++ )
                  {
   	             ix = CBFPOWX[i+icbfind] * jdim ;
	             iy = CBFPOWY[i+icbfind] * jdim ;
	             iz = CBFPOWZ[i+icbfind] * jdim ;
                     for ( j = 0 ; j < ncfuncj ; j++, n++ )
                     {
	                jxix = CBFPOWX[j+jcbfind] + ix ;
	                jyiy = CBFPOWY[j+jcbfind] + iy ;
	                jziz = CBFPOWZ[j+jcbfind] + iz ;
                        denfac = expfac * ibasis->shells[ishell].primitives[ip].ccbf[i] * jbasis->shells[jshell].primitives[jp].ccbf[j] ;
    		        sx[n] += denfac * xd[jxix] * yo[jyiy] * zo[jziz] ;
    		        sy[n] += denfac * xo[jxix] * yd[jyiy] * zo[jziz] ;
    		        sz[n] += denfac * xo[jxix] * yo[jyiy] * zd[jziz] ;
                     }
                  }
               }
            }

            /* . Transform to spherical harmonics. */
/*            if ( QTOSPHERICAL )
            {
                Integral_Block_Transform_M ( sx, iangmom, ncfunci, nfunci, jangmom, ncfuncj, nfuncj ) ;
                Integral_Block_Transform_M ( sy, iangmom, ncfunci, nfunci, jangmom, ncfuncj, nfuncj ) ;
                Integral_Block_Transform_M ( sz, iangmom, ncfunci, nfunci, jangmom, ncfuncj, nfuncj ) ;
            }
*/

            /* . Put the integrals in the correct place. */
            for ( i = 0, n = 0 ; i < ibasis->shells[ishell].nbasisw ; i++ )
            {
                for ( j = 0 ; j < jbasis->shells[jshell].nbasisw ; j++, n++ )
                {
                    Real2DArray_Item ( (*sijx), i+ibasis->shells[ishell].nstartw, j+jbasis->shells[jshell].nstartw ) = sx[n] ;
                    Real2DArray_Item ( (*sijy), i+ibasis->shells[ishell].nstartw, j+jbasis->shells[jshell].nstartw ) = sy[n] ;
                    Real2DArray_Item ( (*sijz), i+ibasis->shells[ishell].nstartw, j+jbasis->shells[jshell].nstartw ) = sz[n] ;
                }
            }
         }
      }
   }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the overlap and its derivatives between two bases.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasis_2OverlapFD ( const GaussianBasis  *ibasis ,
                                const Real           *ri     ,
                                const GaussianBasis  *jbasis ,
                                const Real           *rj     ,
                                      Real2DArray   **sij    ,
                                      Real2DArray   **sijx   ,
                                      Real2DArray   **sijy   ,
                                      Real2DArray   **sijz   )
{
   Integer i, iammax, icbfind, ip, ishell, ix, iy, iz, j, jammax, jcbfind, jdim, jp, jshell, jxix, jyiy, jziz, n, ncfunci, ncfuncj ;
   Real    aa, aainv, ai, aj, arri, denfac, expfac, fac, rij2, xij2, yij2, zij2 ;
   Real    ar[3], ari[3], s[MAXCBF*MAXCBF]       , sx[MAXCBF*MAXCBF]      , sy[MAXCBF*MAXCBF]       , sz[MAXCBF*MAXCBF] ,
                          xd[MAXAMP1*MAXAMP1]    , yd[MAXAMP1*MAXAMP1]    , zd[MAXAMP1*MAXAMP1]     ,
                          xo[MAXAMP1*(MAXAMP1+1)], yo[MAXAMP1*(MAXAMP1+1)], zo[MAXAMP1*(MAXAMP1+1)] ;
   if ( ( ibasis == NULL ) || ( jbasis == NULL ) )
   {
      (*sij) = NULL ; (*sijx) = NULL ; (*sijy) = NULL ; (*sijz) = NULL ;
   }
   else
   {
      /* . Allocate space. */
      (*sij ) = Real2DArray_Allocate ( ibasis->nbasisw, jbasis->nbasisw, NULL ) ; Real2DArray_Set ( (*sij ), 0.0e+00 ) ;
      (*sijx) = Real2DArray_Allocate ( ibasis->nbasisw, jbasis->nbasisw, NULL ) ; Real2DArray_Set ( (*sijx), 0.0e+00 ) ;
      (*sijy) = Real2DArray_Allocate ( ibasis->nbasisw, jbasis->nbasisw, NULL ) ; Real2DArray_Set ( (*sijy), 0.0e+00 ) ;
      (*sijz) = Real2DArray_Allocate ( ibasis->nbasisw, jbasis->nbasisw, NULL ) ; Real2DArray_Set ( (*sijz), 0.0e+00 ) ;

      /* . Get the distance squared between atom centers. */
      xij2 = ri[0] - rj[0] ;
      yij2 = ri[1] - rj[1] ;
      zij2 = ri[2] - rj[2] ;
      rij2 = xij2 * xij2 + yij2 * yij2 + zij2 * zij2 ;

      /* . Outer loop over shells. */
      for ( ishell = 0 ; ishell < ibasis->nshells ; ishell++ )
      {

         /* . Get information about the shell. */
         iammax  = ibasis->shells[ishell].type->angularmomentum_high ;
         icbfind = ibasis->shells[ishell].type->cbfindex ;
         ncfunci = ibasis->shells[ishell].type->ncbf     ;

         /* . Inner loop over shells. */
         for ( jshell = 0 ; jshell < jbasis->nshells ; jshell++ )
         {

            /* . Get information about the shell. */
            jammax  = jbasis->shells[jshell].type->angularmomentum_high ;
            jdim    = jammax + 1 ;
            jcbfind = jbasis->shells[jshell].type->cbfindex ;
            ncfuncj = jbasis->shells[jshell].type->ncbf     ;

            /* . Initialize the integral block. */
            for ( i = 0 ; i < ( ncfunci * ncfuncj ) ; i++ ) { s[i] = 0.0e+00 ; sx[i] = 0.0e+00 ; sy[i] = 0.0e+00 ; sz[i] = 0.0e+00 ; }

            /* . Outer loop over primitives. */
            for ( ip = 0 ; ip < ibasis->shells[ishell].nprimitives ; ip++ )
            {
               /* . Get some information for the primitive. */
	       ai   = ibasis->shells[ishell].primitives[ip].exponent ;
	       arri = ai * rij2 ;
	       for ( i = 0 ; i < 3 ; i++ ) ari[i] = ai * ri[i] ;

               /* . Inner loop over primitives. */
               for ( jp = 0 ; jp < jbasis->shells[jshell].nprimitives ; jp++ )
               {
                  /* . Get some information for the primitive. */
	          aj    = jbasis->shells[jshell].primitives[jp].exponent ;
	          aa    = ai + aj ;
	          aainv = 1.0e+00 / aa ;
	          fac   = aj * arri * aainv ;
	          if ( fac > PRIMITIVE_OVERLAP_TOLERANCE ) continue ;
                  expfac = exp ( - fac ) ;
	          for ( i = 0 ; i < 3 ; i++ ) ar[i] = ( ari[i] + aj * rj[i] ) * aainv ;

                  /* . Calculate the overlap integrals and derivatives. */
                  Subsidiary_Integral_Overlap2    ( xo, yo, zo, aa, ar, ri, rj, iammax+1, jammax ) ;
                  Subsidiary_Integral_Derivative2 ( xo, yo, zo, ai, iammax, jammax, jdim, xd, yd, zd ) ;

                  /* . Add in the contributions to the full integrals. */
                  for ( i = 0, n = 0 ; i < ncfunci ; i++ )
                  {
   	             ix = CBFPOWX[i+icbfind] * jdim ;
	             iy = CBFPOWY[i+icbfind] * jdim ;
	             iz = CBFPOWZ[i+icbfind] * jdim ;
                     for ( j = 0 ; j < ncfuncj ; j++, n++ )
                     {
	                jxix = CBFPOWX[j+jcbfind] + ix ;
	                jyiy = CBFPOWY[j+jcbfind] + iy ;
	                jziz = CBFPOWZ[j+jcbfind] + iz ;
                        denfac = expfac * ibasis->shells[ishell].primitives[ip].ccbf[i] * jbasis->shells[jshell].primitives[jp].ccbf[j] ;
                        s [n] += denfac * xo[jxix] * yo[jyiy] * zo[jziz] ;
    		        sx[n] += denfac * xd[jxix] * yo[jyiy] * zo[jziz] ;
    		        sy[n] += denfac * xo[jxix] * yd[jyiy] * zo[jziz] ;
    		        sz[n] += denfac * xo[jxix] * yo[jyiy] * zd[jziz] ;
                     }
                  }
               }
            }

            /* . Transform to spherical harmonics. */
/*            if ( QTOSPHERICAL )
            {
                Integral_Block_Transform_M ( sx, iangmom, ncfunci, nfunci, jangmom, ncfuncj, nfuncj ) ;
                Integral_Block_Transform_M ( sy, iangmom, ncfunci, nfunci, jangmom, ncfuncj, nfuncj ) ;
                Integral_Block_Transform_M ( sz, iangmom, ncfunci, nfunci, jangmom, ncfuncj, nfuncj ) ;
            }
*/

            /* . Put the integrals in the correct place. */
            for ( i = 0, n = 0 ; i < ibasis->shells[ishell].nbasisw ; i++ )
            {
                for ( j = 0 ; j < jbasis->shells[jshell].nbasisw ; j++, n++ )
                {
                    Real2DArray_Item ( (*sij ), i+ibasis->shells[ishell].nstartw, j+jbasis->shells[jshell].nstartw ) = s [n] ;
                    Real2DArray_Item ( (*sijx), i+ibasis->shells[ishell].nstartw, j+jbasis->shells[jshell].nstartw ) = sx[n] ;
                    Real2DArray_Item ( (*sijy), i+ibasis->shells[ishell].nstartw, j+jbasis->shells[jshell].nstartw ) = sy[n] ;
                    Real2DArray_Item ( (*sijz), i+ibasis->shells[ishell].nstartw, j+jbasis->shells[jshell].nstartw ) = sz[n] ;
                }
            }
         }
      }
   }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check normalization.
!---------------------------------------------------------------------------------------------------------------------------------*/
# ifdef CHECKNORMALIZATION
static void CheckNormalization ( GaussianBasis *self )
{
    if ( self != NULL )
    {
        auto Real  maxv, r[3] = { 0.0e+00, 0.0e+00, 0.0e+00 } ;
        auto Integer     i, j, n ;
        auto Real2DArray          *integrals  = NULL ;
        auto SymmetricMatrix *sintegrals = NULL ;

        /* . Calculate the normalization integrals. */
        if ( self->type == GaussianBasisType_Coulomb ) integrals = GaussianBasis_2Coulomb ( self, r, self, r ) ;
        else                                           integrals = GaussianBasis_2Overlap ( self, r, self, r ) ;

        /* . Create a symmetric matrix with the integrals. */
        sintegrals = SymmetricMatrix_Allocate ( self->nbasisw ) ;
        SymmetricMatrix_Set ( sintegrals, 0.0e+00 ) ;
        for ( i = n = 0 ; i < self->nbasisw ; i++ )
        {
            for ( j = 0 ; j < self->nbasisw ; j++, n++ ) SymmetricMatrix_Set_Component ( sintegrals, i, j, Real2DArray_Item ( integrals, i, j ) ) ;
        }

        /* . Transform the integrals. */
        SymmetricMatrix_Transform_In_Place ( sintegrals, self->c2o ) ;

        /* . Find the maximum deviation of an element from the identity. */
        maxv = 0.0e+00 ;
        for ( i = 0 ; i < sintegrals->dimension ; i++ )
        {
            for ( j = 0 ; j < i ; j++ ) maxv = Maximum ( maxv, fabs ( SymmetricMatrix_Get_Component ( sintegrals, i, j ) ) ) ;
            maxv = Maximum ( maxv, fabs ( SymmetricMatrix_Get_Component ( sintegrals, i, i ) - 1.0e+00 ) ) ;
        }
        printf ( "\n\nMaximum Deviation = %25.15f %10d\n", maxv, self->type ) ;
        if ( maxv > 1.0e-10 )
        {
            printf ( "\nSymmetric Matrix:\n" ) ;
            for ( i = 0 ; i < sintegrals->dimension ; i++ )
            {
                for ( j = 0 ; j <= i ; j++ ) printf ( "%15.10f", SymmetricMatrix_Get_Component ( sintegrals, i, j ) ) ;
                printf ( "\n" ) ;
            }
            printf ( "\n" ) ;
        }

        /* . Free the integrals. */
        Real2DArray_Deallocate ( &integrals ) ;
        SymmetricMatrix_Deallocate ( &sintegrals ) ;
    }
}
# endif
