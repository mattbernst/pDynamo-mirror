/*------------------------------------------------------------------------------
! . File      : GaussianBasisIntegrals.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==============================================================================
! . Procedures for calculating integrals over Gaussians.
!=============================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "GaussianBasisIntegrals.h"
# include "GaussianBasisSubsidiary.h"
# include "Memory.h"
# include "RysQuadrature.h"

/* # define PRINTSMALLEIGENVALUES */
/* # define PRINTDIAGONALELEMENTS */

/*------------------------------------------------------------------------------
! . Dipole integrals.
!-----------------------------------------------------------------------------*/
void GaussianBasisIntegrals_Dipole ( const QCAtomContainer *qcAtoms, const QCParameter *qcParameters, Coordinates3 *qccoordinates3, const Real *center,
                                                                                SymmetricMatrix **dipx, SymmetricMatrix **dipy , SymmetricMatrix **dipz )
{
    if ( ( qcAtoms != NULL ) && ( qcParameters != NULL ) && ( qccoordinates3 != NULL ) && ( center != NULL ) )
    {
        auto Boolean                QDIAGONAL, QIQMJQM ;
        auto Real              aa, aainv, ai, aj, arri, expfac, fac, rij2, ti, tij, xij, yij, zij ;
        auto Real              ar[3], ari[3], sx[MAXCBF*MAXCBF],   sy[MAXCBF*MAXCBF],   sz[MAXCBF*MAXCBF],
                                                xo[MAXAMP1*MAXAMP1], yo[MAXAMP1*MAXAMP1], zo[MAXAMP1*MAXAMP1],
                                                xd[MAXAMP1*MAXAMP1], yd[MAXAMP1*MAXAMP1], zd[MAXAMP1*MAXAMP1] ;
        auto Real             *ri, *rj ;
        auto Integer                 i, iammax, ic, icbfind, iqm, ip, ishell, istart, ix, iy, iz,
                                 j, jammax, jc, jcbfind, jdim, jp, jqm, jshell, jstart, jupper, jxix, jyiy, jziz, n, ncfunci, ncfuncj ;
        auto GaussianBasis *ibasis, *jbasis ;

        /* . Initialization. */
        (*dipx)    = SymmetricMatrix_Allocate ( qcAtoms->nobasisw ) ;
        (*dipy)    = SymmetricMatrix_Allocate ( qcAtoms->nobasisw ) ;
        (*dipz)    = SymmetricMatrix_Allocate ( qcAtoms->nobasisw ) ;
        SymmetricMatrix_Set_Zero ( (*dipx) ) ;
        SymmetricMatrix_Set_Zero ( (*dipy) ) ;
        SymmetricMatrix_Set_Zero ( (*dipz) ) ;

        /* . Outer loop over centers. */
        for ( iqm = 0 ; iqm < qcAtoms->natoms ; iqm++ )
        {

            /* . Get data for the center. */
            ic     = qcAtoms->data[iqm].center ;
            ibasis = qcParameters->centers[ic].orbitalbasis ;
            istart = qcAtoms->data[iqm].ostartw ;
            ri     = Coordinates3_RowPointer ( qccoordinates3, iqm ) ;

            /* . Inner loop over centers. */
            for ( jqm = 0 ; jqm <= iqm ; jqm++ )
            {

                /* . Get data for the center. */
                jc     = qcAtoms->data[jqm].center ;
                jbasis = qcParameters->centers[jc].orbitalbasis ;
                jstart = qcAtoms->data[jqm].ostartw ;
                rj     = Coordinates3_RowPointer ( qccoordinates3, jqm ) ;

                /* . Set the diagonal atom flag. */
                QIQMJQM = ( iqm == jqm ) ;

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

                    /* . Set the upper limit for the JSHELL loops. */
                    if ( QIQMJQM ) jupper = ishell + 1 ;
                    else           jupper = jbasis->nshells ;

                    /* . Inner loop over shells. */
                    for ( jshell = 0 ; jshell < jupper ; jshell++ )
                    {

                        /* . Get information about the shell. */
                        jammax  = jbasis->shells[jshell].type->angularmomentum_high ;
                        jdim    = jammax + 1 ;
                        jcbfind = jbasis->shells[jshell].type->cbfindex ;
                        ncfuncj = jbasis->shells[jshell].type->ncbf     ;

                        /* . Set the diagonal block flag. */
                        QDIAGONAL = QIQMJQM && ( ishell == jshell ) ;

                        /* . Initialize the integral blocks. */
                        for ( i = 0 ; i < ( ncfunci * ncfuncj ) ; i++ )
                        {
                            sx[i] = 0.0e+00 ;
                            sy[i] = 0.0e+00 ;
                            sz[i] = 0.0e+00 ;
                        }

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

                                /* . Calculate the subsidiary integrals. */
                                Subsidiary_Integral_Overlap2 ( xo, yo, zo, aa, ar, ri, rj,         iammax, jammax ) ;
                                Subsidiary_Integral_Dipole   ( xd, yd, zd, aa, ar, ri, rj, center, iammax, jammax ) ;

                                /* . Add in the contributions to the full integrals. */
                                for ( i = 0, n = 0 ; i < ncfunci ; i++ )
                                {
   	                            ix = CBFPOWX[i+icbfind] * jdim ;
	                            iy = CBFPOWY[i+icbfind] * jdim ;
	                            iz = CBFPOWZ[i+icbfind] * jdim ;
                                    ti = expfac * ibasis->shells[ishell].primitives[ip].ccbf[i] ;
                                    for ( j = 0 ; j < ncfuncj ; j++, n++ )
                                    {
	                                jxix = CBFPOWX[j+jcbfind] + ix ;
	                                jyiy = CBFPOWY[j+jcbfind] + iy ;
	                                jziz = CBFPOWZ[j+jcbfind] + iz ;
                                        tij  = ti * jbasis->shells[jshell].primitives[jp].ccbf[j] ;
    		                        sx[n] += tij * xd[jxix] * yo[jyiy] * zo[jziz] ;
    		                        sy[n] += tij * xo[jxix] * yd[jyiy] * zo[jziz] ;
    		                        sz[n] += tij * xo[jxix] * yo[jyiy] * zd[jziz] ;
                                    }
                                }
                            }
                        }

                        /* . Transform the integrals. */
                        /*
                        if ( qcAtoms->QTOSPHERICAL )
                        {
	                    Integral_Block_Transform_M ( sx, iangmom, ncfunci, nfunci, jangmom, ncfuncj, nfuncj ) ;
	                    Integral_Block_Transform_M ( sy, iangmom, ncfunci, nfunci, jangmom, ncfuncj, nfuncj ) ;
	                    Integral_Block_Transform_M ( sz, iangmom, ncfunci, nfunci, jangmom, ncfuncj, nfuncj ) ;
                        }
                        */

                        /* . Put the integrals in the proper place. */
                        if ( QDIAGONAL )
                        {
                            SymmetricMatrix_Set_DBlockM ( (*dipx), istart + ibasis->shells[ishell].nstartw, ibasis->shells[ishell].nbasisw, sx ) ;
                            SymmetricMatrix_Set_DBlockM ( (*dipy), istart + ibasis->shells[ishell].nstartw, ibasis->shells[ishell].nbasisw, sy ) ;
                            SymmetricMatrix_Set_DBlockM ( (*dipz), istart + ibasis->shells[ishell].nstartw, ibasis->shells[ishell].nbasisw, sz ) ;
                        }
                        else
                        {
                            SymmetricMatrix_Set_OBlockM ( (*dipx), istart + ibasis->shells[ishell].nstartw, ibasis->shells[ishell].nbasisw, jstart + jbasis->shells[jshell].nstartw, jbasis->shells[jshell].nbasisw, sx ) ;
                            SymmetricMatrix_Set_OBlockM ( (*dipy), istart + ibasis->shells[ishell].nstartw, ibasis->shells[ishell].nbasisw, jstart + jbasis->shells[jshell].nstartw, jbasis->shells[jshell].nbasisw, sy ) ;
                            SymmetricMatrix_Set_OBlockM ( (*dipz), istart + ibasis->shells[ishell].nstartw, ibasis->shells[ishell].nbasisw, jstart + jbasis->shells[jshell].nstartw, jbasis->shells[jshell].nbasisw, sz ) ;
                        }
                    }
                }
            }
        }
    }
}

/*------------------------------------------------------------------------------
! . Calculate the electron-fit integrals.
!-----------------------------------------------------------------------------*/
# define FITINTEGRALS_BLOCKSIZE 1024
# define FITINTEGRALS_UNDERFLOW 1.0e-12

void GaussianBasis_Electron_Fit ( const QCAtomContainer *qcAtoms, const QCParameter *qcParameters, Coordinates3 *qccoordinates3, BlockStorage **fitintegrals )
{
    if ( ( qcAtoms != NULL ) && ( qcParameters != NULL ) && ( qccoordinates3 != NULL ) )
    {
        auto Boolean                QDIAGONAL, QF0, QF1, QIJ0, QIJ1, QIQMJQM, QSKIP ;
        auto Integer16               indices16[MAXCBF*MAXCBF*MAXCBF] ;
        auto Integer32               indices32[MAXCBF*MAXCBF*MAXCBF] ;
        auto Real              aa, aandb, aainv, ab, ai, aj, arri, axac, ayac, azac, bp01, b00, b10, c1x, c1y, c1z, c3x, c3y, c3z, c4x, c4y, c4z,
                                 dnuc, dxijt, dyijt, dzijt, expfac, expf, fac, fac2, f00, rho, rij2, ti, tij, tijf, u2,
                                 xc00, xcp00, xij, yc00, ycp00, yij, zc00, zcp00, zij ;
        auto Real              ar[3], ari[3], g[MAXCBF*MAXCBF*MAXCBF],
                                 xint[MAXAMP1*MAXAMP1*MAXAMP1*MAXRYS], yint[MAXAMP1*MAXAMP1*MAXAMP1*MAXRYS], zint[MAXAMP1*MAXAMP1*MAXAMP1*MAXRYS] ;
        auto Real             *rc, *rf, *ri, *rj ;
        auto Integer                 dim1, dim2, dim3,
                                 f, fammax,          fc, fcbfind, fijx, fijy, fijz, fqm, fp, fshell, fstart,
                                 i, iammax, iammaxt, ic, icbfind, ii, ij, iqm, ip, ishell, istart, ix, iy, iz,
                                 j, jammax, jammaxt, jc, jcbfind, jix, jiy, jiz, jj, jp, jqm, jshell, jstart,
                                 jupper, m, n, ncfuncf, ncfunci, ncfuncj, nroots ;
        auto GaussianBasis *ibasis, *jbasis, *fbasis ;
        auto RysQuadrature  roots ;
# ifdef PRINTINTEGRALS
auto Integer ntotal = 0 ;
printf ( "\nElectron-Fit Integrals:\n\n" ) ;
# endif
        /* . Initialization. */
        (*fitintegrals) = BlockStorage_Allocate ( ) ;
        (*fitintegrals)->blocksize  = FITINTEGRALS_BLOCKSIZE ;
        (*fitintegrals)->nindices16 = 1 ;
        (*fitintegrals)->nindices32 = 1 ;
        (*fitintegrals)->QUNDERFLOW = True ;
        (*fitintegrals)->underflow  = FITINTEGRALS_UNDERFLOW ;

        /*----------------------------------------------------------------------
        ! . Triple loop over centers.
        !---------------------------------------------------------------------*/
        /* . Outer loop over centers. */
        for ( iqm = 0 ; iqm < qcAtoms->natoms ; iqm++ )
        {

            /* . Get data for the center. */
            ic     = qcAtoms->data[iqm].center ;
            ibasis = qcParameters->centers[ic].orbitalbasis ;
            istart = qcAtoms->data[iqm].ostartw ;
            ri     = Coordinates3_RowPointer ( qccoordinates3, iqm ) ;

            /* . Inner loop over centers. */
            for ( jqm = 0 ; jqm <= iqm ; jqm++ )
            {

                /* . Get data for the center. */
                jc     = qcAtoms->data[jqm].center ;
                jbasis = qcParameters->centers[jc].orbitalbasis ;
                jstart = qcAtoms->data[jqm].ostartw ;
                rj     = Coordinates3_RowPointer ( qccoordinates3, jqm ) ;

                /* . Set the diagonal atom flag. */
                QIQMJQM = ( iqm == jqm ) ;

                /* . Calculate some distance factors. */
	        xij  = ri[0] - rj[0] ;
	        yij  = ri[1] - rj[1] ;
	        zij  = ri[2] - rj[2] ;
	        rij2 = xij * xij + yij * yij + zij * zij ;

                /* . Loop over fitting atoms. */
                for ( fqm = 0 ; fqm < qcAtoms->natoms ; fqm++ )
                {

                    /* . Get information about the atom. */
                    fc     = qcAtoms->data[fqm].center ;
                    fbasis = qcParameters->centers[fc].densitybasis ;
                    fstart = qcAtoms->data[fqm].dstartw ;
                    rf     = Coordinates3_RowPointer ( qccoordinates3, fqm ) ;

                    /* . Get some other distances. */
/*
                    for ( i = 0, rfi2 = rfj2 = 0.0e+00 ; i < 3 ; i++ )
                    {
                        faci  = rf[i] - ri[i] ;
                        facj  = rf[i] - rj[i] ;
                        rfi2 += faci * faci ;
                        rfj2 += facj * facj ;
                    }
*/

                    /*----------------------------------------------------------
                    ! . Triple loop over shells.
                    !---------------------------------------------------------*/
                    /* . Outer loop over shells. */
                    for ( ishell = 0 ; ishell < ibasis->nshells ; ishell++ )
                    {

                        /* . Get information about the shell. */
                        iammax  = ibasis->shells[ishell].type->angularmomentum_high ;
                        icbfind = ibasis->shells[ishell].type->cbfindex ;
                        ncfunci = ibasis->shells[ishell].type->ncbf     ;

                        /* . Set the upper limit for the JSHELL loops. */
                        if ( QIQMJQM ) jupper = ishell + 1 ;
                        else           jupper = jbasis->nshells ;

                        /* . Inner loop over shells. */
                        for ( jshell = 0 ; jshell < jupper ; jshell++ )
                        {

                            /* . Get information about the shell. */
                            jammax  = jbasis->shells[jshell].type->angularmomentum_high ;
                            jcbfind = jbasis->shells[jshell].type->cbfindex ;
                            ncfuncj = jbasis->shells[jshell].type->ncbf     ;

                            /* . Set the diagonal block flag. */
                            QDIAGONAL = QIQMJQM && ( ishell == jshell ) ;

                            /* . Set some flags. */
                            QIJ0 = ( iammax + jammax == 0 ) ;
                            QIJ1 = ( iammax + jammax <= 1 ) ;

                            /* . Select the expansion center for the recurrence relations. */
                            if ( iammax >= jammax )
                            {
                               iammaxt = iammax ;
                               jammaxt = jammax ;
                               dxijt   = xij ;
                               dyijt   = yij ;
                               dzijt   = zij ;
                               rc      = ri ;
                            }
                            else
                            {
                               iammaxt = jammax ;
                               jammaxt = iammax ;
                               dxijt   = - xij ;
                               dyijt   = - yij ;
                               dzijt   = - zij ;
                               rc      = rj ;
                            }


                            /* . Loop over fitting shells. */
                            for ( fshell = 0 ; fshell < fbasis->nshells ; fshell++ )
                            {

                                /* . Get information about the shell. */
                                fammax  = fbasis->shells[fshell].type->angularmomentum_high ;
                                fcbfind = fbasis->shells[fshell].type->cbfindex ;
                                ncfuncf = fbasis->shells[fshell].type->ncbf     ;

                                /* . Set some flags. */
                                QF0 = ( fammax == 0 ) ;
                                QF1 = ( fammax <= 1 ) ;

                                /* . Get the number of roots. */
                                nroots = ( fammax + iammax + jammax ) / 2 + 1 ;

                                /* . Initialize the integral blocks. */
                                for ( i = 0 ; i < ( ncfuncf * ncfunci * ncfuncj ) ; i++ ) g[i] = 0.0e+00 ;

                                /* . Set the array dimensions. */
                                dim1 = fammax + 1 ;
                                if ( iammax >= jammax ) dim2 = dim1 * ( jammax + 1 ) ;
                                else                    dim2 = dim1 * ( iammax + 1 ) ;
                                dim3 = dim1 * ( iammax + 1 ) * ( jammax + 1 ) ;

                                /*----------------------------------------------
                                ! . Triple loop over primitives.
                                !---------------------------------------------*/
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
                                        expfac = exp ( - fac ) * PI252 * aainv ;
	                                for ( i = 0 ; i < 3 ; i++ ) ar[i] = ( ari[i] + aj * rj[i] ) * aainv ;

                                        /* . Loop over fitting primitives. */
                                        for ( fp = 0 ; fp < fbasis->shells[fshell].nprimitives ; fp++ )
                                        {
                                            /* . Get some information for the primitive. */
	                                    expf = fbasis->shells[fshell].primitives[fp].exponent ;

                                            /* . Calculate some factors. */
                                            ab    = aa * expf ;
                                            aandb = aa + expf ;
                                            rho   = ab / aandb ;
                                            dnuc  = expfac / ( expf * sqrt ( aandb ) ) ;

                                            /* . Calculate the rys polynomial roots. */
                                            c1x = ( ar[0] - rf[0] ) ;
                                            c1y = ( ar[1] - rf[1] ) ;
                                            c1z = ( ar[2] - rf[2] ) ;
                                            RysQuadrature_Roots ( &roots, nroots, rho * ( c1x * c1x + c1y * c1y + c1z * c1z ) ) ;

                                            /* . Calculate some displacements. */
                                            axac = aa * ( ar[0] - rc[0] ) ;
                                            ayac = aa * ( ar[1] - rc[1] ) ;
                                            azac = aa * ( ar[2] - rc[2] ) ;
                                            c1x *= aa ;
                                            c1y *= aa ;
                                            c1z *= aa ;
                                            c3x  = expf * ( rf[0] - rc[0] ) + axac ;
                                            c3y  = expf * ( rf[1] - rc[1] ) + ayac ;
                                            c3z  = expf * ( rf[2] - rc[2] ) + azac ;
                                            c4x  = expf * axac ;
                                            c4y  = expf * ayac ;
                                            c4z  = expf * azac ;

                                            /* . Loop over the roots and construct the subsidiary integrals. */
                                            for ( m = 0 ; m < nroots ; m++ )
                                            {
                                                u2    = roots.roots[m] * rho ;
                                                f00   = roots.weights[m] ;
                                                fac   = 1.0e+00 / ( ab + u2 * aandb ) ;
                                                fac2  = 0.5e+00 * fac ;
                                                bp01  = ( aa   + u2 ) * fac2 ;
                                                b00   =          u2   * fac2 ;
                                                b10   = ( expf + u2 ) * fac2 ;
                                                xcp00 = u2 * c1x * fac ;
                                                ycp00 = u2 * c1y * fac ;
                                                zcp00 = u2 * c1z * fac ;
                                                xc00  = ( u2 * c3x + c4x ) * fac ;
                                                yc00  = ( u2 * c3y + c4y ) * fac ;
                                                zc00  = ( u2 * c3z + c4z ) * fac ;
                                                Subsidiary_Integral_Nuclear3C ( iammaxt, jammaxt, fammax, QIJ0, QIJ1, QF0, QF1, b00, b10, bp01, dxijt, dyijt, dzijt, f00,
                                                                                xc00, xcp00, yc00, ycp00, zc00, zcp00, dim1, dim2, &xint[m*dim3], &yint[m*dim3], &zint[m*dim3] ) ;
                                            }

                                            /* . Assemble the integrals. */
                                            if ( iammax >= jammax )
                                            {
                                                for ( i = 0, n = 0 ; i < ncfunci ; i++ )
                                                {
   	                                            ix = CBFPOWX[i+icbfind] * dim2 ;
	                                            iy = CBFPOWY[i+icbfind] * dim2 ;
	                                            iz = CBFPOWZ[i+icbfind] * dim2 ;
                                                    ti = dnuc * ibasis->shells[ishell].primitives[ip].ccbf[i] ;
                                                    for ( j = 0 ; j < ncfuncj ; j++ )
                                                    {
	                                                jix = CBFPOWX[j+jcbfind] * dim1 + ix ;
	                                                jiy = CBFPOWY[j+jcbfind] * dim1 + iy ;
	                                                jiz = CBFPOWZ[j+jcbfind] * dim1 + iz ;
                                                        tij = ti * jbasis->shells[jshell].primitives[jp].ccbf[j] ;
                                                        for ( f = 0 ; f < ncfuncf ; f++, n++ )
                                                        {
                                                            fijx = CBFPOWX[f+fcbfind] + jix ;
                                                            fijy = CBFPOWY[f+fcbfind] + jiy ;
                                                            fijz = CBFPOWZ[f+fcbfind] + jiz ;
                                                            for ( m = 0, fac = 0.0e+00 ; m < nroots ; m++ ) fac += xint[fijx+m*dim3] * yint[fijy+m*dim3] * zint[fijz+m*dim3] ;
                                                            tijf = tij * fbasis->shells[fshell].primitives[fp].ccbf[f] ;
                                                            g[n] += tijf * fac ;
                                                        }
                                                    }
                                                }
                                            }
                                            else
                                            {
                                                for ( j = 0 ; j < ncfuncj ; j++ )
                                                {
   	                                            ix = CBFPOWX[j+jcbfind] * dim2 ;
	                                            iy = CBFPOWY[j+jcbfind] * dim2 ;
	                                            iz = CBFPOWZ[j+jcbfind] * dim2 ;
                                                    ti = dnuc * jbasis->shells[jshell].primitives[jp].ccbf[j] ;
                                                    for ( i = 0 ; i < ncfunci ; i++ )
                                                    {
                                                        n   = j * ncfuncf + i * ( ncfuncf * ncfuncj ) ;
	                                                jix = CBFPOWX[i+icbfind] * dim1 + ix ;
	                                                jiy = CBFPOWY[i+icbfind] * dim1 + iy ;
	                                                jiz = CBFPOWZ[i+icbfind] * dim1 + iz ;
                                                        tij = ti * ibasis->shells[ishell].primitives[ip].ccbf[i] ;
                                                        for ( f = 0 ; f < ncfuncf ; f++, n++ )
                                                        {
                                                            fijx = CBFPOWX[f+fcbfind] + jix ;
                                                            fijy = CBFPOWY[f+fcbfind] + jiy ;
                                                            fijz = CBFPOWZ[f+fcbfind] + jiz ;
                                                            for ( m = 0, fac = 0.0e+00 ; m < nroots ; m++ ) fac += xint[fijx+m*dim3] * yint[fijy+m*dim3] * zint[fijz+m*dim3] ;
                                                            tijf = tij * fbasis->shells[fshell].primitives[fp].ccbf[f] ;
                                                            g[n] += tijf * fac ;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }

                                /* . Transform the integrals. */
                                /* if ( qcAtoms->QTOSPHERICAL ) Integral_Block_O3Transform ( g, fangmom, ncfuncf, nsfuncf, iangmom, ncfunci, nfunci, jangmom, ncfuncj, nfuncj ) ; */

                                /* . Generate the final integral and index arrays. */
                                for ( i = 0, m = 0, n = 0 ; i < ibasis->shells[ishell].nbasisw ; i++ )
                                {
                                    ii = istart + ibasis->shells[ishell].nstartw + i ;
                                    for ( j = 0 ; j < jbasis->shells[jshell].nbasisw ; j++ )
                                    {
                                        QSKIP = QDIAGONAL && ( j > i ) ;
                                        jj    = jstart + jbasis->shells[jshell].nstartw + j ;
                                        ij    = ( ii * ( ii + 1 ) ) / 2 + jj ;
                                        for ( f = 0 ; f < fbasis->shells[fshell].nbasisw ; f++, n++ )
                                        {
                                            if ( ! QSKIP )
                                            {
                                                indices16[m] = fstart + fbasis->shells[fshell].nstartw + f ;
                                                indices32[m] = ij ;
                                                g[m] = g[n] ;
# ifdef PRINTINTEGRALS
ntotal += 1 ;
printf ( "%10d %25.15f %10d %10d\n", ntotal, g[m], indices16[m], indices32[m] ) ;
# endif
                                                m++ ;
                                            }
                                        }
                                    }
                                }

                                /* . Save the integrals. */
                                BlockStorage_Data_Add ( (*fitintegrals), m, g, indices16, indices32 ) ;
                            }
                        }
                    }
                }
            }
        }
    }
}

# undef FITINTEGRALS_BLOCKSIZE
# undef FITINTEGRALS_UNDERFLOW

/*------------------------------------------------------------------------------
! . Calculate the electron-nuclear integrals.
!-----------------------------------------------------------------------------*/
void GaussianBasis_Electron_Nuclear ( const QCAtomContainer *qcAtoms, const QCParameter *qcParameters, Coordinates3 *qccoordinates3, SymmetricMatrix *oneelectronmatrix )
{
    if ( ( qcAtoms != NULL ) && ( qcParameters != NULL ) && ( qccoordinates3 != NULL ) && ( oneelectronmatrix != NULL ) )
    {
        auto Boolean                QDIAGONAL, QIJ0, QIJ1, QIQMJQM ;
        auto Real              aa, aandb, aainv, ab, ai, aj, arri, axac, ayac, azac, bp01, b00, b10, c1x, c1y, c1z, c3x, c3y, c3z, c4x, c4y, c4z,
                                 dnuc, dxijt, dyijt, dzijt, expfac, expn, fac, facn, fac2, f00, qn, rho, rij2, ti, tij, u2,
                                 xc00, xcp00, xij, yc00, ycp00, yij, zc00, zcp00, zij ;
        auto Real              ar[3], ari[3], g[MAXCBF*MAXCBF], xint[MAXAMP1*MAXAMP1*MAXRYS], yint[MAXAMP1*MAXAMP1*MAXRYS], zint[MAXAMP1*MAXAMP1*MAXRYS] ;
        auto Real             *rc, *ri, *rj, *rn ;
        auto Integer                 i, iammax, iammaxt, ic, icbfind, iqm, ip, ishell, istart, ix, iy, iz,
                                 j, jammax, jammaxt, jc, jcbfind, jdim, jdimm, jp, jqm, jshell, jstart,
                                 jupper, jxix, jyiy, jziz, kqm, m, n, ncfunci, ncfuncj, nroots ;
        auto GaussianBasis *ibasis, *jbasis ;
        auto RysQuadrature  roots ;

        /* . Outer loop over centers. */
        for ( iqm = 0 ; iqm < qcAtoms->natoms ; iqm++ )
        {

            /* . Get data for the center. */
            ic     = qcAtoms->data[iqm].center ;
            ibasis = qcParameters->centers[ic].orbitalbasis ;
            istart = qcAtoms->data[iqm].ostartw ;
            ri     = Coordinates3_RowPointer ( qccoordinates3, iqm ) ;

            /* . Inner loop over centers. */
            for ( jqm = 0 ; jqm <= iqm ; jqm++ )
            {

                /* . Get data for the center. */
                jc     = qcAtoms->data[jqm].center ;
                jbasis = qcParameters->centers[jc].orbitalbasis ;
                jstart = qcAtoms->data[jqm].ostartw ;
                rj     = Coordinates3_RowPointer ( qccoordinates3, jqm ) ;

                /* . Set the diagonal atom flag. */
                QIQMJQM = ( iqm == jqm ) ;

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

                    /* . Set the upper limit for the JSHELL loops. */
                    if ( QIQMJQM ) jupper = ishell + 1 ;
                    else           jupper = jbasis->nshells ;

                    /* . Inner loop over shells. */
                    for ( jshell = 0 ; jshell < jupper ; jshell++ )
                    {

                        /* . Get information about the shell. */
                        jammax  = jbasis->shells[jshell].type->angularmomentum_high ;
                        jdimm   = ( iammax + 1 ) * ( jammax + 1 ) ;
                        jcbfind = jbasis->shells[jshell].type->cbfindex ;
                        ncfuncj = jbasis->shells[jshell].type->ncbf     ;

                        /* . Set the diagonal block flag. */
                        QDIAGONAL = QIQMJQM && ( ishell == jshell ) ;

                        /* . Get the number of roots. */
                        nroots = ( iammax + jammax ) / 2 + 1 ;

                        /* . Set some flags. */
                        QIJ0 = ( iammax + jammax == 0 ) ;
                        QIJ1 = ( iammax + jammax <= 1 ) ;

                        /* . Initialize the integral blocks. */
                        for ( i = 0 ; i < ( ncfunci * ncfuncj ) ; i++ ) g[i] = 0.0e+00 ;

                        /* . Select the expansion center for the recurrence relations. */
                        if ( iammax >= jammax )
                        {
                           iammaxt = iammax ;
                           jammaxt = jammax ;
                           jdim    = jammax + 1 ;
                           dxijt   = xij ;
                           dyijt   = yij ;
                           dzijt   = zij ;
                           rc      = ri ;
                        }
                        else
                        {
                           iammaxt = jammax ;
                           jammaxt = iammax ;
                           jdim    = iammax + 1 ;
                           dxijt   = - xij ;
                           dyijt   = - yij ;
                           dzijt   = - zij ;
                           rc      = rj ;
                        }

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
                                expfac = exp ( - fac ) * PI252 * aainv ;
	                        for ( i = 0 ; i < 3 ; i++ ) ar[i] = ( ari[i] + aj * rj[i] ) * aainv ;

                                /* . Loop over the nuclear densities. */
                                for ( kqm = 0 ; kqm < qcAtoms->natoms ; kqm++ )
                                {

                                    /* . Get information about the atom. */
                                    expn = qcAtoms->data[kqm].widthe ;
                                    facn = qcAtoms->data[kqm].widthn ;
                                    qn   = - ( Real ) qcAtoms->data[kqm].atomicNumber ;
                                    rn   = Coordinates3_RowPointer ( qccoordinates3, kqm ) ;

                                    /* . Calculate some factors. */
                                    ab     = aa * expn ;
                                    aandb  = aa + expn ;
                                    rho    = ab / aandb ;
                                    dnuc   = expfac * ( facn * qn ) / ( expn * sqrt ( aandb ) ) ;

                                    /* . Calculate the rys polynomial roots. */
                                    c1x = ( ar[0] - rn[0] ) ;
                                    c1y = ( ar[1] - rn[1] ) ;
                                    c1z = ( ar[2] - rn[2] ) ;
                                    RysQuadrature_Roots ( &roots, nroots, rho * ( c1x * c1x + c1y * c1y + c1z * c1z ) ) ;

                                    /* . Calculate some displacements. */
                                    axac = aa * ( ar[0] - rc[0] ) ;
                                    ayac = aa * ( ar[1] - rc[1] ) ;
                                    azac = aa * ( ar[2] - rc[2] ) ;
                                    c1x *= aa ;
                                    c1y *= aa ;
                                    c1z *= aa ;
                                    c3x  = expn * ( rn[0] - rc[0] ) + axac ;
                                    c3y  = expn * ( rn[1] - rc[1] ) + ayac ;
                                    c3z  = expn * ( rn[2] - rc[2] ) + azac ;
                                    c4x  = expn * axac ;
                                    c4y  = expn * ayac ;
                                    c4z  = expn * azac ;

                                    /* . Loop over the roots and construct the subsidiary integrals. */
                                    for ( m = 0 ; m < nroots ; m++ )
                                    {
                                        u2    = roots.roots[m] * rho ;
                                        f00   = roots.weights[m] ;
                                        fac   = 1.0e+00 / ( ab + u2 * aandb ) ;
                                        fac2  = 0.5e+00 * fac ;
                                        bp01  = ( aa   + u2 ) * fac2 ;
                                        b00   =          u2   * fac2 ;
                                        b10   = ( expn + u2 ) * fac2 ;
                                        xcp00 = u2 * c1x * fac ;
                                        ycp00 = u2 * c1y * fac ;
                                        zcp00 = u2 * c1z * fac ;
                                        xc00  = ( u2 * c3x + c4x ) * fac ;
                                        yc00  = ( u2 * c3y + c4y ) * fac ;
                                        zc00  = ( u2 * c3z + c4z ) * fac ;
                                        Subsidiary_Integral_Nuclear3C ( iammaxt, jammaxt, 0, QIJ0, QIJ1, True, True, b00, b10, bp01, dxijt, dyijt,dzijt, f00,
                                                                        xc00, xcp00, yc00, ycp00, zc00, zcp00, 1, jdim, &xint[m*jdimm], &yint[m*jdimm], &zint[m*jdimm] ) ;
                                    }

                                    /* . Assemble the integrals. */
                                    if ( iammax >= jammax )
                                    {
                                        for ( i = 0, n = 0 ; i < ncfunci ; i++ )
                                        {
   	                                    ix = CBFPOWX[i+icbfind] * jdim ;
	                                    iy = CBFPOWY[i+icbfind] * jdim ;
	                                    iz = CBFPOWZ[i+icbfind] * jdim ;
                                            ti = dnuc * ibasis->shells[ishell].primitives[ip].ccbf[i] ;
                                            for ( j = 0 ; j < ncfuncj ; j++, n++ )
                                            {
	                                        jxix = CBFPOWX[j+jcbfind] + ix ;
	                                        jyiy = CBFPOWY[j+jcbfind] + iy ;
	                                        jziz = CBFPOWZ[j+jcbfind] + iz ;
                                                for ( m = 0, fac = 0.0e+00 ; m < nroots ; m++ ) fac += xint[jxix+m*jdimm] * yint[jyiy+m*jdimm] * zint[jziz+m*jdimm] ;
                                                tij = ti * jbasis->shells[jshell].primitives[jp].ccbf[j] ;
                                                g[n] += tij * fac ;
                                            }
                                        }
                                    }
                                    else
                                    {
                                        for ( j = 0 ; j < ncfuncj ; j++ )
                                        {
   	                                    ix = CBFPOWX[j+jcbfind] * jdim ;
	                                    iy = CBFPOWY[j+jcbfind] * jdim ;
	                                    iz = CBFPOWZ[j+jcbfind] * jdim ;
                                            ti = dnuc * jbasis->shells[jshell].primitives[jp].ccbf[j] ;
                                            for ( i = 0, n = j ; i < ncfunci ; i++, n+= ncfuncj )
                                            {
	                                        jxix = CBFPOWX[i+icbfind] + ix ;
	                                        jyiy = CBFPOWY[i+icbfind] + iy ;
	                                        jziz = CBFPOWZ[i+icbfind] + iz ;
                                                for ( m = 0, fac = 0.0e+00 ; m < nroots ; m++ ) fac += xint[jxix+m*jdimm] * yint[jyiy+m*jdimm] * zint[jziz+m*jdimm] ;
                                                tij = ti * ibasis->shells[ishell].primitives[ip].ccbf[i] ;
                                                g[n] += tij * fac ;
                                            }
                                        }
                                    }
                                }
                            }
                        }

                        /* . Transform the integrals. */
	                /* if ( qcAtoms->QTOSPHERICAL ) Integral_Block_Transform_M ( g, iangmom, ncfunci, nfunci, jangmom, ncfuncj, nfuncj ) ; */

                        /* . Put the integrals into the proper place. */
                        if ( QDIAGONAL )
                        {
                            SymmetricMatrix_Increment_DBlockM ( oneelectronmatrix, istart + ibasis->shells[ishell].nstartw, ibasis->shells[ishell].nbasisw, g ) ;
                        }
                        else
                        {
                            SymmetricMatrix_Increment_OBlockM ( oneelectronmatrix, istart + ibasis->shells[ishell].nstartw, ibasis->shells[ishell].nbasisw, jstart + jbasis->shells[jshell].nstartw, jbasis->shells[jshell].nbasisw, g ) ;
                        }
                    }
                }
            }
        }
    }
}

/*------------------------------------------------------------------------------
! . Calculate the fit-fit integrals.
!-----------------------------------------------------------------------------*/
# define DOCOTRANSFORMATION /* . This option is for checking whether c->o transformation is really necessary when calculating the inverse. */
# define INVERSEFITTOL 1.0e-5
void GaussianBasis_Fit_Fit ( const QCAtomContainer *qcAtoms, const QCParameter *qcParameters, Coordinates3 *qccoordinates3, const Real1DArray *fselfoverlap, SymmetricMatrix **inversefitmatrix )
{
    if ( ( qcAtoms != NULL ) && ( qcParameters != NULL ) && ( qccoordinates3 != NULL ) && ( fselfoverlap != NULL ) )
    {
        auto Boolean                  QDIAGONAL, QIQMJQM ;
        auto Real                aandb, ab, ai, aj, bp01, b00, b10, c1x, c1y, c1z, c3x, c3y, c3z, dfi, dfij, fac, fac2, f00,
                                   rho, rij2, ti, tij, u2, xc00, xcp00, xij, yc00, ycp00, yij, zc00, zcp00, zij ;
        auto Real                g[MAXCBF*MAXCBF], xint[MAXAMP1*MAXAMP1*MAXRYS], yint[MAXAMP1*MAXAMP1*MAXRYS], zint[MAXAMP1*MAXAMP1*MAXRYS] ;
        auto Real               *ri, *rj ;
        auto Integer                   i, iammax, ic, icbfind, ij, iqm, ip, ishell, istart, ivec, ix, iy, iz,
                                   j, jammax, jc, jcbfind, jdim, jdimm, jp, jqm, jshell, jstart, jupper, jxix, jyiy, jziz, m, n, ncfunci, ncfuncj, ndim, nroots ;
        auto GaussianBasis   *ibasis, *jbasis ;
        auto Real2DArray          *c2o, *eigenvectors ;
        auto RysQuadrature    roots ;
        auto SymmetricMatrix *eri, *temporary ;
        auto Real1DArray          *eigenvalues ;

        /* . Initialization. */
        ndim = qcAtoms->nfbasisw + 1 ;
        eri  = SymmetricMatrix_Allocate ( ndim ) ;
        SymmetricMatrix_Set_Zero ( eri ) ;

        /* . Outer loop over centers. */
        for ( iqm = 0 ; iqm < qcAtoms->natoms ; iqm++ )
        {

            /* . Get data for the center. */
            ic     = qcAtoms->data[iqm].center ;
            ibasis = qcParameters->centers[ic].densitybasis ;
            istart = qcAtoms->data[iqm].dstartw ;
            ri     = Coordinates3_RowPointer ( qccoordinates3, iqm ) ;

            /* . Inner loop over centers. */
            for ( jqm = 0 ; jqm <= iqm ; jqm++ )
            {

                /* . Get data for the center. */
                jc     = qcAtoms->data[jqm].center ;
                jbasis = qcParameters->centers[jc].densitybasis ;
                jstart = qcAtoms->data[jqm].dstartw ;
                rj     = Coordinates3_RowPointer ( qccoordinates3, jqm ) ;

                /* . Set the diagonal atom flag. */
                QIQMJQM = ( iqm == jqm ) ;

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

                    /* . Set the upper limit for the JSHELL loops. */
                    if ( QIQMJQM ) jupper = ishell + 1 ;
                    else           jupper = jbasis->nshells ;

                    /* . Inner loop over shells. */
                    for ( jshell = 0 ; jshell < jupper ; jshell++ )
                    {

                        /* . Get information about the shell. */
                        jammax  = jbasis->shells[jshell].type->angularmomentum_high ;
                        jdim    = jammax + 1 ;
                        jdimm   = ( iammax + 1 ) * ( jammax + 1 ) ;
                        jcbfind = jbasis->shells[jshell].type->cbfindex ;
                        ncfuncj = jbasis->shells[jshell].type->ncbf     ;

                        /* . Set the diagonal block flag. */
                        QDIAGONAL = QIQMJQM && ( ishell == jshell ) ;

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
                                    Subsidiary_Integral_Nuclear2C ( iammax, jammax, b00, b10, bp01, f00, xc00, xcp00, yc00, ycp00, zc00, zcp00,
                                                                    jdim, &xint[m*jdimm], &yint[m*jdimm], &zint[m*jdimm] ) ;
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

                        /* . Put the integrals into the proper place. */
                        if ( QDIAGONAL )
                        {
                            SymmetricMatrix_Set_DBlockM ( eri, istart + ibasis->shells[ishell].nstartw, ibasis->shells[ishell].nbasisw, g ) ;
                        }
                        else
                        {
                            SymmetricMatrix_Set_OBlockM ( eri, istart + ibasis->shells[ishell].nstartw, ibasis->shells[ishell].nbasisw, jstart + jbasis->shells[jshell].nstartw, jbasis->shells[jshell].nbasisw, g ) ;
                        }
                    }
                }
            }
        }
# ifdef PRINTDIAGONALELEMENTS
{
auto Integer i ;
printf ( "\nDIAGONAL FIT/FIT ELEMENTS:\n" ) ;
for ( i = 0 ; i < ndim - 1 ; i++ )
{
    printf ( "FF> %10d %25.15f %25.15f\n", i, SymmetricMatrix_Get_Component ( eri, i, i ), SymmetricMatrix_Get_Component ( eri, i, i ) - 1.0e+00 ) ;
}
}
# endif
        /*----------------------------------------------------------------------
        ! . Get the inverse fit matrix.
        !---------------------------------------------------------------------*/
        /* . Complete the specification of the eri matrix for the charge constraint. */
        for ( i = 0 ; i < qcAtoms->nfbasisw ; i++ ) eri->data[(qcAtoms->nfbasisw*(qcAtoms->nfbasisw+1))/2+i] = -fselfoverlap->data[i] ;
        eri->data[eri->size-1] = 0.0e+00 ;
# ifdef PRINTINTEGRALS
{
auto Integer ij ;
printf ( "\nFit-Fit Matrix:\n\n" ) ;
for ( i = 0, ij = 0 ; i < ndim ; i++ )
{
for ( j = 0 ; j <= i ; j++, ij++ ) printf ( "%10d %10d %25.15f\n", i, j, eri->data[ij] ) ;
}
}
# endif


# ifdef DOCOTRANSFORMATION
        /* . Get the transformation c->o. */
        {
            auto Integer i, ic, iqm, istart, istartw, j ;
            auto GaussianBasis *ibasis ;
            auto Real2DArray *m ;

            /* . Reallocate space. */
            c2o = Real2DArray_Allocate ( ndim, qcAtoms->nfbasis + 1, NULL ) ;
            Real2DArray_Set ( c2o, 0.0e+00 ) ;

            /* . Loop over the orbital bases for each atom. */
            for ( iqm = 0 ; iqm < qcAtoms->natoms ; iqm++ )
            {
                ic      = qcAtoms->data[iqm].center ;
                ibasis  = qcParameters->centers[ic].densitybasis ;
                istart  = qcAtoms->data[iqm].fstart  ;
                istartw = qcAtoms->data[iqm].fstartw ;
                m = ibasis->c2o ;
                for ( i = 0 ; i < m->length0 ; i++ )
                {
                    for ( j = 0 ; j < m->length1 ; j++ ) Real2DArray_Item ( c2o, i + istartw, j + istart ) =  Real2DArray_Item ( m, i, j ) ;
                }
            }
        }

        /* . Get the new dimension for the problem. */
        ndim = qcAtoms->nfbasis + 1 ;

        /* . Allocate space for the reduced matrix. */
        temporary = SymmetricMatrix_Allocate ( ndim ) ;
        SymmetricMatrix_Set ( temporary, 0.0e+00 ) ;

        /* . Forward transform. */
        SymmetricMatrix_Transform ( eri, c2o, False, temporary ) ;
# else
        /* . Alias. */
        temporary = eri ;
# endif

        /* . Get the inverse of the matrix in temporary. */
        /* . Allocate more space. */
        eigenvalues  = Real1DArray_Allocate ( ndim      , NULL ) ;
        eigenvectors = Real2DArray_Allocate ( ndim, ndim, NULL ) ;

        /* . Diagonalize the matrix. */
        SymmetricMatrix_Diagonalize ( temporary, eigenvalues, eigenvectors, NULL ) ;

        /* . Construct the inverse matrix. */
        SymmetricMatrix_Set_Zero ( temporary ) ;

        /* . Double loop over eigenvectors. */
        for ( ivec = 0 ; ivec < ndim ; ivec++ )
        {
            if ( fabs ( eigenvalues->data[ivec] ) > INVERSEFITTOL )
            {
                fac = 1.0e+00 / eigenvalues->data[ivec] ;
                for ( i = 0, ij = 0 ; i < ndim ; i++ )
                {
                    for ( j = 0 ; j <= i ; ij++, j++ ) temporary->data[ij] += fac * eigenvectors->data[ivec+i*ndim] * eigenvectors->data[ivec+j*ndim] ;
                }
            }
# ifdef PRINTSMALLEIGENVALUES
            else
            {
                printf ( "SMALL FIT-FIT EIGENVALUE %10d %25.15f\n", ivec,  eigenvalues->data[ivec] ) ;
            }
# endif
        }

        /* . Deallocate space. */
        Real2DArray_Deallocate ( &eigenvectors ) ;
        Real1DArray_Deallocate ( &eigenvalues  ) ;

# ifdef DOCOTRANSFORMATION
        /* . Alias. */
        (*inversefitmatrix) = eri ;
        SymmetricMatrix_Set_Zero ( (*inversefitmatrix) ) ;

        /* . Back transform. */
        SymmetricMatrix_Transform ( temporary, c2o, True, (*inversefitmatrix) ) ;

        /* . Deallocate space. */
        Real2DArray_Deallocate          ( &c2o       ) ;
        SymmetricMatrix_Deallocate ( &temporary ) ;
# else
        /* . Alias. */
        (*inversefitmatrix) = temporary ;
# endif
    }
}
# undef INVERSEFITTOL
# undef DOCOTRANSFORMATION

/*------------------------------------------------------------------------------
! . Form the two-electron parts of the Fock matrices.
! . These parts are the same for alpha and beta spins.
! . May need to change the size of the block indices here.
!-----------------------------------------------------------------------------*/
double GaussianBasis_Fock ( const Integer nfbasis, QCOnePDM *density, BlockStorage *fitintegrals, Real1DArray *fpotential, SymmetricMatrix *inversefitmatrix )
{
    Real etei = 0.0e+00 ;
    if ( ( density != NULL ) )
    {
        auto Block *block ;
        auto Integer         i, ndim ;
        auto Real1DArray *b ;

        /* . Initialization. */
        SymmetricMatrix_Set_Zero ( density->fock ) ;

        /* . Scale the diagonal elements of the density by 1/2. */
        SymmetricMatrix_Scale_Diagonal ( density->density, 0, density->density->dimension, 0.5e+00 ) ;

        /* . Allocate space. */
        ndim = nfbasis + 1 ;
        b = Real1DArray_Allocate ( ndim, NULL ) ;
        Real1DArray_Set ( b, 0.0e+00 ) ;
        b->data[nfbasis] = - density->totalCharge ;

        /* . Loop over the integral blocks. */
        List_Iterate_Initialize ( fitintegrals->blocks ) ;
        while ( ( block = BlockStorage_Iterate ( fitintegrals ) ) != NULL )
        {
            for ( i = 0 ; i < block->ndata ; i++ ) b->data[block->indices16[i]] += density->density->data[block->indices32[i]] * block->data[i] ;
        }

        /* . Scale b by 2. */
        for ( i = 0 ; i < nfbasis ; i++ ) b->data[i] *= 2.0e+00 ;

        /* . Find the fitted potential by multiplying by the inverse matrix. */
        SymmetricMatrix_VectorMultiply ( inversefitmatrix, b, fpotential, NULL ) ;

        /* . Loop over the integral blocks. */
        List_Iterate_Initialize ( fitintegrals->blocks ) ;
        while ( ( block = BlockStorage_Iterate ( fitintegrals ) ) != NULL )
        {
            for ( i = 0 ; i < block->ndata ; i++ ) density->fock->data[block->indices32[i]] += fpotential->data[block->indices16[i]] * block->data[i] ;
        }

        /* . Get the electronic energy. */
        etei  = Real1DArray_Dot ( b, fpotential, NULL ) ;
        etei *= 0.5e+00 ;

        /* . Deallocate space. */
        Real1DArray_Deallocate ( &b ) ;

        /* . Scale the diagonal elements of the density by 2. */
        SymmetricMatrix_Scale_Diagonal ( density->density, 0, density->density->dimension, 2.0e+00 ) ;
    }
    return etei ;
}

/*------------------------------------------------------------------------------
! . Calculate the kinetic energy and overlap integrals.
!-----------------------------------------------------------------------------*/
void GaussianBasis_Kinetic_2Overlap ( const QCAtomContainer *qcAtoms, const QCParameter *qcParameters, Coordinates3 *qccoordinates3, SymmetricMatrix **kinetic, SymmetricMatrix **overlap )
{
    if ( ( qcAtoms != NULL ) && ( qcParameters != NULL ) && ( qccoordinates3 != NULL ) )
    {
        auto Boolean                QDIAGONAL, QIQMJQM ;
        auto Real              aa, aainv, ai, aj, arri, expfac, fac, rij2, ti, tij, xij, yij, zij ;
        auto Real              ar[3], ari[3], s[MAXCBF*MAXCBF], t[MAXCBF*MAXCBF], xo[MAXAMP1*MAXAMP3], yo[MAXAMP1*MAXAMP3], zo[MAXAMP1*MAXAMP3],
                                                                                    xt[MAXAMP1*MAXAMP1], yt[MAXAMP1*MAXAMP1], zt[MAXAMP1*MAXAMP1] ;
        auto Real             *ri, *rj ;
        auto Integer                 i, iammax, ic, icbfind, iqm, ip, ishell, istart, ixo, ixt, iyo, iyt, izo, izt,
                                 j, jammax, jc, jcbfind, jdimo, jdimt, jp, jqm, jshell, jstart, jupper, jxixo, jyiyo, jzizo, jxixt, jyiyt, jzizt, n, ncfunci, ncfuncj ;
        auto GaussianBasis *ibasis, *jbasis ;

        /* . Initialization. */
        (*kinetic) = SymmetricMatrix_Allocate ( qcAtoms->nobasisw ) ;
        (*overlap) = SymmetricMatrix_Allocate ( qcAtoms->nobasisw ) ;
        SymmetricMatrix_Set_Zero ( (*kinetic) ) ;
        SymmetricMatrix_Set_Zero ( (*overlap) ) ;

        /* . Outer loop over centers. */
        for ( iqm = 0 ; iqm < qcAtoms->natoms ; iqm++ )
        {

            /* . Get data for the center. */
            ic     = qcAtoms->data[iqm].center ;
            ibasis = qcParameters->centers[ic].orbitalbasis ;
            istart = qcAtoms->data[iqm].ostartw ;
            ri     = Coordinates3_RowPointer ( qccoordinates3, iqm ) ;

            /* . Inner loop over centers. */
            for ( jqm = 0 ; jqm <= iqm ; jqm++ )
            {

                /* . Get data for the center. */
                jc     = qcAtoms->data[jqm].center ;
                jbasis = qcParameters->centers[jc].orbitalbasis ;
                jstart = qcAtoms->data[jqm].ostartw ;
                rj     = Coordinates3_RowPointer ( qccoordinates3, jqm ) ;

                /* . Set the diagonal atom flag. */
                QIQMJQM = ( iqm == jqm ) ;

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

                    /* . Set the upper limit for the JSHELL loops. */
                    if ( QIQMJQM ) jupper = ishell + 1 ;
                    else           jupper = jbasis->nshells ;

                    /* . Inner loop over shells. */
                    for ( jshell = 0 ; jshell < jupper ; jshell++ )
                    {

                        /* . Get information about the shell. */
                        jammax  = jbasis->shells[jshell].type->angularmomentum_high ;
                        jdimo   = jammax + 3 ;
                        jdimt   = jammax + 1 ;
                        jcbfind = jbasis->shells[jshell].type->cbfindex ;
                        ncfuncj = jbasis->shells[jshell].type->ncbf     ;

                        /* . Set the diagonal block flag. */
                        QDIAGONAL = QIQMJQM && ( ishell == jshell ) ;

                        /* . Initialize the integral blocks. */
                        for ( i = 0 ; i < ( ncfunci * ncfuncj ) ; i++ )
                        {
                            s[i] = 0.0e+00 ;
                            t[i] = 0.0e+00 ;
                        }

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

                                /* . Calculate the subsidiary integrals. */
                                Subsidiary_Integral_Overlap2 ( xo, yo, zo, aa, ar, ri, rj, iammax, jammax + 2 ) ;
                                Subsidiary_Integral_Kinetic  ( xo, yo, zo, xt, yt, zt, aj, iammax, jammax, jdimo, jdimt ) ;

                                /* . Add in the contributions to the full integrals. */
                                for ( i = 0, n = 0 ; i < ncfunci ; i++ )
                                {
   	                            ixo = CBFPOWX[i+icbfind] * jdimo ;
	                            iyo = CBFPOWY[i+icbfind] * jdimo ;
	                            izo = CBFPOWZ[i+icbfind] * jdimo ;
   	                            ixt = CBFPOWX[i+icbfind] * jdimt ;
	                            iyt = CBFPOWY[i+icbfind] * jdimt ;
	                            izt = CBFPOWZ[i+icbfind] * jdimt ;
                                    ti = expfac * ibasis->shells[ishell].primitives[ip].ccbf[i] ;
                                    for ( j = 0 ; j < ncfuncj ; j++, n++ )
                                    {
	                                jxixo = CBFPOWX[j+jcbfind] + ixo ;
	                                jyiyo = CBFPOWY[j+jcbfind] + iyo ;
	                                jzizo = CBFPOWZ[j+jcbfind] + izo ;
	                                jxixt = CBFPOWX[j+jcbfind] + ixt ;
	                                jyiyt = CBFPOWY[j+jcbfind] + iyt ;
	                                jzizt = CBFPOWZ[j+jcbfind] + izt ;
                                        tij  = ti * jbasis->shells[jshell].primitives[jp].ccbf[j] ;
    		                        s[n] += tij * xo[jxixo] * yo[jyiyo] * zo[jzizo] ;
                                        t[n] += tij * ( xt[jxixt] * yo[jyiyo] * zo[jzizo] +
                                                        xo[jxixo] * yt[jyiyt] * zo[jzizo] +
                                                        xo[jxixo] * yo[jyiyo] * zt[jzizt] ) ;
                                    }
                                }
                            }
                        }

                        /* . Transform the integrals. */
                        /* if ( qcAtoms->QTOSPHERICAL )
                        {
	                    Integral_Block_Transform_M ( s, iangmom, ncfunci, nfunci, jangmom, ncfuncj, nfuncj ) ;
	                    Integral_Block_Transform_M ( t, iangmom, ncfunci, nfunci, jangmom, ncfuncj, nfuncj ) ;
                        }
                        */

                        /* . Put the integrals into the proper place. */
                        if ( QDIAGONAL )
                        {
                            SymmetricMatrix_Set_DBlockM ( (*overlap), istart + ibasis->shells[ishell].nstartw, ibasis->shells[ishell].nbasisw, s ) ;
                            SymmetricMatrix_Set_DBlockM ( (*kinetic), istart + ibasis->shells[ishell].nstartw, ibasis->shells[ishell].nbasisw, t ) ;
                        }
                        else
                        {
                            SymmetricMatrix_Set_OBlockM ( (*overlap), istart + ibasis->shells[ishell].nstartw, ibasis->shells[ishell].nbasisw, jstart + jbasis->shells[jshell].nstartw, jbasis->shells[jshell].nbasisw, s ) ;
                            SymmetricMatrix_Set_OBlockM ( (*kinetic), istart + ibasis->shells[ishell].nstartw, ibasis->shells[ishell].nbasisw, jstart + jbasis->shells[jshell].nstartw, jbasis->shells[jshell].nbasisw, t ) ;
                        }
                    }
                }
            }
        }
# ifdef PRINTDIAGONALELEMENTS
{
auto Integer i ;
printf ( "\nDIAGONAL OVERLAP ELEMENTS:\n" ) ;
for ( i = 0 ; i < qcAtoms->nobasisw ; i++ )
{
    printf ( "OV> %10d %25.15f %25.15f\n", i, SymmetricMatrix_Get_Component ( (*overlap), i, i ), SymmetricMatrix_Get_Component ( (*overlap), i, i ) - 1.0e+00 ) ;
}
}
# endif
    }
}

/*------------------------------------------------------------------------------
! . Calculate the nuclear-nuclear terms.
!-----------------------------------------------------------------------------*/
double GaussianBasis_Nuclear_Nuclear ( const QCAtomContainer *qcAtoms, Coordinates3 *qccoordinates3 )
{
    Real enuclear = 0.0e+00 ;
    if ( ( qcAtoms != NULL ) && ( qccoordinates3 != NULL ) )
    {
        auto Integer                 iqm, jqm ;
        auto Real              iandj, ij, rho, rij2, xij, yij, zij, zi, zj ;
        auto Real             *ri, *rj ;
        auto QCAtom        *idata, *jdata ;
        auto RysQuadrature  roots ;

        /* . Outer loop over centers. */
        for ( iqm = 0 ; iqm < qcAtoms->natoms ; iqm++ )
        {

            /* . Get data for the center. */
            idata = &qcAtoms->data[iqm] ;
            ri    = Coordinates3_RowPointer ( qccoordinates3, iqm ) ;
            zi    = ( Real ) idata->atomicNumber ;

            /* . Inner loop over centers. */
            for ( jqm = 0 ; jqm < iqm ; jqm++ )
            {

                /* . Get data for the center. */
                jdata = &qcAtoms->data[jqm] ;
                rj    = Coordinates3_RowPointer ( qccoordinates3, jqm ) ;
                zj    = ( Real ) jdata->atomicNumber ;

                /* . Calculate some distance factors. */
	        xij  = ri[0] - rj[0] ;
	        yij  = ri[1] - rj[1] ;
	        zij  = ri[2] - rj[2] ;
	        rij2 = xij * xij + yij * yij + zij * zij ;

                /* . Nuclear repulsion. */
                iandj  = idata->widthe + jdata->widthe ;
                ij     = idata->widthe * jdata->widthe ;
                rho    = ij / iandj ;
                RysQuadrature_Roots ( &roots, 1, ( rho * rij2 ) ) ;
                enuclear += ( zi * zj * PI252 * idata->widthn * jdata->widthn * roots.weights[0] / ( ij * sqrt ( iandj ) ) ) ;
            }
        }
    }
    return enuclear ;
}

/*------------------------------------------------------------------------------
! . Calculate the electron-point integrals.
! . These procedures need to be rationalized with the nuclear ones!
!-----------------------------------------------------------------------------*/
/* . Point charge distribution parameters - like nuclei. */
# define POINT_WIDTH  1.0e-4
# define POINT_WIDTHE ( 4.0e+00 / ( POINT_WIDTH * POINT_WIDTH * M_PI ) )
# define POINT_WIDTHN sqrt ( pow ( ( POINT_WIDTHE / M_PI ), 3 ) )
void GaussianBasis_Point_Electron ( const QCAtomContainer *qcAtoms, const QCParameter *qcParameters, const Coordinates3 *qccoordinates3, const Coordinates3 *points,
                                                                                                         const SymmetricMatrix *odensityt, Real1DArray *potentials )
{
    if ( ( qcAtoms != NULL ) && ( qcParameters != NULL ) && ( qccoordinates3 != NULL ) && ( points != NULL ) && ( odensityt != NULL ) && ( potentials != NULL ) && ( points->length0 <= potentials->length ) )
    {
        auto Boolean                QDIAGONAL, QIJ0, QIJ1, QIQMJQM ;
        auto Real              aa, aandb, aainv, ab, ai, aj, arri, axac, ayac, azac, bp01, b00, b10, c1x, c1y, c1z, c3x, c3y, c3z, c4x, c4y, c4z,
                                 dnuc, dxijt, dyijt, dzijt, expfac, expn, fac, facn, fac2, f00, pot, rho, rij2, scale, ti, tij, u2,
                                 xc00, xcp00, xij, yc00, ycp00, yij, zc00, zcp00, zij ;
        auto Real              ar[3], ari[3], g[MAXCBF*MAXCBF], xint[MAXAMP1*MAXAMP1*MAXRYS], yint[MAXAMP1*MAXAMP1*MAXRYS], zint[MAXAMP1*MAXAMP1*MAXRYS] ;
        auto Real             *rc, *ri, *rj, *rn ;
        auto Integer                 i, iammax, iammaxt, ic, icbfind, ii, ij, iqm, ip, ishell, istart, ix, iy, iz,
                                 j, jammax, jammaxt, jc, jcbfind, jdim, jdimm, jj, jp, jqm, jshell, jstart,
                                 jupper, jxix, jyiy, jziz, k, m, n, ncfunci, ncfuncj, nroots ;
        auto GaussianBasis *ibasis, *jbasis ;
        auto RysQuadrature  roots ;

        /* . Loop over the points. */
        for ( k = 0 ; k < points->length0 ; k++ )
        {
            /* . Get information about the atom. */
            expn = POINT_WIDTHE ;
            facn = POINT_WIDTHN ;
            rn   = Coordinates3_RowPointer ( points, k ) ;

            /* . Initialize potential. */
            pot = 0.0e+00 ;

            /* . Outer loop over centers. */
            for ( iqm = 0 ; iqm < qcAtoms->natoms ; iqm++ )
            {

                /* . Get data for the center. */
                ic     = qcAtoms->data[iqm].center ;
                ibasis = qcParameters->centers[ic].orbitalbasis ;
                istart = qcAtoms->data[iqm].ostartw ;
                ri     = Coordinates3_RowPointer ( qccoordinates3, iqm ) ;

                /* . Inner loop over centers. */
                for ( jqm = 0 ; jqm <= iqm ; jqm++ )
                {

                    /* . Get data for the center. */
                    jc     = qcAtoms->data[jqm].center ;
                    jbasis = qcParameters->centers[jc].orbitalbasis ;
                    jstart = qcAtoms->data[jqm].ostartw ;
                    rj     = Coordinates3_RowPointer ( qccoordinates3, jqm ) ;

                    /* . Set the diagonal atom flag. */
                    QIQMJQM = ( iqm == jqm ) ;

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

                        /* . Set the upper limit for the JSHELL loops. */
                        if ( QIQMJQM ) jupper = ishell + 1 ;
                        else           jupper = jbasis->nshells ;

                        /* . Inner loop over shells. */
                        for ( jshell = 0 ; jshell < jupper ; jshell++ )
                        {

                            /* . Get information about the shell. */
                            jammax  = jbasis->shells[jshell].type->angularmomentum_high ;
                            jdimm   = ( iammax + 1 ) * ( jammax + 1 ) ;
                            jcbfind = jbasis->shells[jshell].type->cbfindex ;
                            ncfuncj = jbasis->shells[jshell].type->ncbf     ;

                            /* . Set the diagonal block flag. */
                            QDIAGONAL = QIQMJQM && ( ishell == jshell ) ;

                            /* . Get the number of roots. */
                            nroots = ( iammax + jammax ) / 2 + 1 ;

                            /* . Set some flags. */
                            QIJ0 = ( iammax + jammax == 0 ) ;
                            QIJ1 = ( iammax + jammax <= 1 ) ;

                            /* . Initialize the integral blocks. */
                            for ( i = 0 ; i < ( ncfunci * ncfuncj ) ; i++ ) g[i] = 0.0e+00 ;

                            /* . Select the expansion center for the recurrence relations. */
                            if ( iammax >= jammax )
                            {
                               iammaxt = iammax ;
                               jammaxt = jammax ;
                               jdim    = jammax + 1 ;
                               dxijt   = xij ;
                               dyijt   = yij ;
                               dzijt   = zij ;
                               rc      = ri ;
                            }
                            else
                            {
                               iammaxt = jammax ;
                               jammaxt = iammax ;
                               jdim    = iammax + 1 ;
                               dxijt   = - xij ;
                               dyijt   = - yij ;
                               dzijt   = - zij ;
                               rc      = rj ;
                            }

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
                                    expfac = exp ( - fac ) * PI252 * aainv ;
	                            for ( i = 0 ; i < 3 ; i++ ) ar[i] = ( ari[i] + aj * rj[i] ) * aainv ;

                                    /* . Start of point-specific code. */
                                    /* . Calculate some factors. */
                                    ab     = aa * expn ;
                                    aandb  = aa + expn ;
                                    rho    = ab / aandb ;
                                    dnuc   = expfac * facn / ( expn * sqrt ( aandb ) ) ;

                                    /* . Calculate the rys polynomial roots. */
                                    c1x = ( ar[0] - rn[0] ) ;
                                    c1y = ( ar[1] - rn[1] ) ;
                                    c1z = ( ar[2] - rn[2] ) ;
                                    RysQuadrature_Roots ( &roots, nroots, rho * ( c1x * c1x + c1y * c1y + c1z * c1z ) ) ;

                                    /* . Calculate some displacements. */
                                    axac = aa * ( ar[0] - rc[0] ) ;
                                    ayac = aa * ( ar[1] - rc[1] ) ;
                                    azac = aa * ( ar[2] - rc[2] ) ;
                                    c1x *= aa ;
                                    c1y *= aa ;
                                    c1z *= aa ;
                                    c3x  = expn * ( rn[0] - rc[0] ) + axac ;
                                    c3y  = expn * ( rn[1] - rc[1] ) + ayac ;
                                    c3z  = expn * ( rn[2] - rc[2] ) + azac ;
                                    c4x  = expn * axac ;
                                    c4y  = expn * ayac ;
                                    c4z  = expn * azac ;

                                    /* . Loop over the roots and construct the subsidiary integrals. */
                                    for ( m = 0 ; m < nroots ; m++ )
                                    {
                                        u2    = roots.roots[m] * rho ;
                                        f00   = roots.weights[m] ;
                                        fac   = 1.0e+00 / ( ab + u2 * aandb ) ;
                                        fac2  = 0.5e+00 * fac ;
                                        bp01  = ( aa   + u2 ) * fac2 ;
                                        b00   =          u2   * fac2 ;
                                        b10   = ( expn + u2 ) * fac2 ;
                                        xcp00 = u2 * c1x * fac ;
                                        ycp00 = u2 * c1y * fac ;
                                        zcp00 = u2 * c1z * fac ;
                                        xc00  = ( u2 * c3x + c4x ) * fac ;
                                        yc00  = ( u2 * c3y + c4y ) * fac ;
                                        zc00  = ( u2 * c3z + c4z ) * fac ;
                                        Subsidiary_Integral_Nuclear3C ( iammaxt, jammaxt, 0, QIJ0, QIJ1, True, True, b00, b10, bp01, dxijt, dyijt,dzijt, f00,
                                                                        xc00, xcp00, yc00, ycp00, zc00, zcp00, 1, jdim, &xint[m*jdimm], &yint[m*jdimm], &zint[m*jdimm] ) ;
                                    }

                                    /* . Assemble the integrals. */
                                    if ( iammax >= jammax )
                                    {
                                        for ( i = 0, n = 0 ; i < ncfunci ; i++ )
                                        {
   	                                    ix = CBFPOWX[i+icbfind] * jdim ;
	                                    iy = CBFPOWY[i+icbfind] * jdim ;
	                                    iz = CBFPOWZ[i+icbfind] * jdim ;
                                            ti = dnuc * ibasis->shells[ishell].primitives[ip].ccbf[i] ;
                                            for ( j = 0 ; j < ncfuncj ; j++, n++ )
                                            {
	                                        jxix = CBFPOWX[j+jcbfind] + ix ;
	                                        jyiy = CBFPOWY[j+jcbfind] + iy ;
	                                        jziz = CBFPOWZ[j+jcbfind] + iz ;
                                                for ( m = 0, fac = 0.0e+00 ; m < nroots ; m++ ) fac += xint[jxix+m*jdimm] * yint[jyiy+m*jdimm] * zint[jziz+m*jdimm] ;
                                                tij = ti * jbasis->shells[jshell].primitives[jp].ccbf[j] ;
                                                g[n] += tij * fac ;
                                            }
                                        }
                                    }
                                    else
                                    {
                                        for ( j = 0 ; j < ncfuncj ; j++ )
                                        {
   	                                    ix = CBFPOWX[j+jcbfind] * jdim ;
	                                    iy = CBFPOWY[j+jcbfind] * jdim ;
	                                    iz = CBFPOWZ[j+jcbfind] * jdim ;
                                            ti = dnuc * jbasis->shells[jshell].primitives[jp].ccbf[j] ;
                                            for ( i = 0, n = j ; i < ncfunci ; i++, n+= ncfuncj )
                                            {
	                                        jxix = CBFPOWX[i+icbfind] + ix ;
	                                        jyiy = CBFPOWY[i+icbfind] + iy ;
	                                        jziz = CBFPOWZ[i+icbfind] + iz ;
                                                for ( m = 0, fac = 0.0e+00 ; m < nroots ; m++ ) fac += xint[jxix+m*jdimm] * yint[jyiy+m*jdimm] * zint[jziz+m*jdimm] ;
                                                tij = ti * ibasis->shells[ishell].primitives[ip].ccbf[i] ;
                                                g[n] += tij * fac ;
                                            }
                                        }
                                    }
                                }
                            }


                            /* . Transform the integrals. */
	                    /* if ( qcAtoms->QTOSPHERICAL ) Integral_Block_Transform_M ( g, iangmom, ncfunci, nfunci, jangmom, ncfuncj, nfuncj ) ; */

                            /* . Determine scaling. */
                            if ( QDIAGONAL ) scale = 1.0e+00 ;
                            else             scale = 2.0e+00 ;

                            /* .  Add in the block of integrals to the potential - i is usually greater than j. */
                            for ( i = 0, n = 0 ; i < ibasis->shells[ishell].nbasisw ; i++ )
                            {
                                ii = istart + ibasis->shells[ishell].nstartw + i ;
                                for ( j = 0 ; j < jbasis->shells[jshell].nbasisw ; j++, n++ )
                                {
                                    jj = jstart + jbasis->shells[jshell].nstartw + j ;
                                    if ( ii >= jj ) ij = ( ii * ( ii + 1 ) ) / 2 + jj ;
                                    else            ij = ( jj * ( jj + 1 ) ) / 2 + ii ;
                                    pot += odensityt->data[ij] * g[n] * scale ;
	                        }
                            }
                        }
                    }
                }
            }
            /* . Sum in the contribution to the potential. */
            Real1DArray_Item ( potentials, k ) -= pot ;
        }
    }
}

/*------------------------------------------------------------------------------
! . Calculate nuclear-point potentials.
! . Potentials should be initialized before entry!
!-----------------------------------------------------------------------------*/
void GaussianBasis_Point_Nuclear ( const QCAtomContainer *qcAtoms, const Real1DArray *znuclear, const Coordinates3 *qccoordinates3, const Coordinates3 *points, Real1DArray *potentials )
{
    if ( ( qcAtoms != NULL ) && ( znuclear != NULL ) && ( qccoordinates3 != NULL ) && ( points != NULL ) && ( potentials != NULL ) && ( points->length0 <= potentials->length ) )
    {
        auto Integer                 iqm, j ;
        auto Real              iandj, ij, rho, rij2, xij, yij, zij, zi ;
        auto Real             *ri, *rj ;
        auto QCAtom        *idata ;
        auto RysQuadrature  roots ;

        /* . Loop over centers. */
        for ( iqm = 0 ; iqm < qcAtoms->natoms ; iqm++ )
        {

            /* . Get data for the center. */
            idata = &qcAtoms->data[iqm] ;
            ri    = Coordinates3_RowPointer ( qccoordinates3, iqm ) ;
            zi    = Real1DArray_Item ( znuclear, iqm ) ;

            /* . Loop over points. */
            for ( j = 0 ; j < points->length0 ; j++ )
            {
                /* . Get data for the center. */
                rj = Coordinates3_RowPointer ( points, j ) ;

                /* . Calculate some distance factors. */
	        xij  = ri[0] - rj[0] ;
	        yij  = ri[1] - rj[1] ;
	        zij  = ri[2] - rj[2] ;
	        rij2 = xij * xij + yij * yij + zij * zij ;

                /* . Potential. */
                iandj  = idata->widthe + POINT_WIDTHE ;
                ij     = idata->widthe * POINT_WIDTHE ;
                rho    = ij / iandj ;
                RysQuadrature_Roots ( &roots, 1, ( rho * rij2 ) ) ;
                Real1DArray_Item ( potentials, j ) += ( zi * PI252 * idata->widthn * POINT_WIDTHN * roots.weights[0] / ( ij * sqrt ( iandj ) ) ) ;
            }
        }
    }
}
# undef POINT_WIDTH
# undef POINT_WIDTHE
# undef POINT_WIDTHN

/*------------------------------------------------------------------------------
! . Calculate the selfoverlap integrals for a basis.
! . All zero or even 1-D polynomials are non-zero in Cartesians.
! . 1-D integrals for x^n (n even) are (n-1)!!/(2 a)^(n/2) (Pi/a)^(1/2).
! . For spherical harmonics all zero except for s-functions.
!-----------------------------------------------------------------------------*/
void GaussianBasis_SelfOverlap ( const QCAtomContainer *qcAtoms, const QCParameter *qcParameters, Real1DArray **fselfoverlap )
{
    if ( ( qcAtoms != NULL ) && ( qcParameters != NULL ) )
    {
        auto Real              ci, ei, si ;
        auto Integer                 i, ic, icbfind, ifstart, ip, iqm, ishell, ix, iy, iz, ncfunci, t ;
        auto GaussianBasis *ibasis ;

        /* . Initialization. */
        (*fselfoverlap) = Real1DArray_Allocate ( qcAtoms->nfbasisw, NULL ) ;
        Real1DArray_Set ( (*fselfoverlap), 0.0e+00 ) ;
/*printf ( "\nSelf Overlaps:\n" ) ;*/
        /* . Loop over atoms. */
        for ( iqm = 0 ; iqm < qcAtoms->natoms ; iqm++ )
        {
            ic      = qcAtoms->data[iqm].center  ;
            ifstart = qcAtoms->data[iqm].fstartw ;
            ibasis  = qcParameters->centers[ic].densitybasis ;
            /* . Loop over shells. */
            for ( ishell = 0 ; ishell < ibasis->nshells ; ishell++ )
            {
                /* . Get information about the shell. */
                icbfind = ibasis->shells[ishell].type->cbfindex ;
                ncfunci = ibasis->shells[ishell].type->ncbf     ;
                /* . Loop over functions in the shell. */
                for ( i = 0 ; i < ncfunci ; i++ )
                {
   	            ix = CBFPOWX[i+icbfind] ;
	            iy = CBFPOWY[i+icbfind] ;
	            iz = CBFPOWZ[i+icbfind] ;
                    /* . Zero or even polynomials only. */
                    if ( IsEven ( ix ) && IsEven ( iy ) && IsEven ( iz ) )
                    {
                        si = 0.0e+00 ;
                        for ( ip = 0 ; ip < ibasis->shells[ishell].nprimitives ; ip++ )
                        {
                            ci = ibasis->shells[ishell].primitives[ip].ccbf[i]  ;
                            ei = ibasis->shells[ishell].primitives[ip].exponent ;
                            for ( t = 1 ; t <= ( ix / 2 ) ; t++ ) ci *= ( Real ) ( 2 * t - 1 ) / ( 2.0e+00 * ei ) ;
                            for ( t = 1 ; t <= ( iy / 2 ) ; t++ ) ci *= ( Real ) ( 2 * t - 1 ) / ( 2.0e+00 * ei ) ;
                            for ( t = 1 ; t <= ( iz / 2 ) ; t++ ) ci *= ( Real ) ( 2 * t - 1 ) / ( 2.0e+00 * ei ) ;
                            si += ci / ( ei * sqrt ( ei ) ) ;
                        }
                        (*fselfoverlap)->data[ifstart+ibasis->shells[ishell].nstartw+i] = PI32 * si ;
                    }
/*printf ( "%10d %10d %10d %10d %25.15f\n", ifstart+ibasis->shells[ishell].nstartw+i, ix, iy, iz, (*fselfoverlap)->data[ifstart+ibasis->shells[ishell].nstartw+i] ) ;*/
                }
            }
        }
    }
}
