/*------------------------------------------------------------------------------
! . File      : GaussianBasisDerivatives.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==============================================================================
! . Procedures for calculating the derivatives of integrals over Gaussians.
!=============================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "GaussianBasisDerivatives.h"
# include "GaussianBasisSubsidiary.h"
# include "RysQuadrature.h"
# include "Memory.h"

/*------------------------------------------------------------------------------
! . Calculate the electron-fit integrals.
!-----------------------------------------------------------------------------*/
void GaussianBasis_Electron_FitD ( const QCAtomContainer *qcAtoms, const QCParameter *qcParameters, Coordinates3 *qccoordinates3,
                                           SymmetricMatrix *odensityt, const Real1DArray *fpotential, const Real1DArray *wvector,
                                                                                                     Coordinates3 *qcgradients3 )
{
    if ( ( qcAtoms != NULL ) && ( qcParameters != NULL ) && ( qccoordinates3 != NULL ) )
    {
        auto Boolean                  QDIAGONAL, QF0, QF1, QIQMJQM ;
        auto Real                aa, aandb, aainv, ab, ag, ah, ai, aj, arri, axac, ayac, azac, bp01, b00, b10, c1x, c1y, c1z, c3x, c3y, c3z, c4x, c4y, c4z,
                                   dgx, dgy, dgz, dhx, dhy, dhz, dnuc, dxijt, dyijt, dzijt, expfac, expf, fac, facf, facgx, facgy, facgz, fachx, fachy, fachz,
                                   faci, facj, fac2, f00, rho, rfi2, rfj2, rij2, scale, ti, tij, tijf, u2, xc00, xcp00, xij, yc00, ycp00, yij, zc00, zcp00, zij ;
        auto Real                ar[3], ari[3],
                                   gx[MAXCBF*MAXCBF*MAXCBF], gy[MAXCBF*MAXCBF*MAXCBF], gz[MAXCBF*MAXCBF*MAXCBF],
                                   hx[MAXCBF*MAXCBF*MAXCBF], hy[MAXCBF*MAXCBF*MAXCBF], hz[MAXCBF*MAXCBF*MAXCBF],
                                   xidg[MAXAMP1*MAXAMP1*MAXAMP1*MAXRYS], yidg[MAXAMP1*MAXAMP1*MAXAMP1*MAXRYS], zidg[MAXAMP1*MAXAMP1*MAXAMP1*MAXRYS],
                                   xidh[MAXAMP1*MAXAMP1*MAXAMP1*MAXRYS], yidh[MAXAMP1*MAXAMP1*MAXAMP1*MAXRYS], zidh[MAXAMP1*MAXAMP1*MAXAMP1*MAXRYS],
                                   xint[MAXAMP1*MAXAMP2*MAXAMP2*MAXRYS], yint[MAXAMP1*MAXAMP2*MAXAMP2*MAXRYS], zint[MAXAMP1*MAXAMP2*MAXAMP2*MAXRYS] ;
        auto Real               *rc, *rf, *ri, *rj ;
        auto Integer                   dim1, ddim2, dim2, ddim3, dim3,
                                   f, fammax,          fc, fcbfind, fijx, fijxd, fijy, fijyd, fijz, fijzd, fqm, fp, fshell, fstart,
                                   i, iammax, iammaxt, ic, icbfind, ii, ij, iqm, ip, ishell, istart, ix, ixd, iy, iyd, iz, izd,
                                   j, jammax, jammaxt, jc, jcbfind, jix, jixd, jiy, jiyd, jiz, jizd, jj, jp, jqm, jshell, jstart,
                                   jupper, m, n, ncfuncf, ncfunci, ncfuncj, nroots ;
        auto GaussianBasis   *ibasis, *jbasis, *fbasis ;
        auto RysQuadrature    roots ;

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
                    for ( i = 0, rfi2 = rfj2 = 0.0e+00 ; i < 3 ; i++ )
                    {
                        faci  = rf[i] - ri[i] ;
                        facj  = rf[i] - rj[i] ;
                        rfi2 += faci * faci ;
                        rfj2 += facj * facj ;
                    }

                    /* . Initialize some accumulators. */
                    dgx = dgy = dgz = dhx = dhy = dhz = 0.0e+00 ;

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
                                nroots = ( fammax + iammax + jammax + 2 ) / 2 + 1 ;

                                /* . Initialize the integral blocks. */
                                for ( i = 0 ; i < ( ncfuncf * ncfunci * ncfuncj ) ; i++ ) gx[i] = gy[i] = gz[i] = hx[i] = hy[i] = hz[i] = 0.0e+00 ;

                                /* . Set the array dimensions. */
                                dim1 = fammax + 1 ;
                                if ( iammax >= jammax )
                                {
                                    ddim2 = dim1 * ( jammax + 1 ) ;
                                    dim2  = dim1 * ( jammax + 2 ) ;
                                }
                                else
                                {
                                    ddim2 = dim1 * ( iammax + 1 ) ;
                                    dim2  = dim1 * ( iammax + 2 ) ;
                                }
                                ddim3 = dim1 * ( iammax + 1 ) * ( jammax + 1 ) ;
                                dim3  = dim1 * ( iammax + 2 ) * ( jammax + 2 ) ;

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

                                        /* . Set the primitive exponents. */
                                        if ( iammax >= jammax )
                                        {
                                           ag = ai ;
                                           ah = aj ;
                                        }
                                        else
                                        {
                                           ag = aj ;
                                           ah = ai ;
                                        }

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
                                                Subsidiary_Integral_Nuclear3C ( iammaxt+1, jammaxt+1, fammax, False, False, QF0, QF1, b00, b10, bp01, dxijt, dyijt, dzijt, f00,
                                                                                xc00, xcp00, yc00, ycp00, zc00, zcp00, dim1, dim2, &xint[m*dim3], &yint[m*dim3], &zint[m*dim3] ) ;
                                                Subsidiary_Integral_Derivative3 ( &xint[m*dim3], &yint[m*dim3], &zint[m*dim3],
                                                                                  &xidg[m*ddim3], &yidg[m*ddim3], &zidg[m*ddim3],
                                                                                  &xidh[m*ddim3], &yidh[m*ddim3], &zidh[m*ddim3],
                                                                                  ag, ah, iammaxt, jammaxt, fammax, dim1, dim2, dim1, ddim2 ) ;
                                            }

                                            /* . Assemble the integrals. */
                                            if ( iammax >= jammax )
                                            {
                                                for ( i = 0, n = 0 ; i < ncfunci ; i++ )
                                                {
   	                                            ix  = CBFPOWX[i+icbfind] * dim2  ;
	                                            iy  = CBFPOWY[i+icbfind] * dim2  ;
	                                            iz  = CBFPOWZ[i+icbfind] * dim2  ;
   	                                            ixd = CBFPOWX[i+icbfind] * ddim2 ;
	                                            iyd = CBFPOWY[i+icbfind] * ddim2 ;
	                                            izd = CBFPOWZ[i+icbfind] * ddim2 ;
                                                    ti  = dnuc * ibasis->shells[ishell].primitives[ip].ccbf[i] ;
                                                    for ( j = 0 ; j < ncfuncj ; j++ )
                                                    {
	                                                jix  = CBFPOWX[j+jcbfind] * dim1 + ix  ;
	                                                jiy  = CBFPOWY[j+jcbfind] * dim1 + iy  ;
	                                                jiz  = CBFPOWZ[j+jcbfind] * dim1 + iz  ;
	                                                jixd = CBFPOWX[j+jcbfind] * dim1 + ixd ;
	                                                jiyd = CBFPOWY[j+jcbfind] * dim1 + iyd ;
	                                                jizd = CBFPOWZ[j+jcbfind] * dim1 + izd ;
                                                        tij  = ti * jbasis->shells[jshell].primitives[jp].ccbf[j] ;
                                                        for ( f = 0 ; f < ncfuncf ; f++, n++ )
                                                        {
                                                            fijx  = CBFPOWX[f+fcbfind] + jix  ;
                                                            fijy  = CBFPOWY[f+fcbfind] + jiy  ;
                                                            fijz  = CBFPOWZ[f+fcbfind] + jiz  ;
                                                            fijxd = CBFPOWX[f+fcbfind] + jixd ;
                                                            fijyd = CBFPOWY[f+fcbfind] + jiyd ;
                                                            fijzd = CBFPOWZ[f+fcbfind] + jizd ;
                                                            for ( m = 0, facgx = facgy = facgz = fachx = fachy = fachz = 0.0e+00 ; m < nroots ; m++ )
                                                            {
                                                                facgx += xidg[fijxd+m*ddim3] * yint[fijy+m*dim3] * zint[fijz+m*dim3] ;
                                                                facgy += xint[fijx+m*dim3] * yidg[fijyd+m*ddim3] * zint[fijz+m*dim3] ;
                                                                facgz += xint[fijx+m*dim3] * yint[fijy+m*dim3] * zidg[fijzd+m*ddim3] ;
                                                                fachx += xidh[fijxd+m*ddim3] * yint[fijy+m*dim3] * zint[fijz+m*dim3] ;
                                                                fachy += xint[fijx+m*dim3] * yidh[fijyd+m*ddim3] * zint[fijz+m*dim3] ;
                                                                fachz += xint[fijx+m*dim3] * yint[fijy+m*dim3] * zidh[fijzd+m*ddim3] ;
                                                            }
                                                            tijf = tij * fbasis->shells[fshell].primitives[fp].ccbf[f] ;
                                                            gx[n] += tijf * facgx ;
                                                            gy[n] += tijf * facgy ;
                                                            gz[n] += tijf * facgz ;
                                                            hx[n] += tijf * fachx ;
                                                            hy[n] += tijf * fachy ;
                                                            hz[n] += tijf * fachz ;
/*printf ( "\na %6d%6d%6d%6d%6d%6d%25.15f%25.15f%25.15f%25.15f%25.15f", ishell, jshell, fshell, i, j, f, gx[n], hx[n], tijf, facgx, fachx ) ;*/
                                                        }
                                                    }
                                                }
                                            }
                                            else
                                            {
                                                for ( j = 0 ; j < ncfuncj ; j++ )
                                                {
   	                                            ix  = CBFPOWX[j+jcbfind] * dim2  ;
	                                            iy  = CBFPOWY[j+jcbfind] * dim2  ;
	                                            iz  = CBFPOWZ[j+jcbfind] * dim2  ;
    	                                            ixd = CBFPOWX[j+jcbfind] * ddim2 ;
	                                            iyd = CBFPOWY[j+jcbfind] * ddim2 ;
	                                            izd = CBFPOWZ[j+jcbfind] * ddim2 ;
                                                    ti  = dnuc * jbasis->shells[jshell].primitives[jp].ccbf[j] ;
                                                    for ( i = 0 ; i < ncfunci ; i++ )
                                                    {
                                                        n    = j * ncfuncf + i * ( ncfuncf * ncfuncj ) ;
	                                                jix  = CBFPOWX[i+icbfind] * dim1 + ix  ;
	                                                jiy  = CBFPOWY[i+icbfind] * dim1 + iy  ;
	                                                jiz  = CBFPOWZ[i+icbfind] * dim1 + iz  ;
	                                                jixd = CBFPOWX[i+icbfind] * dim1 + ixd ;
	                                                jiyd = CBFPOWY[i+icbfind] * dim1 + iyd ;
	                                                jizd = CBFPOWZ[i+icbfind] * dim1 + izd ;
                                                        tij  = ti * ibasis->shells[ishell].primitives[ip].ccbf[i] ;
                                                        for ( f = 0 ; f < ncfuncf ; f++, n++ )
                                                        {
                                                            fijx  = CBFPOWX[f+fcbfind] + jix  ;
                                                            fijy  = CBFPOWY[f+fcbfind] + jiy  ;
                                                            fijz  = CBFPOWZ[f+fcbfind] + jiz  ;
                                                            fijxd = CBFPOWX[f+fcbfind] + jixd ;
                                                            fijyd = CBFPOWY[f+fcbfind] + jiyd ;
                                                            fijzd = CBFPOWZ[f+fcbfind] + jizd ;
                                                            for ( m = 0, facgx = facgy = facgz = fachx = fachy = fachz = 0.0e+00 ; m < nroots ; m++ )
                                                            {
                                                                facgx += xidh[fijxd+m*ddim3] * yint[fijy+m*dim3] * zint[fijz+m*dim3] ;
                                                                facgy += xint[fijx+m*dim3] * yidh[fijyd+m*ddim3] * zint[fijz+m*dim3] ;
                                                                facgz += xint[fijx+m*dim3] * yint[fijy+m*dim3] * zidh[fijzd+m*ddim3] ;
                                                                fachx += xidg[fijxd+m*ddim3] * yint[fijy+m*dim3] * zint[fijz+m*dim3] ;
                                                                fachy += xint[fijx+m*dim3] * yidg[fijyd+m*ddim3] * zint[fijz+m*dim3] ;
                                                                fachz += xint[fijx+m*dim3] * yint[fijy+m*dim3] * zidg[fijzd+m*ddim3] ;
                                                            }
                                                            tijf = tij * fbasis->shells[fshell].primitives[fp].ccbf[f] ;
                                                            gx[n] += tijf * facgx ;
                                                            gy[n] += tijf * facgy ;
                                                            gz[n] += tijf * facgz ;
                                                            hx[n] += tijf * fachx ;
                                                            hy[n] += tijf * fachy ;
                                                            hz[n] += tijf * fachz ;
/*printf ( "\nb %6d%6d%6d%6d%6d%6d%25.15f%25.15f%25.15f%25.15f%25.15f", ishell, jshell, fshell, i, j, f, gx[n], hx[n], tijf, facgx, fachx ) ;*/
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }

                                /* . Transform the integrals. */
                                /* if ( qcAtoms->QTOSPHERICAL )
                                {
   	                            Integral_Block_O3Transform ( gx, fangmom, ncfuncf, nsfuncf, iangmom, ncfunci, nfunci, jangmom, ncfuncj, nfuncj ) ;
   	                            Integral_Block_O3Transform ( gy, fangmom, ncfuncf, nsfuncf, iangmom, ncfunci, nfunci, jangmom, ncfuncj, nfuncj ) ;
   	                            Integral_Block_O3Transform ( gz, fangmom, ncfuncf, nsfuncf, iangmom, ncfunci, nfunci, jangmom, ncfuncj, nfuncj ) ;
   	                            Integral_Block_O3Transform ( hx, fangmom, ncfuncf, nsfuncf, iangmom, ncfunci, nfunci, jangmom, ncfuncj, nfuncj ) ;
   	                            Integral_Block_O3Transform ( hy, fangmom, ncfuncf, nsfuncf, iangmom, ncfunci, nfunci, jangmom, ncfuncj, nfuncj ) ;
   	                            Integral_Block_O3Transform ( hz, fangmom, ncfuncf, nsfuncf, iangmom, ncfunci, nfunci, jangmom, ncfuncj, nfuncj ) ;
                                }
                                */

                                /* . Get the scale factor. */
                                if ( QDIAGONAL ) scale = 1.0e+00 ;
                                else             scale = 2.0e+00 ;

                                /* . Determine the gradient contributions. */
                                for ( i = 0, n = 0 ; i < ibasis->shells[ishell].nbasisw ; i++ )
                                {
                                    ii = istart + ibasis->shells[ishell].nstartw + i ;
                                    for ( j = 0 ; j < jbasis->shells[jshell].nbasisw ; j++ )
                                    {
                                        jj = jstart + jbasis->shells[jshell].nstartw + j ;
                                        if ( ii >= jj ) ij = ( ii * ( ii + 1 ) ) / 2 + jj ;
                                        else            ij = ( jj * ( jj + 1 ) ) / 2 + ii ;
                                        fac = scale * odensityt->data[ij] ;
                                        for ( f = 0 ; f < fbasis->shells[fshell].nbasisw ; f++, n++ )
                                        {
                                            facf = fac * ( fpotential->data[fstart + fbasis->shells[fshell].nstartw + f] + wvector->data[fstart + fbasis->shells[fshell].nstartw + f] ) ;
                                            dgx += facf * gx[n] ;
                                            dgy += facf * gy[n] ;
                                            dgz += facf * gz[n] ;
                                            dhx += facf * hx[n] ;
                                            dhy += facf * hy[n] ;
                                            dhz += facf * hz[n] ;
                                        }
                                    }
                                }
                            }
                        }
                    }

                    /* . Sum in the contribution to the gradient. */
                    Coordinates3_IncrementRow ( qcgradients3, iqm, dgx, dgy, dgz ) ;
                    Coordinates3_IncrementRow ( qcgradients3, jqm, dhx, dhy, dhz ) ;
                    Coordinates3_DecrementRow ( qcgradients3, fqm, dgx, dgy, dgz ) ;
                    Coordinates3_DecrementRow ( qcgradients3, fqm, dhx, dhy, dhz ) ;
                }
            }
        }
    }
}

/*------------------------------------------------------------------------------
! . Calculate the electron-nuclear integrals.
!-----------------------------------------------------------------------------*/
void GaussianBasis_Electron_NuclearD ( const QCAtomContainer *qcAtoms, const QCParameter *qcParameters, Coordinates3 *qccoordinates3,
                                                                                  SymmetricMatrix *odensityt, Coordinates3 *qcgradients3 )
{
    if ( ( qcAtoms != NULL ) && ( qcParameters != NULL ) && ( qccoordinates3 != NULL ) && ( odensityt != NULL ) && ( qcgradients3 != NULL ) )
    {
        auto Boolean                  QDIAGONAL, QIQMJQM ;
        auto Real                aa, aandb, aainv, ab, ag, ah, ai, aj, arri, axac, ayac, azac, bp01, b00, b10, c1x, c1y, c1z, c3x, c3y, c3z, c4x, c4y, c4z,
                                   dgx, dgy, dgz, dhx, dhy, dhz, dnuc, dxijt, dyijt, dzijt, expfac, expn, fac, facgx, facgy, facgz, fachx, fachy, fachz,
                                   facn, fac2, f00, qn, rho, rij2, scale, ti, tij, u2,
                                   xc00, xcp00, xij, yc00, ycp00, yij, zc00, zcp00, zij ;
        auto Real                ar[3], ari[3],
                                   gx[MAXCBF*MAXCBF], gy[MAXCBF*MAXCBF], gz[MAXCBF*MAXCBF],
                                   hx[MAXCBF*MAXCBF], hy[MAXCBF*MAXCBF], hz[MAXCBF*MAXCBF],
                                   xidg[MAXAMP1*MAXAMP1*MAXRYS], yidg[MAXAMP1*MAXAMP1*MAXRYS], zidg[MAXAMP1*MAXAMP1*MAXRYS],
                                   xidh[MAXAMP1*MAXAMP1*MAXRYS], yidh[MAXAMP1*MAXAMP1*MAXRYS], zidh[MAXAMP1*MAXAMP1*MAXRYS],
                                   xint[MAXAMP2*MAXAMP2*MAXRYS], yint[MAXAMP2*MAXAMP2*MAXRYS], zint[MAXAMP2*MAXAMP2*MAXRYS] ;
        auto Real               *rc, *ri, *rj, *rn ;
        auto Integer                   ddim1, ddim2,
                                   i, iammax, iammaxt, ic, icbfind, ii, ij, iqm, ip, ishell, istart, ix, ixd, iy, iyd, iz, izd,
                                   j, jammax, jammaxt, jc, jcbfind, jdim, jdimm, jj, jp, jqm, jshell, jstart,
                                   jupper, jxix, jxixd, jyiy, jyiyd, jziz, jzizd, kqm, kc, m, n, ncfunci, ncfuncj, nroots ;
        auto GaussianBasis   *ibasis, *jbasis ;
        auto RysQuadrature    roots ;

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

                /* . Loop over the nuclear densities. */
                for ( kqm = 0 ; kqm < qcAtoms->natoms ; kqm++ )
                {

                    /* . Get information about the atom. */
                    kc   = qcAtoms->data[kqm].center ;
                    expn = qcAtoms->data[kqm].widthe ;
                    facn = qcAtoms->data[kqm].widthn ;
                    qn   = - ( Real ) qcAtoms->data[kqm].atomicNumber ;
                    rn   = Coordinates3_RowPointer ( qccoordinates3, kqm ) ;

                    /* . Initialize some accumulators. */
                    dgx = dgy = dgz = dhx = dhy = dhz = 0.0e+00 ;

                    /*----------------------------------------------------------
                    ! . Loops over shells.
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
                            ddim2   = ( iammax + 1 ) * ( jammax + 1 ) ;
                            jdimm   = ( iammax + 2 ) * ( jammax + 2 ) ;
                            jcbfind = jbasis->shells[jshell].type->cbfindex ;
                            ncfuncj = jbasis->shells[jshell].type->ncbf     ;

                            /* . Set the diagonal block flag. */
                            QDIAGONAL = QIQMJQM && ( ishell == jshell ) ;

                            /* . Get the number of roots. */
                            nroots = ( iammax + jammax + 2 ) / 2 + 1 ;

                            /* . Initialize the integral blocks. */
                            for ( i = 0 ; i < ( ncfunci * ncfuncj ) ; i++ ) gx[i] = gy[i] = gz[i] = hx[i] = hy[i] = hz[i] = 0.0e+00 ;

                            /* . Select the expansion center for the recurrence relations. */
                            if ( iammax >= jammax )
                            {
                               iammaxt = iammax ;
                               jammaxt = jammax ;
                               ddim1   = jammax + 1 ;
                               jdim    = jammax + 2 ;
                               dxijt   = xij ;
                               dyijt   = yij ;
                               dzijt   = zij ;
                               rc      = ri ;
                            }
                            else
                            {
                               iammaxt = jammax ;
                               jammaxt = iammax ;
                               ddim1   = iammax + 1 ;
                               jdim    = iammax + 2 ;
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

                                    /* . Set the primitive exponents. */
                                    if ( iammax >= jammax )
                                    {
                                       ag = ai ;
                                       ah = aj ;
                                    }
                                    else
                                    {
                                       ag = aj ;
                                       ah = ai ;
                                    }

                                    /* . Calculate some factors. */
                                    ab    = aa * expn ;
                                    aandb = aa + expn ;
                                    rho   = ab / aandb ;
                                    dnuc  = expfac * ( facn * qn ) / ( expn * sqrt ( aandb ) ) ;

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
                                        Subsidiary_Integral_Nuclear3C ( iammaxt+1, jammaxt+1, 0, False, False, True, True, b00, b10, bp01, dxijt, dyijt, dzijt, f00,
                                                                        xc00, xcp00, yc00, ycp00, zc00, zcp00, 1, jdim, &xint[m*jdimm], &yint[m*jdimm], &zint[m*jdimm] ) ;
                                        Subsidiary_Integral_Derivative3 ( &xint[m*jdimm], &yint[m*jdimm], &zint[m*jdimm],
                                                                          &xidg[m*ddim2], &yidg[m*ddim2], &zidg[m*ddim2],
                                                                          &xidh[m*ddim2], &yidh[m*ddim2], &zidh[m*ddim2],
                                                                          ag, ah, iammaxt, jammaxt, 0, 1, jdim, 1, ddim1 ) ;
                                    }

                                    /* . Assemble the integrals. */
                                    if ( iammax >= jammax )
                                    {
                                        for ( i = 0, n = 0 ; i < ncfunci ; i++ )
                                        {
   	                                    ix  = CBFPOWX[i+icbfind] * jdim  ;
	                                    iy  = CBFPOWY[i+icbfind] * jdim  ;
	                                    iz  = CBFPOWZ[i+icbfind] * jdim  ;
   	                                    ixd = CBFPOWX[i+icbfind] * ddim1 ;
	                                    iyd = CBFPOWY[i+icbfind] * ddim1 ;
	                                    izd = CBFPOWZ[i+icbfind] * ddim1 ;
                                            ti  = dnuc * ibasis->shells[ishell].primitives[ip].ccbf[i] ;
                                            for ( j = 0 ; j < ncfuncj ; j++, n++ )
                                            {
	                                        jxix  = CBFPOWX[j+jcbfind] + ix  ;
	                                        jyiy  = CBFPOWY[j+jcbfind] + iy  ;
	                                        jziz  = CBFPOWZ[j+jcbfind] + iz  ;
	                                        jxixd = CBFPOWX[j+jcbfind] + ixd ;
	                                        jyiyd = CBFPOWY[j+jcbfind] + iyd ;
	                                        jzizd = CBFPOWZ[j+jcbfind] + izd ;
                                                for ( m = 0, facgx = facgy = facgz = fachx = fachy = fachz = 0.0e+00 ; m < nroots ; m++ )
                                                {
                                                    facgx += xidg[jxixd+m*ddim2] * yint[jyiy+m*jdimm] * zint[jziz+m*jdimm] ;
                                                    facgy += xint[jxix+m*jdimm] * yidg[jyiyd+m*ddim2] * zint[jziz+m*jdimm] ;
                                                    facgz += xint[jxix+m*jdimm] * yint[jyiy+m*jdimm] * zidg[jzizd+m*ddim2] ;
                                                    fachx += xidh[jxixd+m*ddim2] * yint[jyiy+m*jdimm] * zint[jziz+m*jdimm] ;
                                                    fachy += xint[jxix+m*jdimm] * yidh[jyiyd+m*ddim2] * zint[jziz+m*jdimm] ;
                                                    fachz += xint[jxix+m*jdimm] * yint[jyiy+m*jdimm] * zidh[jzizd+m*ddim2] ;
/*printf ( "\na %6d%6d%6d%6d%6d%25.15f%25.15f%25.15f%25.15f", ishell, jshell, i, j, m,
          xidg[jxixd+m*ddim2], xidh[jxixd+m*ddim2], yint[jyiy+m*jdimm], zint[jziz+m*jdimm] ) ;*/
                                                }
                                                tij = ti * jbasis->shells[jshell].primitives[jp].ccbf[j] ;
                                                gx[n] += tij * facgx ;
                                                gy[n] += tij * facgy ;
                                                gz[n] += tij * facgz ;
                                                hx[n] += tij * fachx ;
                                                hy[n] += tij * fachy ;
                                                hz[n] += tij * fachz ;
/*printf ( "\na %6d%6d%6d%6d%25.15f%25.15f%25.15f%25.15f%25.15f", ishell, jshell, i, j, gx[n], hx[n], tij, facgx, fachx ) ;*/
                                            }
                                        }
                                    }
                                    else
                                    {
                                        for ( j = 0 ; j < ncfuncj ; j++ )
                                        {
   	                                    ix  = CBFPOWX[j+jcbfind] * jdim  ;
	                                    iy  = CBFPOWY[j+jcbfind] * jdim  ;
	                                    iz  = CBFPOWZ[j+jcbfind] * jdim  ;
   	                                    ixd = CBFPOWX[j+jcbfind] * ddim1 ;
	                                    iyd = CBFPOWY[j+jcbfind] * ddim1 ;
	                                    izd = CBFPOWZ[j+jcbfind] * ddim1 ;
                                            ti  = dnuc * jbasis->shells[jshell].primitives[jp].ccbf[j] ;
                                            for ( i = 0, n = j ; i < ncfunci ; i++, n+= ncfuncj )
                                            {
	                                        jxix  = CBFPOWX[i+icbfind] + ix  ;
	                                        jyiy  = CBFPOWY[i+icbfind] + iy  ;
	                                        jziz  = CBFPOWZ[i+icbfind] + iz  ;
	                                        jxixd = CBFPOWX[i+icbfind] + ixd ;
	                                        jyiyd = CBFPOWY[i+icbfind] + iyd ;
	                                        jzizd = CBFPOWZ[i+icbfind] + izd ;
                                                for ( m = 0, facgx = facgy = facgz = fachx = fachy = fachz = 0.0e+00 ; m < nroots ; m++ )
                                                {
                                                    facgx += xidh[jxixd+m*ddim2] * yint[jyiy+m*jdimm] * zint[jziz+m*jdimm] ;
                                                    facgy += xint[jxix+m*jdimm] * yidh[jyiyd+m*ddim2] * zint[jziz+m*jdimm] ;
                                                    facgz += xint[jxix+m*jdimm] * yint[jyiy+m*jdimm] * zidh[jzizd+m*ddim2] ;
                                                    fachx += xidg[jxixd+m*ddim2] * yint[jyiy+m*jdimm] * zint[jziz+m*jdimm] ;
                                                    fachy += xint[jxix+m*jdimm] * yidg[jyiyd+m*ddim2] * zint[jziz+m*jdimm] ;
                                                    fachz += xint[jxix+m*jdimm] * yint[jyiy+m*jdimm] * zidg[jzizd+m*ddim2] ;
/*printf ( "\nb %6d%6d%6d%6d%6d%25.15f%25.15f%25.15f%25.15f", ishell, jshell, i, j, m,
          xidg[jxixd+m*ddim2], xidh[jxixd+m*ddim2], yint[jyiy+m*jdimm], zint[jziz+m*jdimm] ) ;*/
                                                }
                                                tij = ti * ibasis->shells[ishell].primitives[ip].ccbf[i] ;
                                                gx[n] += tij * facgx ;
                                                gy[n] += tij * facgy ;
                                                gz[n] += tij * facgz ;
                                                hx[n] += tij * fachx ;
                                                hy[n] += tij * fachy ;
                                                hz[n] += tij * fachz ;
/*printf ( "\nb %6d%6d%6d%6d%25.15f%25.15f%25.15f%25.15f%25.15f", ishell, jshell, i, j, gx[n], hx[n], tij, facgx, fachx ) ;*/
                                            }
                                        }
                                    }
                                }
                            }

                            /* . Transform the integrals. */
                            /* if ( qcAtoms->QTOSPHERICAL )
                            {
	                        Integral_Block_Transform_M ( gx, iangmom, ncfunci, nfunci, jangmom, ncfuncj, nfuncj ) ;
	                        Integral_Block_Transform_M ( gy, iangmom, ncfunci, nfunci, jangmom, ncfuncj, nfuncj ) ;
	                        Integral_Block_Transform_M ( gz, iangmom, ncfunci, nfunci, jangmom, ncfuncj, nfuncj ) ;
	                        Integral_Block_Transform_M ( hx, iangmom, ncfunci, nfunci, jangmom, ncfuncj, nfuncj ) ;
	                        Integral_Block_Transform_M ( hy, iangmom, ncfunci, nfunci, jangmom, ncfuncj, nfuncj ) ;
	                        Integral_Block_Transform_M ( hz, iangmom, ncfunci, nfunci, jangmom, ncfuncj, nfuncj ) ;
                            }
                            */

                            /* . Get the scale factor. */
                            if ( QDIAGONAL ) scale = 1.0e+00 ;
                            else             scale = 2.0e+00 ;

                            /* .  Add in the blocks of integrals to the derivatives. */
                            /* . i is usually greater than j. */
                            for ( i = 0, n = 0 ; i < ibasis->shells[ishell].nbasisw ; i++ )
                            {
                                ii = istart + ibasis->shells[ishell].nstartw + i ;
                                for ( j = 0 ; j < jbasis->shells[jshell].nbasisw ; j++, n++ )
                                {
                                    jj = jstart + jbasis->shells[jshell].nstartw + j ;
                                    if ( ii >= jj ) ij = ( ii * ( ii + 1 ) ) / 2 + jj ;
                                    else            ij = ( jj * ( jj + 1 ) ) / 2 + ii ;
                                    fac  = scale * odensityt->data[ij] ;
                                    dgx += fac * gx[n] ;
                                    dgy += fac * gy[n] ;
                                    dgz += fac * gz[n] ;
                                    dhx += fac * hx[n] ;
                                    dhy += fac * hy[n] ;
                                    dhz += fac * hz[n] ;
/*printf ( "\n%6d%6d%6d%6d%25.15f%25.15f%25.15f%25.15f%25.15f", ishell, jshell, i, j, dgx, dhx, fac, gx[n], hx[n] ) ;*/
	                        }
                            }
                        }
                    }

                    /* . Sum in the contribution to the gradient. */
                    Coordinates3_IncrementRow ( qcgradients3, iqm, dgx, dgy, dgz ) ;
                    Coordinates3_IncrementRow ( qcgradients3, jqm, dhx, dhy, dhz ) ;
                    Coordinates3_DecrementRow ( qcgradients3, kqm, dgx, dgy, dgz ) ;
                    Coordinates3_DecrementRow ( qcgradients3, kqm, dhx, dhy, dhz ) ;
                }
            }
        }
    }
}

/*------------------------------------------------------------------------------
! . Calculate the fit-fit derivatives.
!-----------------------------------------------------------------------------*/
void GaussianBasis_Fit_FitD ( const QCAtomContainer *qcAtoms, const QCParameter *qcParameters, Coordinates3 *qccoordinates3, const Real1DArray *fpotential, const Real1DArray *wvector, Coordinates3 *qcgradients3 )
{
    if ( ( qcAtoms != NULL ) && ( qcParameters != NULL ) && ( qccoordinates3 != NULL ) && ( fpotential != NULL ) && ( qcgradients3 != NULL ) )
    {
        auto Real                aandb, ab, ai, aj, bp01, b00, b10, c1x, c1y, c1z, c3x, c3y, c3z, dfi, dfij, dfx, dfy, dfz, fac, facx, facy, facz, fac2, fi, fj, f00,
                                   rho, rij2, ti, tij, u2, wi, wj, xc00, xcp00, xij, yc00, ycp00, yij, zc00, zcp00, zij ;
        auto Real                gx[MAXCBF*MAXCBF], gy[MAXCBF*MAXCBF], gz[MAXCBF*MAXCBF],
                                   xind[MAXAMP1*MAXAMP1*MAXRYS], yind[MAXAMP1*MAXAMP1*MAXRYS], zind[MAXAMP1*MAXAMP1*MAXRYS],
                                   xint[MAXAMP1*MAXAMP2*MAXRYS], yint[MAXAMP1*MAXAMP2*MAXRYS], zint[MAXAMP1*MAXAMP2*MAXRYS] ;
        auto Real               *ri, *rj ;
        auto Integer                   i, iammax, ic, icbfind, iqm, ip, ishell, istart, ix, iy, iz,
                                   j, jammax, jc, jcbfind, jdim, jdimd,  jdimm, jp, jqm, jshell, jstart, jxix, jyiy, jziz, m, n, ncfunci, ncfuncj, nroots ;
        auto GaussianBasis   *ibasis, *jbasis ;
        auto RysQuadrature    roots ;

        /* . Outer loop over centers. */
        for ( iqm = 0 ; iqm < qcAtoms->natoms ; iqm++ )
        {

            /* . Get data for the center. */
            ic     = qcAtoms->data[iqm].center ;
            ibasis = qcParameters->centers[ic].densitybasis ;
            istart = qcAtoms->data[iqm].dstartw ;
            ri     = Coordinates3_RowPointer ( qccoordinates3, iqm ) ;

            /* . Inner loop over centers - no diagonal terms. */
            for ( jqm = 0 ; jqm <= iqm ; jqm++ )
            {

                /* . Get data for the center. */
                jc     = qcAtoms->data[jqm].center ;
                jbasis = qcParameters->centers[jc].densitybasis ;
                jstart = qcAtoms->data[jqm].dstartw ;
                rj     = Coordinates3_RowPointer ( qccoordinates3, jqm ) ;

                /* . Calculate some distance factors. */
	        xij  = ri[0] - rj[0] ;
	        yij  = ri[1] - rj[1] ;
	        zij  = ri[2] - rj[2] ;
	        rij2 = xij * xij + yij * yij + zij * zij ;

                /* . Initialize the accumulators. */
                dfx = dfy = dfz = 0.0e+00 ;

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
                        jdimd   = ( iammax + 1 ) * ( jammax + 1 ) ;
                        jdimm   = ( iammax + 2 ) * ( jammax + 1 ) ;
                        jcbfind = jbasis->shells[jshell].type->cbfindex ;
                        ncfuncj = jbasis->shells[jshell].type->ncbf     ;

                        /* . Get the number of roots. */
                        nroots = ( iammax + jammax + 1 ) / 2 + 1 ;

                        /* . Initialize the integral blocks. */
                        for ( i = 0 ; i < ( ncfunci * ncfuncj ) ; i++ ) gx[i] = gy[i] = gz[i] = 0.0e+00 ;

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
                                    Subsidiary_Integral_Nuclear2C ( iammax+1, jammax, b00, b10, bp01, f00, xc00, xcp00, yc00, ycp00, zc00, zcp00,
                                                                    jdim, &xint[m*jdimm], &yint[m*jdimm], &zint[m*jdimm] ) ;
                                    Subsidiary_Integral_Derivative2 ( &xint[m*jdimm], &yint[m*jdimm], &zint[m*jdimm], ai, iammax, jammax,
                                                                      jdim, &xind[m*jdimd], &yind[m*jdimd], &zind[m*jdimd] ) ;
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
                                        for ( m = 0, facx = facy = facz = 0.0e+00 ; m < nroots ; m++ )
                                        {
                                            facx += xind[jxix+m*jdimd] * yint[jyiy+m*jdimm] * zint[jziz+m*jdimm] ;
                                            facy += xint[jxix+m*jdimm] * yind[jyiy+m*jdimd] * zint[jziz+m*jdimm] ;
                                            facz += xint[jxix+m*jdimm] * yint[jyiy+m*jdimm] * zind[jziz+m*jdimd] ;
                                        }
                                        tij = ti * jbasis->shells[jshell].primitives[jp].ccbf[j] ;
                                        gx[n] += tij * facx ;
                                        gy[n] += tij * facy ;
                                        gz[n] += tij * facz ;
                                        n++ ;
                                    }
                                }
                            }
                        }

                        /* . Transform the integrals. */
                        /* if ( qcAtoms->QTOSPHERICAL )
                        {
	                    Integral_Block_Transform_M ( gx, iangmom, ncfunci, nfunci, jangmom, ncfuncj, nfuncj ) ;
	                    Integral_Block_Transform_M ( gy, iangmom, ncfunci, nfunci, jangmom, ncfuncj, nfuncj ) ;
	                    Integral_Block_Transform_M ( gz, iangmom, ncfunci, nfunci, jangmom, ncfuncj, nfuncj ) ;
                        }
                        */

                        /* .  Add in the blocks of integrals to the derivatives. */
                        for ( i = 0, n = 0 ; i < ibasis->shells[ishell].nbasisw ; i++ )
                        {
                            fi = fpotential->data[istart + ibasis->shells[ishell].nstartw+i] ;
                            wi = wvector->data[istart + ibasis->shells[ishell].nstartw+i] ;
                            for ( j = 0 ; j < jbasis->shells[jshell].nbasisw ; j++, n++ )
                            {
                                fj   = fpotential->data[jstart + jbasis->shells[jshell].nstartw+j] ;
                                wj   = wvector->data[jstart + jbasis->shells[jshell].nstartw+j]    ;
                                fac  = - ( ( fi + wi ) * ( fj + wj ) - wi * wj ) ;
                                dfx += fac * gx[n] ;
                                dfy += fac * gy[n] ;
                                dfz += fac * gz[n] ;
	                    }
                        }
                    }
                }

                /* . Sum in the contribution to the gradient. */
                Coordinates3_IncrementRow ( qcgradients3, iqm, dfx, dfy, dfz ) ;
                Coordinates3_DecrementRow ( qcgradients3, jqm, dfx, dfy, dfz ) ;
            }
        }
    }
}

/*------------------------------------------------------------------------------
! . Calculate the kinetic energy and overlap integrals.
!-----------------------------------------------------------------------------*/
void GaussianBasis_Kinetic_2OverlapD ( const QCAtomContainer *qcAtoms, const QCParameter *qcParameters, Coordinates3 *qccoordinates3,
                                       const SymmetricMatrix *odensityt, const SymmetricMatrix *odensityw, Coordinates3 *qcgradients3 )
{
    if ( ( qcAtoms != NULL ) && ( qcParameters != NULL ) && ( qccoordinates3 != NULL ) && ( odensityt != NULL ) && ( odensityw != NULL ) && ( qcgradients3 != NULL ) )
    {
        auto Real              aa, aainv, ai, aj, arri, dfx, dfy, dfz, expfac, fac, facs, fact, rij2, ti, tij, xij, yij, zij ;
        auto Real              ar[3], ari[3], sx[MAXCBF*MAXCBF], sy[MAXCBF*MAXCBF], sz[MAXCBF*MAXCBF],
                                                tx[MAXCBF*MAXCBF], ty[MAXCBF*MAXCBF], tz[MAXCBF*MAXCBF],
                                                xo[MAXAMP2*MAXAMP3], yo[MAXAMP2*MAXAMP3], zo[MAXAMP2*MAXAMP3],
                                                xt[MAXAMP1*MAXAMP2], yt[MAXAMP1*MAXAMP2], zt[MAXAMP1*MAXAMP2],
                                                xod[MAXAMP1*MAXAMP3], yod[MAXAMP1*MAXAMP3], zod[MAXAMP1*MAXAMP3],
                                                xtd[MAXAMP1*MAXAMP1], ytd[MAXAMP1*MAXAMP1], ztd[MAXAMP1*MAXAMP1] ;
        auto Real             *ri, *rj ;
        auto Integer                 i, iammax, ic, icbfind, ij, iqm, ip, ishell, istart, ixo, ixt, iyo, iyt, izo, izt,
                                 j, jammax, jc, jcbfind, jdimo, jdimt, jp, jqm, jshell, jstart, jxixo, jyiyo, jzizo, jxixt, jyiyt, jzizt, n, ncfunci, ncfuncj ;
        auto GaussianBasis *ibasis, *jbasis ;

        /* . Outer loop over centers. */
        for ( iqm = 0 ; iqm < qcAtoms->natoms ; iqm++ )
        {

            /* . Get data for the center. */
            ic     = qcAtoms->data[iqm].center ;
            ibasis = qcParameters->centers[ic].orbitalbasis ;
            istart = qcAtoms->data[iqm].ostartw ;
            ri     = Coordinates3_RowPointer ( qccoordinates3, iqm ) ;

            /* . Inner loop over centers - there are no diagonal terms. */
            for ( jqm = 0 ; jqm < iqm ; jqm++ )
            {

                /* . Get data for the center. */
                jc     = qcAtoms->data[jqm].center ;
                jbasis = qcParameters->centers[jc].orbitalbasis ;
                jstart = qcAtoms->data[jqm].ostartw ;
                rj     = Coordinates3_RowPointer ( qccoordinates3, jqm ) ;

                /* . Calculate some distance factors. */
	        xij  = ri[0] - rj[0] ;
	        yij  = ri[1] - rj[1] ;
	        zij  = ri[2] - rj[2] ;
	        rij2 = xij * xij + yij * yij + zij * zij ;

                /* . Initialize the accumulators. */
                dfx = dfy = dfz = 0.0e+00 ;

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
                        jdimo   = jammax + 3 ;
                        jdimt   = jammax + 1 ;
                        jcbfind = jbasis->shells[jshell].type->cbfindex ;
                        ncfuncj = jbasis->shells[jshell].type->ncbf     ;

                        /* . Initialize the integral blocks. */
                        for ( i = 0 ; i < ( ncfunci * ncfuncj ) ; i++ )
                        {
                            sx[i] = sy[i] = sz[i] = 0.0e+00 ;
                            tx[i] = ty[i] = tz[i] = 0.0e+00 ;
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
                                Subsidiary_Integral_Overlap2 ( xo, yo, zo, aa, ar, ri, rj, iammax + 1, jammax + 2 ) ;
                                Subsidiary_Integral_Kinetic  ( xo, yo, zo, xt, yt, zt, aj, iammax + 1, jammax, jdimo, jdimt ) ;
                                Subsidiary_Integral_Derivative2 ( xo, yo, zo, ai, iammax, jammax, jdimo, xod, yod, zod ) ;
                                Subsidiary_Integral_Derivative2 ( xt, yt, zt, ai, iammax, jammax, jdimt, xtd, ytd, ztd ) ;

                                /* . Add in the contributions to the full integrals. */
                                for ( i = 0, n = 0 ; i < ncfunci ; i++ )
                                {
   	                            ixo = CBFPOWX[i+icbfind] * jdimo ;
	                            iyo = CBFPOWY[i+icbfind] * jdimo ;
	                            izo = CBFPOWZ[i+icbfind] * jdimo ;
   	                            ixt = CBFPOWX[i+icbfind] * jdimt ;
	                            iyt = CBFPOWY[i+icbfind] * jdimt ;
	                            izt = CBFPOWZ[i+icbfind] * jdimt ;
                                    ti  = expfac * ibasis->shells[ishell].primitives[ip].ccbf[i] ;
                                    for ( j = 0 ; j < ncfuncj ; j++, n++ )
                                    {
	                                jxixo = CBFPOWX[j+jcbfind] + ixo ;
	                                jyiyo = CBFPOWY[j+jcbfind] + iyo ;
	                                jzizo = CBFPOWZ[j+jcbfind] + izo ;
	                                jxixt = CBFPOWX[j+jcbfind] + ixt ;
	                                jyiyt = CBFPOWY[j+jcbfind] + iyt ;
	                                jzizt = CBFPOWZ[j+jcbfind] + izt ;
                                        tij   = ti * jbasis->shells[jshell].primitives[jp].ccbf[j] ;
    		                        sx[n] += tij * xod[jxixo] * yo[jyiyo] * zo[jzizo] ;
    		                        sy[n] += tij * xo[jxixo] * yod[jyiyo] * zo[jzizo] ;
    		                        sz[n] += tij * xo[jxixo] * yo[jyiyo] * zod[jzizo] ;
                                        tx[n] += tij * ( xtd[jxixt] * yo[jyiyo] * zo[jzizo] +
                                                         xod[jxixo] * yt[jyiyt] * zo[jzizo] +
                                                         xod[jxixo] * yo[jyiyo] * zt[jzizt] ) ;
                                        ty[n] += tij * ( xt[jxixt] * yod[jyiyo] * zo[jzizo] +
                                                         xo[jxixo] * ytd[jyiyt] * zo[jzizo] +
                                                         xo[jxixo] * yod[jyiyo] * zt[jzizt] ) ;
                                        tz[n] += tij * ( xt[jxixt] * yo[jyiyo] * zod[jzizo] +
                                                         xo[jxixo] * yt[jyiyt] * zod[jzizo] +
                                                         xo[jxixo] * yo[jyiyo] * ztd[jzizt] ) ;
/*printf ( "\n%6d%6d%6d%6d%25.15f%25.15f%25.15f%25.15f%25.15f%25.15f", ishell, jshell, i, j, xod[jxixt], yo[jyiyo], zo[jzizo], xtd[jxixt], yt[jyiyt], zt[jzizt] ) ;*/
                                    }
                                }
                            }
                        }

                        /* . Transform the integrals. */
                        /* if ( qcAtoms->QTOSPHERICAL )
                        {
	                    Integral_Block_Transform_M ( sx, iangmom, ncfunci, nfunci, jangmom, ncfuncj, nfuncj ) ;
	                    Integral_Block_Transform_M ( sy, iangmom, ncfunci, nfunci, jangmom, ncfuncj, nfuncj ) ;
	                    Integral_Block_Transform_M ( sz, iangmom, ncfunci, nfunci, jangmom, ncfuncj, nfuncj ) ;
	                    Integral_Block_Transform_M ( tx, iangmom, ncfunci, nfunci, jangmom, ncfuncj, nfuncj ) ;
	                    Integral_Block_Transform_M ( ty, iangmom, ncfunci, nfunci, jangmom, ncfuncj, nfuncj ) ;
	                    Integral_Block_Transform_M ( tz, iangmom, ncfunci, nfunci, jangmom, ncfuncj, nfuncj ) ;
                        } */

                        /* .  Add in the blocks of integrals to the derivatives. */
                        for ( i = 0, n = 0 ; i < ibasis->shells[ishell].nbasisw ; i++ )
                        {
                            ij = ( ( istart + ibasis->shells[ishell].nstartw + i ) * ( ( istart + ibasis->shells[ishell].nstartw + i ) + 1 ) ) / 2 + jstart + jbasis->shells[jshell].nstartw ;
                            for ( j = 0 ; j < jbasis->shells[jshell].nbasisw ; ij++, j++, n++ )
                            {
                                facs  =           odensityw->data[ij] ;
                                fact  = 2.0e+00 * odensityt->data[ij] ;
                                dfx += facs * sx[n] + fact * tx[n] ;
                                dfy += facs * sy[n] + fact * ty[n] ;
                                dfz += facs * sz[n] + fact * tz[n] ;
	                    }
                        }
                    }
                }

                /* . Sum in the contribution to the gradient. */
                Coordinates3_IncrementRow ( qcgradients3, iqm, dfx, dfy, dfz ) ;
                Coordinates3_DecrementRow ( qcgradients3, jqm, dfx, dfy, dfz ) ;
            }
        }
    }
}

/*------------------------------------------------------------------------------
! . Calculate the nuclear-nuclear terms.
!-----------------------------------------------------------------------------*/
void GaussianBasis_Nuclear_NuclearD ( const QCAtomContainer *qcAtoms, Coordinates3 *qccoordinates3, Coordinates3 *qcgradients3 )
{
    if ( ( qcAtoms != NULL ) && ( qccoordinates3 != NULL ) && ( qcgradients3 != NULL ) )
    {
        auto Integer                 i, iqm, jqm ;
        auto Real              df, dr[3], iandj, ij, rho, rij2, temp, u2, zi, zj ;
        auto Real             *gi, *gj, *xi, *xj ;
        auto QCAtom        *idata, *jdata ;
        auto RysQuadrature  roots ;

        /* . Outer loop over centers. */
        for ( iqm = 0 ; iqm < qcAtoms->natoms ; iqm++ )
        {

            /* . Get data for the center. */
            idata = &qcAtoms->data[iqm] ;
            xi = Coordinates3_RowPointer ( qccoordinates3, iqm ) ;
            gi = Coordinates3_RowPointer ( qcgradients3,   iqm ) ;
            zi = ( Real ) idata->atomicNumber ;

            /* . Inner loop over centers. */
            for ( jqm = 0 ; jqm < iqm ; jqm++ )
            {

                /* . Get data for the center. */
                jdata = &qcAtoms->data[jqm] ;
                xj = Coordinates3_RowPointer ( qccoordinates3, jqm ) ;
                gj = Coordinates3_RowPointer ( qcgradients3,   jqm ) ;
                zj = ( Real ) jdata->atomicNumber ;

                /* . Calculate some distance factors. */
	        dr[0] = xi[0] - xj[0] ;
	        dr[1] = xi[1] - xj[1] ;
	        dr[2] = xi[2] - xj[2] ;
	        rij2  = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2] ;

                /* . Nuclear repulsion. */
                iandj  = idata->widthe + jdata->widthe ;
                ij     = idata->widthe * jdata->widthe ;
                rho    = ij / iandj ;
                RysQuadrature_Roots ( &roots, 1, ( rho * rij2 ) ) ;
                u2     = rho * roots.roots[0] ;
                df     = - 2.0e+00 * PI252 * idata->widthn * jdata->widthn * zi * zj * roots.weights[0] * u2 / ( ( ij + u2 * iandj ) * sqrt ( iandj ) ) ;
                for ( i = 0 ; i < 3 ; i++ )
                {
                    temp   = df * dr[i] ;
                    gi[i] += temp ;
                    gj[i] -= temp ;
                }
            }
        }
    }
}
