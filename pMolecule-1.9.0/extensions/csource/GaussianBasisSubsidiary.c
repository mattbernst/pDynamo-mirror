/*------------------------------------------------------------------------------
! . File      : GaussianBasisSubsidiary.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==============================================================================
! . Subsidiary integral procedures.
!=============================================================================*/

# include <math.h>
# include <stdarg.h>
# include <stdlib.h>

# include "GaussianBasis.h"
# include "GaussianBasisSubsidiary.h"

/* # define CHECKGHPOINTS */
# ifdef CHECKGHPOINTS
# include <stdio.h>
# endif

/*------------------------------------------------------------------------------
! . Determine derivative integrals from input integrals.
! . jdim is the jth-dimension of both x, y, z and xd, yd, zd.
!-----------------------------------------------------------------------------*/
void Subsidiary_Integral_Derivative2 ( const Real *x, const Real *y, const Real *z, const Real a, const Integer ni, const Integer nj,
                                                                                 const Integer jdim, Real *xd, Real *yd, Real *zd )
{
   Integer i, j, n ;
   Real  a2 ;
   a2 = a + a ;
   /* . Loop over the values of nj. */
   for ( j = 0 ; j <= nj ; j++ )
   {
      /* . i = 0. */
      xd[j] = a2 * x[j+jdim] ;
      yd[j] = a2 * y[j+jdim] ;
      zd[j] = a2 * z[j+jdim] ;
      /* . Loop over the remaining values of ni. */
      for ( i = 1 ; i <= ni ; i++ )
      {
         n = j + ( i - 1 ) * jdim ;
	 xd[n+jdim] = a2 * x[n+jdim+jdim] - ( Real ) i * x[n] ;
	 yd[n+jdim] = a2 * y[n+jdim+jdim] - ( Real ) i * y[n] ;
	 zd[n+jdim] = a2 * z[n+jdim+jdim] - ( Real ) i * z[n] ;
      }
   }
}

/*------------------------------------------------------------------------------
! . Determine derivative integrals from input integrals.
!-----------------------------------------------------------------------------*/
void Subsidiary_Integral_Derivative3 ( const Real *x, const Real *y, const Real *z, Real *xg, Real *yg, Real *zg,
                                                        Real *xh, Real *yh, Real *zh, const Real ag, const Real ah,
                                                                                    const Integer ni, const Integer nj, const Integer nf,
                                                                                              const Integer dim1, const Integer dim2,
                                                                                           const Integer ddim1, const Integer ddim2 )
{
    Integer    f, i, j ;
    Real ag2, ah2 ;

    /* . Initialization. */
    ag2 = ag + ag ;
    ah2 = ah + ah ;

    /* . Loop over the values of nf. */
    for ( f = 0 ; f <= nf ; f++ )
    {
        /* . Loop over the values of nj. */
        for ( j = 0 ; j <= nj ; j++ )
        {
            /* . i = 0. */
            xg[f+j*ddim1] = ag2 * x[f+j*dim1+dim2] ;
            yg[f+j*ddim1] = ag2 * y[f+j*dim1+dim2] ;
            zg[f+j*ddim1] = ag2 * z[f+j*dim1+dim2] ;
            /* . Loop over the remaining values of ni. */
            for ( i = 1 ; i <= ni ; i++ )
            {
	        xg[f+j*ddim1+i*ddim2] = ag2 * x[f+j*dim1+(i+1)*dim2] - ( Real ) i * x[f+j*dim1+(i-1)*dim2] ;
	        yg[f+j*ddim1+i*ddim2] = ag2 * y[f+j*dim1+(i+1)*dim2] - ( Real ) i * y[f+j*dim1+(i-1)*dim2] ;
	        zg[f+j*ddim1+i*ddim2] = ag2 * z[f+j*dim1+(i+1)*dim2] - ( Real ) i * z[f+j*dim1+(i-1)*dim2] ;
            }
        }
        /* . Loop over the values of ni. */
        for ( i = 0 ; i <= ni ; i++ )
        {
            /* . j = 0.*/
            xh[f+i*ddim2] = ah2 * x[f+dim1+i*dim2] ;
            yh[f+i*ddim2] = ah2 * y[f+dim1+i*dim2] ;
            zh[f+i*ddim2] = ah2 * z[f+dim1+i*dim2] ;
            /* . Loop over the remaining values of nj. */
            for ( j = 1 ; j <= nj ; j++ )
            {
	        xh[f+j*ddim1+i*ddim2] = ah2 * x[f+(j+1)*dim1+i*dim2] - ( Real ) j * x[f+(j-1)*dim1+i*dim2] ;
	        yh[f+j*ddim1+i*ddim2] = ah2 * y[f+(j+1)*dim1+i*dim2] - ( Real ) j * y[f+(j-1)*dim1+i*dim2] ;
	        zh[f+j*ddim1+i*ddim2] = ah2 * z[f+(j+1)*dim1+i*dim2] - ( Real ) j * z[f+(j-1)*dim1+i*dim2] ;
            }
        }
    }
}

/*------------------------------------------------------------------------------
! . Dipole integrals.
!-----------------------------------------------------------------------------*/
void Subsidiary_Integral_Dipole ( Real *x, Real *y, Real *z, const Real aa, const Real *r0, const Real *ri, const Real *rj,
                                                                                       const Real *center, const Integer ni, const Integer nj )
{
   /* . N-point GH quadrature is exact for a polynomial of order 2n-1.
   ! . Therefore, a pth-order polynomial requires (p+1)/2 points.
   ! . To avoid half-integer values this becomes (p+2)/2. */
   Integer i, j, k, n, npts, p, q ;
   Real  ar, br, ptr, rint[3], t, temp, tinv ;
   t    = sqrt ( aa ) ;
   tinv = 1.0e+00 / t ;
   for ( i = 0, n = 0 ; i <= ni ; i++ )
   {
      for ( j = 0 ; j <= nj ; j++, n++ )
      {
	 for ( k = 0 ; k < 3 ; k++ ) rint[k] = 0.0e+00 ;
	 npts = ( i + j + 2 ) / 2 ;
# ifdef CHECKGHPOINTS
if ( npts > GHMAXPT ) printf ( "\nINVALID NUMBER OF POINTS IN GAUSS-HERMITE QUADRATURE = %d\n", npts ) ;
# endif
	 for ( p = GHFIRST[npts] ; p <= GHLAST[npts] ; p++ )
         {
	    for ( k = 0 ; k < 3 ; k++ )
            {
	       ptr  = GHABSCISSAE[p] / t + r0[k] ;
               ar   = ptr - ri[k] ;
               br   = ptr - rj[k] ;
               temp = GHWEIGHTS[p] * ( ptr - center[k] ) ;
	       for ( q = 1 ; q <= i ; q++ ) temp *= ar ;
	       for ( q = 1 ; q <= j ; q++ ) temp *= br ;
	       rint[k] += temp ;
            }
	 }
	 x[n] = rint[0] * tinv ;
	 y[n] = rint[1] * tinv ;
	 z[n] = rint[2] * tinv ;
      }
   }
}

/*------------------------------------------------------------------------------
! . Determine the kinetic energy subsidiary integrals.
! . The jth dimension of x, y and z is assumed to be nj+3.
!-----------------------------------------------------------------------------*/
void Subsidiary_Integral_Kinetic ( const Real *x, const Real *y, const Real *z, Real *xt, Real *yt, Real *zt, const Real aj,
                                                                            const Integer ni, const Integer nj, const Integer jdimo, const Integer jdimt )
{
    Real a2 ;
    Integer    i, io, it, j ;
    a2 = aj + aj ;
    for ( i = 0, io = 0, it = 0 ; i <= ni ; i++, io += jdimo, it += jdimt )
    {
        xt[it] = ( x[io] - x[io+2] * a2 ) * aj ;
        yt[it] = ( y[io] - y[io+2] * a2 ) * aj ;
        zt[it] = ( z[io] - z[io+2] * a2 ) * aj ;
    }
    if ( nj > 0 )
    {
        for ( i = 0, io = 0, it = 0 ; i <= ni ; i++, io += jdimo, it += jdimt )
        {
            xt[it+1] = ( x[io+1] * 3.0e+00 - x[io+3] * a2 ) * aj ;
            yt[it+1] = ( y[io+1] * 3.0e+00 - y[io+3] * a2 ) * aj ;
            zt[it+1] = ( z[io+1] * 3.0e+00 - z[io+3] * a2 ) * aj ;
        }
        if ( nj > 1 )
        {
            for ( j = 2 ; j <= nj ; j++ )
            {
                for ( i = 0, io = 0, it = 0 ; i <= ni ; i++, io += jdimo, it += jdimt )
                {
	            xt[it+j] = ( x[io+j] * ( Real ) ( 2*j+1 ) - x[io+j+2] * a2 ) * aj - x[io+j-2] * ( Real ) ( (j*(j-1))/2 ) ;
	            yt[it+j] = ( y[io+j] * ( Real ) ( 2*j+1 ) - y[io+j+2] * a2 ) * aj - y[io+j-2] * ( Real ) ( (j*(j-1))/2 ) ;
	            zt[it+j] = ( z[io+j] * ( Real ) ( 2*j+1 ) - z[io+j+2] * a2 ) * aj - z[io+j-2] * ( Real ) ( (j*(j-1))/2 ) ;
                }
            }
        }
    }
}

/*------------------------------------------------------------------------------
! . Determine the 2-center nuclear subsidiary integrals.
! . The first dimension is jdim.
!-----------------------------------------------------------------------------*/
void Subsidiary_Integral_Nuclear2C ( const Integer iangmom, const Integer jangmom, const Real b00, const Real b10, const Real bp01, const Real f00, const Real xc00,
                                                                     const Real xcp00, const Real yc00, const Real ycp00, const Real zc00, const Real zcp00,
                                                                                                            const Integer jdim, Real *xint, Real *yint, Real *zint )
{
    Integer    i, j ;
    Real cp01, cp10 = 0.0e+00, c01, c10 ;

    /* . (0,0). */
    xint[0] = 1.0e+00 ;
    yint[0] = 1.0e+00 ;
    zint[0] = f00 ;

    /* . (1,0). */
    if ( iangmom > 0 )
    {
        xint[jdim] = xc00 ;
        yint[jdim] = yc00 ;
        zint[jdim] = zc00 * f00 ;
    }

    /* . (0,1). */
    if ( jangmom > 0 )
    {
        xint[1] = xcp00 ;
        yint[1] = ycp00 ;
        zint[1] = zcp00 * f00 ;
        /* . (1,1). */
        if ( iangmom > 0 )
        {
            cp10 = b00 ;
            xint[jdim+1] = xcp00 * xint[jdim] + cp10 ;
            yint[jdim+1] = ycp00 * yint[jdim] + cp10 ;
            zint[jdim+1] = zcp00 * zint[jdim] + cp10 * f00 ;
        }
    }

    /* . iangmom > 1. */
    if ( iangmom > 1 )
    {
        c10 = 0.0e+00 ;
        for ( i = 2 ; i <= iangmom ; i++ )
        {
            /* . (i,0). */
            c10 += b10 ;
            xint[i*jdim] = c10 * xint[(i-2)*jdim] + xc00 * xint[(i-1)*jdim] ;
            yint[i*jdim] = c10 * yint[(i-2)*jdim] + yc00 * yint[(i-1)*jdim] ;
            zint[i*jdim] = c10 * zint[(i-2)*jdim] + zc00 * zint[(i-1)*jdim] ;
            /* . (i,1). */
            if ( jangmom > 0 )
            {
                cp10 += b00 ;
                xint[i*jdim+1] = xcp00 * xint[i*jdim] + cp10 * xint[(i-1)*jdim] ;
                yint[i*jdim+1] = ycp00 * yint[i*jdim] + cp10 * yint[(i-1)*jdim] ;
                zint[i*jdim+1] = zcp00 * zint[i*jdim] + cp10 * zint[(i-1)*jdim] ;
            }
        }
    }

    /* . jangmom > 1. */
    if ( jangmom > 1 )
    {
        cp01 = 0.0e+00 ;
        c01  = b00 ;
        for ( j = 2 ; j <= jangmom ; j++ )
        {
            /* . (0,j). */
            cp01 += bp01 ;
            xint[j] = cp01 * xint[j-2] + xcp00 * xint[j-1] ;
            yint[j] = cp01 * yint[j-2] + ycp00 * yint[j-1] ;
            zint[j] = cp01 * zint[j-2] + zcp00 * zint[j-1] ;
            /* . (1,j). */
            if ( iangmom > 0 )
            {
                c01 += b00 ;
                xint[j+jdim] = xc00 * xint[j] + c01 * xint[j-1] ;
                yint[j+jdim] = yc00 * yint[j] + c01 * yint[j-1] ;
                zint[j+jdim] = zc00 * zint[j] + c01 * zint[j-1] ;
            }
        }
    }

    /* . iangmom and jangmom > 1. */
    if ( ( iangmom > 1 ) && ( jangmom > 1 ) )
    {
        /* . (j,i). */
        c01 = b00 ;
        for ( j = 2 ; j <= jangmom ; j++ )
        {
            c01 += b00 ;
            c10  = b10 ;
            for ( i = 2 ; i <= iangmom ; i++ )
            {
               xint[j+i*jdim] = c10 * xint[j+(i-2)*jdim] + xc00 * xint[j+(i-1)*jdim] + c01 * xint[j-1+(i-1)*jdim] ;
               yint[j+i*jdim] = c10 * yint[j+(i-2)*jdim] + yc00 * yint[j+(i-1)*jdim] + c01 * yint[j-1+(i-1)*jdim] ;
               zint[j+i*jdim] = c10 * zint[j+(i-2)*jdim] + zc00 * zint[j+(i-1)*jdim] + c01 * zint[j-1+(i-1)*jdim] ;
               c10 += b10 ;
            }
        }
    }
}

/*------------------------------------------------------------------------------
! . Determine the 3-center nuclear subsidiary integrals.
! . The first dimension is jdim1 and the second dimension jdim2 which means
! .that the element [k,j,i] is accessed as k+j*jdim1+i*jdim2.
!-----------------------------------------------------------------------------*/
void Subsidiary_Integral_Nuclear3C ( const Integer ni, const Integer nj, const Integer nk, const Boolean qij0, const Boolean qij1, const Boolean qn0, const Boolean qn1,
                                     const Real b00, const Real b10, const Real bp01, const Real dxij, const Real dyij, const Real dzij,
                                     const Real f00, const Real xc00, const Real xcp00, const Real  yc00, const Real  ycp00, const Real zc00,
                                     const Real zcp00, const Integer jdim1, const Integer jdim2, Real *xint, Real *yint, Real *zint )
{
    Integer    i, j, k, m ;
    Real cp01, cp10 = 0.0e+00, c01, c10 ;

    /* . i(0,0). */
    xint[0] = 1.0e+00 ;
    yint[0] = 1.0e+00 ;
    zint[0] = f00 ;
    if ( qij0 && qn0 ) return ;
    if ( ! qij0 )
    {
        /* . i(1,0). */
        xint[jdim2] = xc00 ;
        yint[jdim2] = yc00 ;
        zint[jdim2] = zc00 * f00 ;
    }
    if ( ! qn0 )
    {
        /* . i(0,1). */
        xint[1] = xcp00 ;
        yint[1] = ycp00 ;
        zint[1] = zcp00 * f00 ;
        if ( ! qij0 )
        {
            /* . i(1,1).*/
            cp10 = b00 ;
            xint[1+jdim2] = xcp00*xint[jdim2] + cp10 ;
            yint[1+jdim2] = ycp00*yint[jdim2] + cp10 ;
            zint[1+jdim2] = zcp00*zint[jdim2] + cp10 * f00 ;
        }
    }
    if ( ! qij1 )
    {
        c10 = 0.0e+00 ;
        for ( i = 2 ; i <= ni ; i++ )
        {
            /* . i(i,0).*/
            c10 += b10 ;
            xint[i*jdim2] = c10*xint[(i-2)*jdim2]+xc00*xint[(i-1)*jdim2] ;
            yint[i*jdim2] = c10*yint[(i-2)*jdim2]+yc00*yint[(i-1)*jdim2] ;
            zint[i*jdim2] = c10*zint[(i-2)*jdim2]+zc00*zint[(i-1)*jdim2] ;
            /* . i(i,1).*/
            if ( ! qn0 )
            {
                cp10 += b00 ;
                xint[1+i*jdim2] = xcp00*xint[i*jdim2]+cp10*xint[(i-1)*jdim2] ;
                yint[1+i*jdim2] = ycp00*yint[i*jdim2]+cp10*yint[(i-1)*jdim2] ;
                zint[1+i*jdim2] = zcp00*zint[i*jdim2]+cp10*zint[(i-1)*jdim2] ;
            }
        }
        for ( j = 1 ; j <= nj ; j++ )
        {
            /* . i(ni,j,0).*/
            c10 += b10 ;
            if ( j == 1 )
            {
                xint[jdim1+ni*jdim2] = c10*xint[(ni-1)*jdim2]+xc00*xint[ni*jdim2] ;
                yint[jdim1+ni*jdim2] = c10*yint[(ni-1)*jdim2]+yc00*yint[ni*jdim2] ;
                zint[jdim1+ni*jdim2] = c10*zint[(ni-1)*jdim2]+zc00*zint[ni*jdim2] ;
            }
            else
            {
                xint[j*jdim1+ni*jdim2] = c10*xint[(j-2)*jdim1+ni*jdim2]+xc00*xint[(j-1)*jdim1+ni*jdim2] ;
                yint[j*jdim1+ni*jdim2] = c10*yint[(j-2)*jdim1+ni*jdim2]+yc00*yint[(j-1)*jdim1+ni*jdim2] ;
                zint[j*jdim1+ni*jdim2] = c10*zint[(j-2)*jdim1+ni*jdim2]+zc00*zint[(j-1)*jdim1+ni*jdim2] ;
            }
            /* . i(ni,j,1).*/
            if ( ! qn0 )
            {
                cp10 += b00 ;
                xint[1+j*jdim1+ni*jdim2] = xcp00*xint[j*jdim1+ni*jdim2]+cp10*xint[(j-1)*jdim1+ni*jdim2] ;
                yint[1+j*jdim1+ni*jdim2] = ycp00*yint[j*jdim1+ni*jdim2]+cp10*yint[(j-1)*jdim1+ni*jdim2] ;
                zint[1+j*jdim1+ni*jdim2] = zcp00*zint[j*jdim1+ni*jdim2]+cp10*zint[(j-1)*jdim1+ni*jdim2] ;
            }
        }
    }
    if ( ! qn1 )
    {
        cp01 = 0.0e+00 ;
        c01  = b00 ;
        for ( k = 2 ; k <= nk ; k++ )
        {
            /* . i(0,k).*/
            cp01 += bp01 ;
            xint[k] = cp01*xint[k-2]+xcp00*xint[k-1] ;
            yint[k] = cp01*yint[k-2]+ycp00*yint[k-1] ;
            zint[k] = cp01*zint[k-2]+zcp00*zint[k-1] ;
            /* . i(1,k).*/
            if ( ! qij0 )
            {
                c01 += b00 ;
                xint[k+jdim2] = xc00*xint[k]+c01*xint[k-1] ;
                yint[k+jdim2] = yc00*yint[k]+c01*yint[k-1] ;
                zint[k+jdim2] = zc00*zint[k]+c01*zint[k-1] ;
            }
        }
    }
    if ( ! ( qij1 || qn1 ) )
    {
        /* . i(n,m).*/
        c01 = b00 ;
        for ( k = 2 ; k <= nk ; k++ )
        {
            c01 += b00 ;
            c10  = b10 ;
            for ( i = 2 ; i <= ni ; i++ )
            {
               xint[k+i*jdim2] = c10*xint[k+(i-2)*jdim2] + xc00*xint[k+(i-1)*jdim2] + c01*xint[k-1+(i-1)*jdim2] ;
               yint[k+i*jdim2] = c10*yint[k+(i-2)*jdim2] + yc00*yint[k+(i-1)*jdim2] + c01*yint[k-1+(i-1)*jdim2] ;
               zint[k+i*jdim2] = c10*zint[k+(i-2)*jdim2] + zc00*zint[k+(i-1)*jdim2] + c01*zint[k-1+(i-1)*jdim2] ;
               c10 += b10 ;
            }
            for ( j = 1 ; j <= nj ; j++ )
            {
                if ( j == 1 )
                {
                    xint[k+jdim1+ni*jdim2] = c10*xint[k+(ni-1)*jdim2] + xc00*xint[k+ni*jdim2] + c01*xint[k-1+ni*jdim2] ;
                    yint[k+jdim1+ni*jdim2] = c10*yint[k+(ni-1)*jdim2] + yc00*yint[k+ni*jdim2] + c01*yint[k-1+ni*jdim2] ;
                    zint[k+jdim1+ni*jdim2] = c10*zint[k+(ni-1)*jdim2] + zc00*zint[k+ni*jdim2] + c01*zint[k-1+ni*jdim2] ;
                }
                else
                {
                    xint[k+j*jdim1+ni*jdim2] = c10*xint[k+(j-2)*jdim1+ni*jdim2] + xc00*xint[k+(j-1)*jdim1+ni*jdim2] + c01*xint[k-1+(j-1)*jdim1+ni*jdim2] ;
                    yint[k+j*jdim1+ni*jdim2] = c10*yint[k+(j-2)*jdim1+ni*jdim2] + yc00*yint[k+(j-1)*jdim1+ni*jdim2] + c01*yint[k-1+(j-1)*jdim1+ni*jdim2] ;
                    zint[k+j*jdim1+ni*jdim2] = c10*zint[k+(j-2)*jdim1+ni*jdim2] + zc00*zint[k+(j-1)*jdim1+ni*jdim2] + c01*zint[k-1+(j-1)*jdim1+ni*jdim2] ;
                }
                c10 += b10 ;
            }
        }
    }
    if ( nj > 0 )
    {
        /* . i(ni,nj,m).*/
        for ( k = 0 ; k <= nk ; k++ )
        {
            for ( m = 0 ; m <= (nj-1) ; m++ )
            {
                for ( j = nj ; j >= (m+1) ; j-- )
                {
                    xint[k+j*jdim1+ni*jdim2] = xint[k+j*jdim1+ni*jdim2]+dxij*xint[k+(j-1)*jdim1+ni*jdim2] ;
                    yint[k+j*jdim1+ni*jdim2] = yint[k+j*jdim1+ni*jdim2]+dyij*yint[k+(j-1)*jdim1+ni*jdim2] ;
                    zint[k+j*jdim1+ni*jdim2] = zint[k+j*jdim1+ni*jdim2]+dzij*zint[k+(j-1)*jdim1+ni*jdim2] ;
                }
            }
            if ( ni > 0 )
            {
                for ( j = 1 ; j <= nj ; j++ )
                {
                    for ( i = 0 ; i <= (ni-1) ; i++ )
                    {
                        xint[k+j*jdim1+i*jdim2] = xint[k+(j-1)*jdim1+(i+1)*jdim2]+dxij*xint[k+(j-1)*jdim1+i*jdim2] ;
                        yint[k+j*jdim1+i*jdim2] = yint[k+(j-1)*jdim1+(i+1)*jdim2]+dyij*yint[k+(j-1)*jdim1+i*jdim2] ;
                        zint[k+j*jdim1+i*jdim2] = zint[k+(j-1)*jdim1+(i+1)*jdim2]+dzij*zint[k+(j-1)*jdim1+i*jdim2] ;
                    }
                }
            }
        }
    }
}

/*------------------------------------------------------------------------------
! . Double overlap integrals.
!-----------------------------------------------------------------------------*/
void Subsidiary_Integral_Overlap2 ( Real *x, Real *y, Real *z, const Real aa, const Real *r0, const Real *ri, const Real *rj, const Integer ni, const Integer nj )
{
   /* . N-point GH quadrature is exact for a polynomial of order 2n-1.
   ! . Therefore, a pth-order polynomial requires (p+1)/2 points.
   ! . To avoid half-integer values this becomes (p+2)/2. */
   Integer i, j, k, n, npts, p, q ;
   Real  ar, br, ptr, rint[3], t, temp, tinv ;
   t    = sqrt ( aa ) ;
   tinv = 1.0e+00 / t ;
   for ( i = 0, n = 0 ; i <= ni ; i++ )
   {
      for ( j = 0 ; j <= nj ; j++, n++ )
      {
	 for ( k = 0 ; k < 3 ; k++ ) rint[k] = 0.0e+00 ;
	 npts = ( i + j + 2 ) / 2 - 1 ;
# ifdef CHECKGHPOINTS
if ( npts > GHMAXPT ) printf ( "\nINVALID NUMBER OF POINTS IN GAUSS-HERMITE QUADRATURE = %d\n", npts ) ;
# endif
	 for ( p = GHFIRST[npts] ; p <= GHLAST[npts] ; p++ )
         {
	    for ( k = 0 ; k < 3 ; k++ )
            {
	       ptr  = GHABSCISSAE[p] / t + r0[k] ;
               ar   = ptr - ri[k] ;
               br   = ptr - rj[k] ;
               temp = GHWEIGHTS[p] ;
	       for ( q = 1 ; q <= i ; q++ ) temp *= ar ;
	       for ( q = 1 ; q <= j ; q++ ) temp *= br ;
	       rint[k] += temp ;
            }
	 }
	 x[n] = rint[0] * tinv ;
	 y[n] = rint[1] * tinv ;
	 z[n] = rint[2] * tinv ;
      }
   }
}

/*------------------------------------------------------------------------------
! . Triple overlap integrals.
!-----------------------------------------------------------------------------*/
void Subsidiary_Integral_Overlap3 ( Real *x, Real *y, Real *z, const Real aa, const Real *r0, const Real *ri, const Real *rj, const Real *rk, const Integer ni, const Integer nj, const Integer nk )
{
    Real ar, br, cr, ptr, rint[3], t, temp, tinv ;
    Integer    i, j, k, m, n, npts, p, q ;
    t    = sqrt ( aa ) ;
    tinv = 1.0e+00 / t ;
    for ( i = 0, n = 0 ; i <= ni ; i++ )
    {
        for ( j = 0 ; j <= nj ; j++ )
        {
	    for ( k = 0 ; k <= nk ; k++, n++ )
            {
                for ( m = 0 ; m < 3 ; m++ ) rint[m] = 0.0e+00 ;
	        npts = ( i + j + k + 2 ) / 2 ;
# ifdef CHECKGHPOINTS
if ( npts > GHMAXPT ) printf ( "\nINVALID NUMBER OF POINTS IN GAUSS-HERMITE QUADRATURE = %d\n", npts ) ;
# endif
	        for ( p = GHFIRST[npts] ; p <= GHLAST[npts] ; p++ )
                {
	            for ( m = 0 ; m < 3 ; m++ )
                    {
	                ptr  = GHABSCISSAE[p] / t + r0[m] ;
                        ar   = ptr - ri[m] ;
                        br   = ptr - rj[m] ;
                        cr   = ptr - rk[m] ;
                        temp = GHWEIGHTS[p] ;
	                for ( q = 1 ; q <= i ; q++ ) temp *= ar ;
	                for ( q = 1 ; q <= j ; q++ ) temp *= br ;
	                for ( q = 1 ; q <= k ; q++ ) temp *= cr ;
	                rint[m] += temp ;
                    }
	        }
	        x[n] = rint[0] * tinv ;
	        y[n] = rint[1] * tinv ;
	        z[n] = rint[2] * tinv ;
            }
        }
    }
}
