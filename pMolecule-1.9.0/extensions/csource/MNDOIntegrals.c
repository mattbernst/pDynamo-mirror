/*------------------------------------------------------------------------------
! . File      : MNDOIntegrals.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . Procedures for calculating the integrals in a MNDO method.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "BlockStorage.h"
# include "DefineStatements.h"
# include "Definitions.h"
# include "Memory.h"
# include "MNDODefinitions.h"
# include "MNDOIntegrals.h"
# include "MNDOIntegralUtilities.h"
# include "MNDOParameters.h"
# include "Units.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Add in the one-center TEIs.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The maximum number of unique integrals - 1 (s), 16 (sp), 155 (spd). */
void MNDOIntegrals_AddInOneCenterTEIs ( const MNDOParameters *self, const Integer istart, BlockStorage *twoelectronintegrals )
{
    /* . Check that there are integrals. */
    if ( self->nocteis > 0 )
    {
        auto Integer   i ;
        auto Integer16 indices[4*N1CTEIS] ;

        /* . Copy and increment the indices. */
        for ( i = 0 ; i < 4*self->nocteis ; i++ ) indices[i] = self->octeiindices[i] + istart ;

        /* . Save the data. */
        BlockStorage_Data_Add ( twoelectronintegrals, self->nocteis, self->octeivalues, indices, NULL ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Form the two-electron parts of the Fock matrices.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* # define __RANGE_CHECKING */
# define BFINDEX(i) ( i * ( i + 1 ) ) / 2
void MNDOIntegrals_MakeFockG ( BlockStorage *twoelectronintegrals, SymmetricMatrix *densitya, SymmetricMatrix *focka, SymmetricMatrix *densityb, SymmetricMatrix *fockb )
{
   if  ( ( twoelectronintegrals != NULL ) && ( densitya != NULL ) && ( focka != NULL ) )
   {
      auto Block *block ;
      auto Boolean        QAB   ;
      auto Real      value ;
      auto Integer         i, i1, i2, i3, i4, n, nij, nik, nil, njk, njl, nkl, t ;

      /* . Initialization. */
      SymmetricMatrix_Set_Zero ( focka ) ;

      if ( ( densityb == NULL ) || ( fockb == NULL ) )   QAB = False ;
      else { QAB = True  ; SymmetricMatrix_Set_Zero ( fockb ) ; }

      /* . Loop over the integral blocks. */
      List_Iterate_Initialize ( twoelectronintegrals->blocks ) ;
      while ( ( block = BlockStorage_Iterate ( twoelectronintegrals ) ) != NULL )
      {
         /* . Loop over the integrals. */
         for ( i = 0, n = 0 ; i < block->ndata ; i++, n += 4 )
         {
            /* . Get the data. */
            i1    = block->indices16[n  ] ;
            i2    = block->indices16[n+1] ;
            i3    = block->indices16[n+2] ;
            i4    = block->indices16[n+3] ;
            value = block->data[i] ;

            /* . Shuffle the index pairs. */
	    if ( i1 < i2 ) { t = i1 ; i1 = i2 ; i2 = t ; }
            if ( i3 < i4 ) { t = i3 ; i3 = i4 ; i4 = t ; }

            /* . Shuffle the indices of both pairs. */
            if ( ( i1 < i3 ) || ( ( i1 == i3 ) && ( i2 < i4  ) ) ) { t = i1 ; i1 = i3 ; i3 = t ; t = i2 ; i2 = i4 ; i4 = t ; }

            /* . Get the appropriate scaling factors. */
	    if ( i1 == i2 ) value *= 0.5e+00 ;
	    if ( i3 == i4 ) value *= 0.5e+00 ;
            if ( ( i1 == i3 ) && ( i2 == i4 ) ) value *= 0.5e+00 ;

            /* . Get the matrix indices. */
            nij = BFINDEX ( i1 ) + i2 ;
            nkl = BFINDEX ( i3 ) + i4 ;
            nik = BFINDEX ( i1 ) + i3 ;
            nil = BFINDEX ( i1 ) + i4 ;
            if ( i2 > i3 ) njk = BFINDEX ( i2 ) + i3 ;
            else           njk = BFINDEX ( i3 ) + i2 ;
            if ( i2 > i4 ) njl = BFINDEX ( i2 ) + i4 ;
            else           njl = BFINDEX ( i4 ) + i2 ;
# ifdef __RANGE_CHECKING
if ( nij >= densitya->size ) printf ( "QCModelMNDO_FockG: nij out of range: %d\n", nij ) ;
if ( nkl >= densitya->size ) printf ( "QCModelMNDO_FockG: nkl out of range: %d\n", nkl ) ;
if ( nik >= densitya->size ) printf ( "QCModelMNDO_FockG: nik out of range: %d\n", nik ) ;
if ( nil >= densitya->size ) printf ( "QCModelMNDO_FockG: nil out of range: %d\n", nil ) ;
if ( njk >= densitya->size ) printf ( "QCModelMNDO_FockG: njk out of range: %d\n", njk ) ;
if ( njl >= densitya->size ) printf ( "QCModelMNDO_FockG: njl out of range: %d\n", njl ) ;
# endif
            /* . Add in the integrals. */
	    if ( QAB )
            {
	       value *= 2.0e+00 ;
               focka->data[nij] += 2.0e+00 * value * ( densitya->data[nkl] + densityb->data[nkl] ) ;
               focka->data[nkl] += 2.0e+00 * value * ( densitya->data[nij] + densityb->data[nij] ) ;
               fockb->data[nij] += 2.0e+00 * value * ( densitya->data[nkl] + densityb->data[nkl] ) ;
               fockb->data[nkl] += 2.0e+00 * value * ( densitya->data[nij] + densityb->data[nij] ) ;
               focka->data[nik] -=           value *   densitya->data[njl] ;
               focka->data[nil] -=           value *   densitya->data[njk] ;
               focka->data[njk] -=           value *   densitya->data[nil] ;
               focka->data[njl] -=           value *   densitya->data[nik] ;
               fockb->data[nik] -=           value *   densityb->data[njl] ;
               fockb->data[nil] -=           value *   densityb->data[njk] ;
               fockb->data[njk] -=           value *   densityb->data[nil] ;
               fockb->data[njl] -=           value *   densityb->data[nik] ;
            }
            else
            {
               focka->data[nij] += 4.0e+00 * value * densitya->data[nkl] ;
               focka->data[nkl] += 4.0e+00 * value * densitya->data[nij] ;
               focka->data[nik] -=           value * densitya->data[njl] ;
               focka->data[nil] -=           value * densitya->data[njk] ;
               focka->data[njk] -=           value * densitya->data[nil] ;
               focka->data[njl] -=           value * densitya->data[nik] ;
            }
         }
      }
      /* . Scale off-diagonal elements of the matrices by 1/2. */
      SymmetricMatrix_Scale_OD ( focka, 0.5e+00 ) ;
      if ( QAB ) SymmetricMatrix_Scale_OD ( fockb, 0.5e+00 ) ;
   }
}
# undef BFINDEX

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the integrals in the molecular frame.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MNDOIntegrals_MolecularFrame2CIntegrals ( const MNDOParameters *idata, const Integer istart, const Real *xi,
                                               const MNDOParameters *jdata, const Integer jstart, const Real *xj,
                                                         Real *enuc, Real1DArray *mfcore1b, Real1DArray *mfcore2a,
                                                                           BlockStorage *twoelectronintegrals )
{
    if ( ( idata != NULL ) && ( jdata != NULL ) && ( enuc != NULL ) && ( mfcore1b != NULL ) && ( mfcore2a != NULL ) && ( twoelectronintegrals != NULL ) )
    {
        auto Real       r ;
        auto Integer          i, j, k, l, n, ni, nj ;
        auto Integer16        indices[4*N2CTEIS] ;
        auto Real2DArray *hfteis = NULL, *itransformation = NULL, *jtransformation = NULL, *lfteis = NULL, *mfteis = NULL ;
        auto Real1DArray *lfcore1b = NULL, *lfcore2a = NULL ;

        /* . Get the transformation matrices. */
        r = MNDOIntegralUtilities_GetTransformationMatrices ( idata, jdata, xi, xj, &itransformation, &jtransformation, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL ) ;

        /* . Allocate space for the integrals. */
        ni = ( idata->norbitals * ( idata->norbitals + 1 ) ) / 2 ;
        nj = ( jdata->norbitals * ( jdata->norbitals + 1 ) ) / 2 ;
        lfteis = Real2DArray_Allocate ( ni, nj, NULL ) ;
        if ( itransformation == NULL ) lfcore1b = mfcore1b ;
        else                           lfcore1b = Real1DArray_Allocate ( ni, NULL ) ;
        if ( jtransformation == NULL ) lfcore2a = mfcore2a ;
        else                           lfcore2a = Real1DArray_Allocate ( nj, NULL ) ;
        /* . Get the integrals in the local frame. */
        MNDOIntegralUtilities_LocalFrame2CTEIs ( idata, jdata, r, lfteis, lfcore1b, lfcore2a, NULL, NULL, NULL ) ;
# ifdef PRINTMOPACTEIS
{
/* . Printing. */
printf ( "\n\nTWO-ELECTRON INTEGRALS - TRANSFORMATIONS:\n%d %d\n", idata->atomicNumber, jdata->atomicNumber ) ;
printf ( "\nTransformation 1:" ) ;
Matrix_Print ( itransformation ) ;
printf ( "\nTransformation 2:" ) ;
Matrix_Print ( jtransformation ) ;
printf ( "\nLocal Frame TEIs:" ) ;
Matrix_Print ( lfteis ) ;
}
# endif
        /* . Transform from the local to molecular frames. */
        /* . OEIs then TEIs. */
        if ( itransformation == NULL )
        {
            lfcore1b = NULL   ;
            hfteis   = lfteis ;
            lfteis   = NULL   ;
        }
        else
        {
            Real2DArray_VectorMultiply ( False, 1.0e+00, itransformation, lfcore1b, 0.0e+00, mfcore1b, NULL ) ;
            hfteis   = Real2DArray_Allocate ( ni, nj, NULL ) ;
            Real2DArray_MatrixMultiply ( False, False, 1.0e+00, itransformation, lfteis, 0.0e+00, hfteis, NULL ) ;
        }

        if ( jtransformation == NULL )
        {
            lfcore2a = NULL     ;
            mfteis   = hfteis   ;
            hfteis   = NULL     ;
        }
        else
        {
            Real2DArray_VectorMultiply ( False, 1.0e+00, jtransformation, lfcore2a, 0.0e+00, mfcore2a, NULL ) ;
            mfteis   = Real2DArray_Allocate ( ni, nj, NULL ) ;
            Real2DArray_MatrixMultiply ( False, True, 1.0e+00, hfteis, jtransformation, 0.0e+00, mfteis, NULL ) ;
        }

        /* . Determine the TEI indices. */
        /* . There is no restriction on i, j, k and l as their order is checked when building the Fock matrices. */
        n = 0 ;
        for ( i = 0 ; i < idata->norbitals ; i++ )
        {
            for ( j = 0 ; j <= i ; j++ )
            {
                for ( k = 0 ; k < jdata->norbitals ; k++ )
                {
                    for ( l = 0 ; l <= k ; l++ )
                    {
                        indices[n  ] = i + istart ;
                        indices[n+1] = j + istart ;
                        indices[n+2] = k + jstart ;
                        indices[n+3] = l + jstart ;
                        n += 4 ;
                    }
                }
            }
        }

# ifdef PRINTMOPACTEIS
{
/* . Printing. */
printf ( "\n\nTWO-ELECTRON INTEGRALS:\n%d %d\n", idata->atomicNumber, jdata->atomicNumber ) ;
n = 0 ;
for ( i = n = 0 ; i < idata->norbitals ; i++ )
{
    for ( j = 0 ; j <= i ; j++ )
    {
        for ( k = 0 ; k < jdata->norbitals ; k++ )
        {
            for ( l = 0 ; l <= k ; l++, n++ ) printf ( "%5d %5d %5d %5d %5d %20.10f\n", n, i, j, k, l, mfteis->data[n] ) ;
        }
    }
}
printf ( "\n\n" ) ;
}
# endif

        /* . Save the data. */
        BlockStorage_Data_Add ( twoelectronintegrals, ( ni * nj ), mfteis->data, indices, NULL ) ;

        /* . Compute the core-core repulsion terms. */
        (*enuc) = MNDOIntegralUtilities_CoreCore ( idata, jdata, r ) ;

        /* . Finish up. */
        Real2DArray_Deallocate ( &hfteis   ) ;
        Real2DArray_Deallocate ( &lfteis   ) ;
        Real2DArray_Deallocate ( &mfteis   ) ;
        Real1DArray_Deallocate ( &lfcore1b ) ;
        Real1DArray_Deallocate ( &lfcore2a ) ;
        /* . The transformations can be identical so be careful about deallocation. */
        if ( itransformation != jtransformation ) Real2DArray_Deallocate ( &itransformation ) ;
        Real2DArray_Deallocate ( &jtransformation ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the derivatives in the molecular frame.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MNDOIntegrals_MolecularFrame2CIntegralsD ( const MNDOParameters *idata, const Integer ifirst, const Real *xi, const MNDOParameters *jdata, const Integer jfirst, const Real *xj,
                                                                                               const Real1DArray *dOneI, const Real1DArray *dOneJ, const Real2DArray *dTwoIJ, Real *eng )
{
    Boolean         doI, doJ ;
    Real       dcore0, dcore0f, doei0, doei0f, doei1, dtei0, dtei0f, dtei1, factor, r, scaling[3], x, y, z ;
    Integer          i, ix, j, ni, nj ;
    Real2DArray *dhfteis = NULL, *dlfteis = NULL, *lfteis = NULL, *LTj = NULL, *temporaryij = NULL, *TiL = NULL ;
    Real2DArray *itransformation = NULL, *itransformationD[3], *itransformationX = NULL, *itransformationY = NULL, *itransformationZ = NULL,
                *jtransformation = NULL, *jtransformationD[3], *jtransformationX = NULL, *jtransformationY = NULL, *jtransformationZ = NULL ;
    Real1DArray *dlfcore1b = NULL, *dlfcore2a = NULL, *lfcore1b = NULL, *lfcore2a = NULL, *temporaryi = NULL, *temporaryj = NULL ;

    /* . Initialization. */
    for ( i = 0 ; i < 3; i++ ) eng[i] = 0.0e+00;

    /* . Get the transformation matrices. */
    r = MNDOIntegralUtilities_GetTransformationMatrices ( idata, jdata, xi, xj, &itransformation, &jtransformation, &itransformationX, &itransformationY, &itransformationZ, &jtransformationX, &jtransformationY, &jtransformationZ, &x, &y, &z ) ;
    if ( r < SMALL_RIJ ) return ;

    /* . Allocate space. */
    ni = ( idata->norbitals * ( idata->norbitals + 1 ) ) / 2 ;
    nj = ( jdata->norbitals * ( jdata->norbitals + 1 ) ) / 2 ;
    dlfteis   = Real2DArray_Allocate ( ni, nj, NULL ) ; Real2DArray_Set ( dlfteis   , 0.0e+00 ) ;
    lfteis    = Real2DArray_Allocate ( ni, nj, NULL ) ; Real2DArray_Set ( lfteis    , 0.0e+00 ) ;
    dlfcore1b = Real1DArray_Allocate ( ni    , NULL ) ; Real1DArray_Set ( dlfcore1b , 0.0e+00 ) ;
    dlfcore2a = Real1DArray_Allocate ( nj    , NULL ) ; Real1DArray_Set ( dlfcore2a , 0.0e+00 ) ;
    lfcore1b  = Real1DArray_Allocate ( ni    , NULL ) ; Real1DArray_Set ( lfcore1b  , 0.0e+00 ) ;
    lfcore2a  = Real1DArray_Allocate ( nj    , NULL ) ; Real1DArray_Set ( lfcore2a  , 0.0e+00 ) ;

    /* . Define some transformation factors. */
    scaling[0] = x / r ; scaling[1] = y / r ; scaling[2] = z / r ;
    itransformationD[0] = itransformationX ; itransformationD[1] = itransformationY ; itransformationD[2] = itransformationZ ;
    jtransformationD[0] = jtransformationX ; jtransformationD[1] = jtransformationY ; jtransformationD[2] = jtransformationZ ;

    /* . Compute the integrals and derivatives in the local frame. */
    MNDOIntegralUtilities_LocalFrame2CTEIs ( idata, jdata, r, lfteis, lfcore1b, lfcore2a, dlfteis, dlfcore1b, dlfcore2a ) ;

    /* . Set some flags. */
    doI = ( itransformation != NULL ) ;
    doJ = ( jtransformation != NULL ) ;

    /* . Local frame terms (only depend on r). */
    /* . Core. */
    dcore0f = MNDOIntegralUtilities_CoreCoreD ( idata, jdata, r ) ;

    /* . OEIs. */
    doei0f = 0.0e+00 ;
    if ( doI )
    {
        temporaryi = Real1DArray_Allocate ( ni, NULL ) ;
        Real2DArray_VectorMultiply ( False, 1.0e+00, itransformation, dlfcore1b, 0.0e+00, temporaryi, NULL ) ;
    }
    else
    {
        temporaryi = dlfcore1b ;
        dlfcore1b  = NULL      ;
    }
    factor = Real1DArray_Dot ( dOneI, temporaryi, NULL ) ; doei0f += factor ;
    if ( doJ )
    {
        temporaryj = Real1DArray_Allocate ( nj, NULL ) ;
        Real2DArray_VectorMultiply ( False, 1.0e+00, jtransformation, dlfcore2a, 0.0e+00, temporaryj, NULL ) ;
    }
    else
    {
        temporaryj = dlfcore2a ;
        dlfcore2a  = NULL      ;
    }
    factor = Real1DArray_Dot ( dOneJ, temporaryj, NULL ) ; doei0f += factor ;

    /* . TEIs. */
    dtei0f = 0.0e+00 ;
    if ( doI )
    {
        dhfteis = Real2DArray_Allocate ( ni, nj, NULL ) ;
        Real2DArray_MatrixMultiply ( False, False, 1.0e+00, itransformation, dlfteis, 0.0e+00, dhfteis, NULL ) ;
    }
    else
    {
        dhfteis = dlfteis ;
        dlfteis = NULL    ;
    }
    if ( doJ )
    {
        temporaryij = Real2DArray_Allocate ( ni, nj, NULL ) ;
        Real2DArray_MatrixMultiply ( False, True, 1.0e+00, dhfteis, jtransformation, 0.0e+00, temporaryij, NULL ) ;
    }
    else
    {
        temporaryij = dhfteis ;
        dhfteis     = NULL   ;
    }

    /* . Coulomb and exchange terms. */
    for ( i = 0 ; i < ni ; i++ )
    {
        for ( j = 0 ; j < nj ; j++ ) dtei0f -= Real2DArray_Item ( dTwoIJ, i, j ) * Real2DArray_Item ( temporaryij, i, j ) ;
    }

    /* . Determine some intermediate matrices for later. */
    if ( doI )
    {
        LTj = Real2DArray_Allocate ( ni, nj, NULL ) ;
        if ( doJ ) Real2DArray_MatrixMultiply ( False, True, 1.0e+00, lfteis, jtransformation, 0.0e+00, LTj, NULL ) ;
        else       Real2DArray_CopyTo ( lfteis, LTj, NULL ) ;
    }
    if ( doJ )
    {
        TiL = Real2DArray_Allocate ( ni, nj, NULL ) ;
        if ( doI ) Real2DArray_MatrixMultiply ( False, False, 1.0e+00, itransformation, lfteis, 0.0e+00, TiL, NULL ) ;
        else       Real2DArray_CopyTo ( lfteis, TiL, NULL ) ;
    }

    /* . Electronic terms. */
    /* . Loop over the Cartesian components. */
    for ( ix = 0 ; ix < 3 ; ++ix )
    {
        /* . Core. */
        dcore0 = - scaling[ix] * dcore0f ;

        /* . OEI. */
        doei0 = - scaling[ix] * doei0f ;
        doei1 = 0.0e+00 ;
        if ( doI )
        {
            Real2DArray_VectorMultiply ( False, 1.0e+00, itransformationD[ix], lfcore1b, 0.0e+00, temporaryi, NULL ) ;
            factor = Real1DArray_Dot ( dOneI, temporaryi, NULL ) ; doei1 -= factor ;
        }
        if ( doJ )
        {
            Real2DArray_VectorMultiply ( False, 1.0e+00, jtransformationD[ix], lfcore2a, 0.0e+00, temporaryj, NULL ) ;
            factor = Real1DArray_Dot ( dOneJ, temporaryj, NULL ) ; doei1 -= factor ;
        }

        /* . TEI. */
        dtei0 = - scaling[ix] * dtei0f ;
        dtei1 = 0.0e+00 ;
        if ( doI || doJ )
        {
            /* . Get the total integrals. */
            factor = 0.0e+00 ;
            if ( doI ) { Real2DArray_MatrixMultiply ( False, False, 1.0e+00, itransformationD[ix], LTj, 0.0e+00, temporaryij, NULL ) ; factor = 1.0e+00 ; }
            if ( doJ ) { Real2DArray_MatrixMultiply ( False, True,   1.0e+00, TiL, jtransformationD[ix], factor,  temporaryij, NULL ) ; }

            /* . Coulomb and exchange terms. */
            for ( i = 0 ; i < ni ; i++ )
            {
                for ( j = 0 ; j < nj ; j++ ) dtei1 += Real2DArray_Item ( dTwoIJ, i, j ) * Real2DArray_Item ( temporaryij, i, j ) ;
            }
        }

        /* . Save the gradient terms. */
        eng[ix] = dcore0 + doei0 + doei1 + dtei0 + dtei1 ;
# ifdef DEBUGMNDOCIGRADIENTS
 printf ( "\nGradient contributions (%d,%d), component %d, nuclear, oei, tei: %20.10f %20.10f %20.10f\n", ifirst, jfirst, ix, dcore0, doei0+doei1, dtei0+dtei1 ) ;
# endif
    }

    /* . Clear up. */
    Real2DArray_Deallocate ( &dhfteis     ) ;
    Real2DArray_Deallocate ( &dlfteis     ) ;
    Real2DArray_Deallocate ( &lfteis      ) ;
    Real2DArray_Deallocate ( &LTj         ) ;
    Real2DArray_Deallocate ( &TiL         ) ;
    Real2DArray_Deallocate ( &temporaryij ) ;
    Real1DArray_Deallocate ( &dlfcore1b   ) ;
    Real1DArray_Deallocate ( &dlfcore2a   ) ;
    Real1DArray_Deallocate ( &lfcore1b    ) ;
    Real1DArray_Deallocate ( &lfcore2a    ) ;
    Real1DArray_Deallocate ( &temporaryi  ) ;
    Real1DArray_Deallocate ( &temporaryj  ) ;

    /* . The transformations can be identical so be careful about deallocation. */
    if ( itransformation  != jtransformation  ) Real2DArray_Deallocate ( &itransformation  ) ;
    if ( itransformationX != jtransformationX ) Real2DArray_Deallocate ( &itransformationX ) ;
    if ( itransformationY != jtransformationY ) Real2DArray_Deallocate ( &itransformationY ) ;
    if ( itransformationZ != jtransformationZ ) Real2DArray_Deallocate ( &itransformationZ ) ;
    Real2DArray_Deallocate ( &jtransformation  ) ;
    Real2DArray_Deallocate ( &jtransformationX ) ;
    Real2DArray_Deallocate ( &jtransformationY ) ;
    Real2DArray_Deallocate ( &jtransformationZ ) ;
}
