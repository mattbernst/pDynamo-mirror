/*------------------------------------------------------------------------------
! . File      : MNDOIntegralUtilities.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . Utility procedures for calculating the integrals in a MNDO method.
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
# include "MNDOIntegralDefinitions.h"
# include "MNDOIntegralUtilities.h"
# include "MNDOParameters.h"
# include "Units.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Real MNDOIntegralUtilities_2CChargeInteraction  ( const Real r, const Integer l1, const Integer l2, const Integer m, const Real da, const Real db, const Real add ) ;
static Real MNDOIntegralUtilities_2CChargeInteractionD ( const Real r, const Integer l1, const Integer l2, const Integer m, const Real da, const Real db, const Real add ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the core-core terms.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real MNDOIntegralUtilities_CoreCore ( const MNDOParameters *idata, const MNDOParameters *jdata, const Real rij )
{
    Real  aij, ax, d, enuc, gam, scale, xij, zaf, zbf ;
    Integer     i, j, nt ;

    /* . Initialization. */
    enuc = 0.0e+00 ;

    /* . Calculate the core integral. */
    gam = 1.0e+00 / sqrt ( rij*rij + pow ( ( idata->po[8] + jdata->po[8] ), 2 ) ) ;

    /* . Diatomic terms. */
    if ( idata->QDIATOMIC && jdata->QDIATOMIC )
    {
        /* . Initialization. */
        scale = 1.0e+00 ;

        /* . Determine aij and xij parameter-dependent terms. */
        if ( ( jdata->atomicNumber < idata->ndiatomic ) && idata->QDIATOMICFLAGS[jdata->atomicNumber] )
        {
            aij =           idata->diatomica[jdata->atomicNumber] ;
            xij = 2.0e+00 * idata->diatomicx[jdata->atomicNumber] ; /* . Factor of 2 - see Stewart's am1/d paper! */

            /* . N-H and O-H. */
            if ( ( ( idata->atomicNumber == 1 ) && ( ( jdata->atomicNumber == 6 ) || ( jdata->atomicNumber == 7 ) || ( jdata->atomicNumber == 8 ) ) ) ||
                 ( ( jdata->atomicNumber == 1 ) && ( ( idata->atomicNumber == 6 ) || ( idata->atomicNumber == 7 ) || ( idata->atomicNumber == 8 ) ) ) )
            {
                scale += xij * exp ( - aij * rij * rij * UNITS_LENGTH_BOHRS_TO_ANGSTROMS ) ;
            }
            /* . All others. */
            else
            {
                scale += xij * exp ( - aij * rij * ( 1.0e+00 + 0.0003e+00 * pow ( rij * UNITS_LENGTH_BOHRS_TO_ANGSTROMS, 5 ) ) ) ;
            }
/*printf ( "TEST> %5d %5d %25.5f %25.5f %25.5f %25.5f\n", idata->atomicNumber, jdata->atomicNumber, rij, aij, xij, scale ) ;*/
        }
/* . Removed as this is already flagged elsewhere. */
/*
        else printf ( "\nMNDOPARAMETERS_CORECORE> Missing diatomic pair parameters: %d %d\n", idata->atomicNumber, jdata->atomicNumber ) ;
*/

        /* . Element-specific extra terms independent of aij and xij. */
        /* . C-C. */
        if ( ( idata->atomicNumber == 6 ) && ( jdata->atomicNumber == 6 ) )
        {
            scale += 9.28e+00 * exp ( - 5.98e+00 * rij * UNITS_LENGTH_BOHRS_TO_ANGSTROMS ) ;
/*printf ( "CCTEST> %5d %5d %25.5f %25.5f\n", idata->atomicNumber, jdata->atomicNumber, rij, 9.28e+00 * exp ( - 5.98e+00 * rij * UNITS_LENGTH_BOHRS_TO_ANGSTROMS ) ) ;*/
        }
        /* . Si-O. */
        if ( ( ( idata->atomicNumber ==  8 ) && ( jdata->atomicNumber == 14 ) ) ||
             ( ( idata->atomicNumber == 14 ) && ( jdata->atomicNumber ==  8 ) ) )
        {
            scale -= 0.0007e+00 * exp ( - pow ( rij * UNITS_LENGTH_BOHRS_TO_ANGSTROMS - 2.9e+00, 2 ) ) ;
        }

        /* . Initial term. */
        enuc = idata->zcore * jdata->zcore * gam * scale ;

        /* . Unpolarizable core. */
        enuc += PM6_UNPOLARIZABLECORE * pow ( ( ( pow ( idata->zcore, 1.0e+00 / 3.0e+00 ) + pow ( jdata->zcore, 1.0e+00 / 3.0e+00 ) ) / ( rij * UNITS_LENGTH_BOHRS_TO_ANGSTROMS ) ), 12 ) ;
/*printf ( "UNTEST> %5d %5d %25.5f %25.5f\n", idata->atomicNumber, jdata->atomicNumber, rij, PM6_UNPOLARIZABLECORE * pow ( ( ( pow ( idata->zcore, 1.0e+00 / 3.0e+00 ) + pow ( jdata->zcore, 1.0e+00 / 3.0e+00 ) ) / ( rij * UNITS_LENGTH_BOHRS_TO_ANGSTROMS ) ), 12 ) ) ;*/
        scale    = 0.0e+00 ;
    }
    /* . Monoatomic terms. */
    else
    {
        scale = exp ( - idata->alp * rij ) + exp( - jdata->alp * rij ) ;
        nt    = idata->atomicNumber + jdata->atomicNumber ;
        if ( ( nt == 8 ) || ( nt == 9 ) )
        {
            if ( ( idata->atomicNumber == 7 ) || ( idata->atomicNumber == 8 ) ) scale += ( UNITS_LENGTH_BOHRS_TO_ANGSTROMS * rij - 1.0e+00 ) * exp( - idata->alp * rij ) ;
            if ( ( jdata->atomicNumber == 7 ) || ( jdata->atomicNumber == 8 ) ) scale += ( UNITS_LENGTH_BOHRS_TO_ANGSTROMS * rij - 1.0e+00 ) * exp( - jdata->alp * rij ) ;
        }
        enuc  = idata->zcore * jdata->zcore * gam ;
        scale = fabs ( scale * enuc ) ;
    }

    /* . Compute the AM1/PM3-specific terms. */
    for ( i = 0 ; i < idata->nam1pm3g ; i++ )
    {
       d  = rij - idata->fn3[i] ;
       ax = idata->fn2[i] * d * d ;
       if ( ax <= EXPONENT_TOLERANCE ) scale += idata->gphot * jdata->gphot * idata->zcore * jdata->zcore / rij * idata->fn1[i] * exp( -ax ) ;
    }
    for ( i = 0 ; i < jdata->nam1pm3g ; i++ )
    {
       d  = rij - jdata->fn3[i] ;
       ax = jdata->fn2[i] * d * d ;
       if ( ax <= EXPONENT_TOLERANCE ) scale += idata->gphot * jdata->gphot * idata->zcore * jdata->zcore / rij * jdata->fn1[i] * exp( -ax ) ;
    }

    /* . Compute the PDDG-specific terms. */
    if ( ( idata->npddg > 0 ) && ( jdata->npddg > 0 ) )
    {
       zaf = idata->zcore / ( idata->zcore + jdata->zcore ) ;
       zbf = jdata->zcore / ( idata->zcore + jdata->zcore ) ;
       for ( i = 0 ; i < idata->npddg ; i++ )
       {
          for ( j = 0 ; j < jdata->npddg ; j++ )
          {
             d  = rij - idata->pddge[i] - jdata->pddge[j] ;
             ax = PDDG_EXPONENT * d * d ;
             scale += ( zaf * idata->pddgc[i] + zbf * jdata->pddgc[j] ) * exp ( -ax ) ;
          }
       }
    }

    /* . Add in the correction factor. */
    enuc += scale ;

    /* . Finish up. */
    return enuc ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the core-core derivative with respect to r.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real MNDOIntegralUtilities_CoreCoreD ( const MNDOParameters *idata, const MNDOParameters *jdata, const Real rij )
{
    Real anam1, aij, ax, dd, dedr, dgdr, dscale, exi, exj, f3, gam, scale, xij, zaf, zbf ;
    Integer    i, j ;

    /* . Initialization. */
    dedr = 0.0e+00 ;

    /* . Calculate the core integral and its derivative. */
    gam  = 1.0e+00 / sqrt ( rij*rij + pow ( ( idata->po[8] + jdata->po[8] ), 2 ) ) ;
    dgdr = - rij * gam * gam * gam ;

    /* . Diatomic terms. */
    if ( idata->QDIATOMIC && jdata->QDIATOMIC )
    {
        /* . Initialization. */
        dscale = 0.0e+00 ;
        scale  = 1.0e+00 ;

        /* . Determine aij and xij parameter-dependent terms. */
        if ( ( jdata->atomicNumber < idata->ndiatomic ) && idata->QDIATOMICFLAGS[jdata->atomicNumber] )
        {
            aij =           idata->diatomica[jdata->atomicNumber] ;
            xij = 2.0e+00 * idata->diatomicx[jdata->atomicNumber] ; /* . Factor of 2 - see Stewart's am1/d paper! */

            /* . C-H, N-H and O-H. */
            if ( ( ( idata->atomicNumber == 1 ) && ( ( jdata->atomicNumber == 6 ) || ( jdata->atomicNumber == 7 ) || ( jdata->atomicNumber == 8 ) ) ) ||
                 ( ( jdata->atomicNumber == 1 ) && ( ( idata->atomicNumber == 6 ) || ( idata->atomicNumber == 7 ) || ( idata->atomicNumber == 8 ) ) ) )
            {
                f3 = xij * exp ( - aij * rij * rij * UNITS_LENGTH_BOHRS_TO_ANGSTROMS ) ;
                dscale -= 2.0e+00 * aij * rij * UNITS_LENGTH_BOHRS_TO_ANGSTROMS * f3 ;
                scale  += f3 ;
            }
            /* . All others. */
            else
            {
                dd = 0.0003e+00 * pow ( rij * UNITS_LENGTH_BOHRS_TO_ANGSTROMS, 5 ) ;
                f3 = xij * exp ( - aij * rij * ( 1.0e+00 + dd ) ) ;
                dscale -= aij * ( 1.0e+00 + 6.0e+00 * dd ) * f3 ;
                scale  += f3 ;
            }
        }

        /* . Element-specific extra terms independent of aij and xij. */
        /* . C-C. */
        if ( ( idata->atomicNumber == 6 ) && ( jdata->atomicNumber == 6 ) )
        {
            f3 = 9.28e+00 * exp ( - 5.98e+00 * rij * UNITS_LENGTH_BOHRS_TO_ANGSTROMS ) ;
            dscale -= 5.98e+00 * UNITS_LENGTH_BOHRS_TO_ANGSTROMS * f3 ;
            scale  += f3 ;
        }
        /* . Si-O. */
        if ( ( ( idata->atomicNumber ==  8 ) && ( jdata->atomicNumber == 14 ) ) ||
             ( ( idata->atomicNumber == 14 ) && ( jdata->atomicNumber ==  8 ) ) )
        {
            dd = rij * UNITS_LENGTH_BOHRS_TO_ANGSTROMS - 2.9e+00 ;
            f3 = - 0.0007e+00 * exp ( - pow ( dd, 2 ) ) ;
            dscale -= 2.0e+00 * dd * UNITS_LENGTH_BOHRS_TO_ANGSTROMS * f3 ;
            scale  += f3 ;
        }

        /* . Initial term. */
        dedr = idata->zcore * jdata->zcore * ( dgdr * scale + gam * dscale ) ;

        /* . Unpolarizable core. */
        dedr -= 12.0e+00 * PM6_UNPOLARIZABLECORE * pow ( ( ( pow ( idata->zcore, 1.0e+00 / 3.0e+00 ) + pow ( jdata->zcore, 1.0e+00 / 3.0e+00 ) ) / ( rij * UNITS_LENGTH_BOHRS_TO_ANGSTROMS ) ), 12 ) / rij ;
    }
    /* . Monoatomic terms. */
    else
    {
        exi = exp ( -idata->alp * rij ) ;
        exj = exp ( -jdata->alp * rij ) ;
	if ( ( idata->atomicNumber == 1 ) && ( ( jdata->atomicNumber == 7 ) || ( jdata->atomicNumber == 8 ) ) )
        {
	    f3 = 1.0e+00 + exi + UNITS_LENGTH_BOHRS_TO_ANGSTROMS * rij * exj ;
	    dd = dgdr * f3 - gam * ( idata->alp * exi + UNITS_LENGTH_BOHRS_TO_ANGSTROMS * ( jdata->alp * rij - 1.0e+00 ) * exj ) ;
	}
        else if ( ( ( idata->atomicNumber == 7 ) || ( idata->atomicNumber == 8 ) ) && ( jdata->atomicNumber == 1 ) )
        {
	    f3 = 1.0e+00 + exj + UNITS_LENGTH_BOHRS_TO_ANGSTROMS * rij * exi ;
	    dd = dgdr * f3 - gam * ( jdata->alp * exj + UNITS_LENGTH_BOHRS_TO_ANGSTROMS * ( idata->alp * rij - 1.0e+00 ) * exi ) ;
	}
        else
        {
	    f3 = 1.0e+00 + exi + exj ;
            dd = dgdr * f3 - gam * ( idata->alp * exi + jdata->alp * exj ) ;
	}
	dedr = idata->zcore * jdata->zcore * dd ;
    }

    /* . Compute the AM1/PM3-specific terms. */
    anam1 = 0.0e+00 ;
    for ( i = 0 ; i < idata->nam1pm3g ; i++ )
    {
        dd = rij - idata->fn3[i] ;
        ax = idata->fn2[i] * dd * dd ;
        if ( ax <= EXPONENT_TOLERANCE ) anam1 += idata->fn1[i] * ( 1.0e+00 / ( rij * rij ) + idata->fn2[i] * 2.0e+00 * dd / rij ) * exp( -ax ) ;
    }
    for ( i = 0 ; i < jdata->nam1pm3g ; i++ )
    {
        dd = rij - jdata->fn3[i] ;
        ax = jdata->fn2[i] * dd * dd ;
        if ( ax <= EXPONENT_TOLERANCE ) anam1 += jdata->fn1[i] * ( 1.0e+00 / ( rij * rij ) + jdata->fn2[i] * 2.0e+00 * dd / rij ) * exp( -ax ) ;
    }
    dedr -= anam1 * idata->gphot * jdata->gphot * idata->zcore * jdata->zcore ;

    /* . Compute the PDDG-specific terms. */
    if ( ( idata->npddg > 0 ) && ( jdata->npddg > 0 ) )
    {
        anam1 = 0.0e+00 ;
        zaf   = idata->zcore / ( idata->zcore + jdata->zcore ) ;
        zbf   = jdata->zcore / ( idata->zcore + jdata->zcore ) ;
        for ( i = 0 ; i < idata->npddg ; i++ )
        {
            for ( j = 0 ; j < jdata->npddg ; j++ )
            {
                dd = rij - idata->pddge[i] - jdata->pddge[j] ;
                ax = PDDG_EXPONENT * dd * dd ;
                anam1 += ( zaf * idata->pddgc[i] + zbf * jdata->pddgc[j] ) * 2.0e+00 * PDDG_EXPONENT * dd * exp ( -ax ) ;
            }
        }
        dedr -= anam1 ;
    }

    /* . Finish up. */
    return dedr ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the core-core derivative with respect to r.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MNDOIntegralUtilities_CoreCoreFD ( const MNDOParameters *idata, const MNDOParameters *jdata, const Real rij, Real *fCore, Real *dCore )
{
    Integer i, j ;
    Real    anam1, aij, ax, dd, dedr, dgdr, dscale, enuc, exi, exj, f3, gam, scale, xij, zaf, zbf ;

    /* . Initialization. */
    enuc = 0.0e+00 ;
    dedr = 0.0e+00 ;

    /* . Calculate the core integral and its derivative. */
    gam  = 1.0e+00 / sqrt ( rij*rij + pow ( ( idata->po[8] + jdata->po[8] ), 2 ) ) ;
    dgdr = - rij * gam * gam * gam ;

    /* . Diatomic terms. */
    if ( idata->QDIATOMIC && jdata->QDIATOMIC )
    {
        /* . Initialization. */
        dscale = 0.0e+00 ;
        scale  = 1.0e+00 ;

        /* . Determine aij and xij parameter-dependent terms. */
        if ( ( jdata->atomicNumber < idata->ndiatomic ) && idata->QDIATOMICFLAGS[jdata->atomicNumber] )
        {
            aij =           idata->diatomica[jdata->atomicNumber] ;
            xij = 2.0e+00 * idata->diatomicx[jdata->atomicNumber] ; /* . Factor of 2 - see Stewart's am1/d paper! */

            /* . C-H, N-H and O-H. */
            if ( ( ( idata->atomicNumber == 1 ) && ( ( jdata->atomicNumber == 6 ) || ( jdata->atomicNumber == 7 ) || ( jdata->atomicNumber == 8 ) ) ) ||
                 ( ( jdata->atomicNumber == 1 ) && ( ( idata->atomicNumber == 6 ) || ( idata->atomicNumber == 7 ) || ( idata->atomicNumber == 8 ) ) ) )
            {
                f3 = xij * exp ( - aij * rij * rij * UNITS_LENGTH_BOHRS_TO_ANGSTROMS ) ;
                dscale -= 2.0e+00 * aij * rij * UNITS_LENGTH_BOHRS_TO_ANGSTROMS * f3 ;
                scale  += f3 ;
            }
            /* . All others. */
            else
            {
                dd = 0.0003e+00 * pow ( rij * UNITS_LENGTH_BOHRS_TO_ANGSTROMS, 5 ) ;
                f3 = xij * exp ( - aij * rij * ( 1.0e+00 + dd ) ) ;
                dscale -= aij * ( 1.0e+00 + 6.0e+00 * dd ) * f3 ;
                scale  += f3 ;
            }
        }

        /* . Element-specific extra terms independent of aij and xij. */
        /* . C-C. */
        if ( ( idata->atomicNumber == 6 ) && ( jdata->atomicNumber == 6 ) )
        {
            f3 = 9.28e+00 * exp ( - 5.98e+00 * rij * UNITS_LENGTH_BOHRS_TO_ANGSTROMS ) ;
            dscale -= 5.98e+00 * UNITS_LENGTH_BOHRS_TO_ANGSTROMS * f3 ;
            scale  += f3 ;
        }
        /* . Si-O. */
        if ( ( ( idata->atomicNumber ==  8 ) && ( jdata->atomicNumber == 14 ) ) ||
             ( ( idata->atomicNumber == 14 ) && ( jdata->atomicNumber ==  8 ) ) )
        {
            dd = rij * UNITS_LENGTH_BOHRS_TO_ANGSTROMS - 2.9e+00 ;
            f3 = - 0.0007e+00 * exp ( - pow ( dd, 2 ) ) ;
            dscale -= 2.0e+00 * dd * UNITS_LENGTH_BOHRS_TO_ANGSTROMS * f3 ;
            scale  += f3 ;
        }

        /* . Initial term. */
        enuc   = idata->zcore * jdata->zcore * gam * scale ;
        dedr   = idata->zcore * jdata->zcore * ( dgdr * scale + gam * dscale ) ;

        /* . Unpolarizable core. */
        enuc  += PM6_UNPOLARIZABLECORE * pow ( ( ( pow ( idata->zcore, 1.0e+00 / 3.0e+00 ) + pow ( jdata->zcore, 1.0e+00 / 3.0e+00 ) ) / ( rij * UNITS_LENGTH_BOHRS_TO_ANGSTROMS ) ), 12 ) ;
        dedr  -= 12.0e+00 * PM6_UNPOLARIZABLECORE * pow ( ( ( pow ( idata->zcore, 1.0e+00 / 3.0e+00 ) + pow ( jdata->zcore, 1.0e+00 / 3.0e+00 ) ) / ( rij * UNITS_LENGTH_BOHRS_TO_ANGSTROMS ) ), 12 ) / rij ;
        scale  = 0.0e+00 ;
    }
    /* . Monoatomic terms. */
    else
    {
        exi   = exp ( -idata->alp * rij ) ;
        exj   = exp ( -jdata->alp * rij ) ;
        scale = exi + exj ;
	if ( ( idata->atomicNumber == 1 ) && ( ( jdata->atomicNumber == 7 ) || ( jdata->atomicNumber == 8 ) ) )
        {
	    f3     = 1.0e+00 + exi + UNITS_LENGTH_BOHRS_TO_ANGSTROMS * rij * exj ;
	    dd     = dgdr * f3 - gam * ( idata->alp * exi + UNITS_LENGTH_BOHRS_TO_ANGSTROMS * ( jdata->alp * rij - 1.0e+00 ) * exj ) ;
            scale += ( UNITS_LENGTH_BOHRS_TO_ANGSTROMS * rij - 1.0e+00 ) * exj ;
	}
        else if ( ( ( idata->atomicNumber == 7 ) || ( idata->atomicNumber == 8 ) ) && ( jdata->atomicNumber == 1 ) )
        {
	    f3     = 1.0e+00 + exj + UNITS_LENGTH_BOHRS_TO_ANGSTROMS * rij * exi ;
	    dd     = dgdr * f3 - gam * ( jdata->alp * exj + UNITS_LENGTH_BOHRS_TO_ANGSTROMS * ( idata->alp * rij - 1.0e+00 ) * exi ) ;
            scale += ( UNITS_LENGTH_BOHRS_TO_ANGSTROMS * rij - 1.0e+00 ) * exi ;
	}
        else
        {
	    f3     = 1.0e+00 + exi + exj ;
            dd     = dgdr * f3 - gam * ( idata->alp * exi + jdata->alp * exj ) ;
	}
	dedr  = idata->zcore * jdata->zcore * dd ;
        enuc  = idata->zcore * jdata->zcore * gam ;
        scale = fabs ( scale * enuc ) ;
    }

    /* . Compute the AM1/PM3-specific terms. */
    anam1 = 0.0e+00 ;
    for ( i = 0 ; i < idata->nam1pm3g ; i++ )
    {
        dd = rij - idata->fn3[i] ;
        ax = idata->fn2[i] * dd * dd ;
        if ( ax <= EXPONENT_TOLERANCE )
        {
            anam1 += idata->fn1[i] * ( 1.0e+00 / ( rij * rij ) + idata->fn2[i] * 2.0e+00 * dd / rij ) * exp( -ax ) ;
            scale += idata->gphot * jdata->gphot * idata->zcore * jdata->zcore / rij * idata->fn1[i] * exp( -ax ) ;
        }
    }
    for ( i = 0 ; i < jdata->nam1pm3g ; i++ )
    {
        dd = rij - jdata->fn3[i] ;
        ax = jdata->fn2[i] * dd * dd ;
        if ( ax <= EXPONENT_TOLERANCE )
        {
            anam1 += jdata->fn1[i] * ( 1.0e+00 / ( rij * rij ) + jdata->fn2[i] * 2.0e+00 * dd / rij ) * exp( -ax ) ;
            scale += idata->gphot * jdata->gphot * idata->zcore * jdata->zcore / rij * jdata->fn1[i] * exp( -ax ) ;
        }
    }
    dedr -= anam1 * idata->gphot * jdata->gphot * idata->zcore * jdata->zcore ;

    /* . Compute the PDDG-specific terms. */
    if ( ( idata->npddg > 0 ) && ( jdata->npddg > 0 ) )
    {
        anam1 = 0.0e+00 ;
        zaf   = idata->zcore / ( idata->zcore + jdata->zcore ) ;
        zbf   = jdata->zcore / ( idata->zcore + jdata->zcore ) ;
        for ( i = 0 ; i < idata->npddg ; i++ )
        {
            for ( j = 0 ; j < jdata->npddg ; j++ )
            {
                dd = rij - idata->pddge[i] - jdata->pddge[j] ;
                ax = PDDG_EXPONENT * dd * dd ;
                anam1 += ( zaf * idata->pddgc[i] + zbf * jdata->pddgc[j] ) * 2.0e+00 * PDDG_EXPONENT * dd * exp ( -ax ) ;
                scale += ( zaf * idata->pddgc[i] + zbf * jdata->pddgc[j] ) * exp ( -ax ) ;
            }
        }
        dedr -= anam1 ;
    }

    /* . Add in the correction factor. */
    enuc += scale ;

    /* . Finish up. */
    (*dCore) = dedr ;
    (*fCore) = enuc ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the transformation matrices for a given atom pair (i-j).
!---------------------------------------------------------------------------------------------------------------------------------*/
/*
! . Phi zero when yji = 0.
! . ca  = cos(phi)     , sa  = sin(phi)
! . cb  = cos(theta)   , sb  = sin(theta)
! . c2a = cos(2*phi)   , s2a = sin(2*phi)
! . c2b = cos(2*theta) , s2b = sin(2*phi)
!
! . There is a problem when atoms are aligned on z-axis as phi is undefined.
! . Therefore, axes are swapped.
!
! . Note that a NULL transformation implies the identity whereas a NULL derivative transformation implies a zero matrix!
*/
# define ALIGNMENT_TOLERANCE 0.99999999e+00
Real MNDOIntegralUtilities_GetTransformationMatrices ( const MNDOParameters *idata, const MNDOParameters *jdata, const Real *coordi, const Real *coordj, Real2DArray **itransformation, Real2DArray **jtransformation,
                                                                                                                                      Real2DArray **itransformationX, Real2DArray **itransformationY, Real2DArray **itransformationZ,
                                                                                                                                      Real2DArray **jtransformationX, Real2DArray **jtransformationY, Real2DArray **jtransformationZ,
                                                                                                                                                                                                    Real *x, Real *y, Real *z )
{
    Real r = 0.0e+00 ;

    /* . Initialization. */
    if ( x != NULL ) (*x) = 0.0e+00 ;
    if ( y != NULL ) (*y) = 0.0e+00 ;
    if ( z != NULL ) (*z) = 0.0e+00 ;

    /* . Check for essential data structures. */
    if ( ( idata != NULL ) && ( jdata != NULL ) && ( coordi != NULL ) && ( coordj != NULL ) && ( itransformation != NULL ) && ( jtransformation != NULL ) )
    {
        auto Boolean         axesSwapped, doGradients ;
        auto Real       r2, xji, yji, zji ;
        auto Integer          ni, nj, norbitals ;
        auto Real2DArray *it = NULL, *itX = NULL, *itY = NULL, *itZ = NULL, *jt = NULL, *jtX = NULL, *jtY = NULL, *jtZ = NULL ;

        /* . Check for gradients. */
        doGradients = ( itransformationX != NULL ) && ( itransformationY != NULL ) && ( itransformationZ != NULL ) && ( jtransformationX != NULL ) && ( jtransformationY != NULL ) && ( jtransformationZ != NULL ) ;

        /* . Get orbital information. */
        ni = idata->norbitals ;
        nj = jdata->norbitals ;
        norbitals = Maximum ( ni, nj ) ;

        /* . Basic distance computation. */
        xji = coordj[0] - coordi[0] ;
        yji = coordj[1] - coordi[1] ;
        zji = coordj[2] - coordi[2] ;
        r2  = xji * xji + yji * yji + zji * zji ;
        r   = sqrt ( r2 ) ;

        /* . Save displacements. */
        if ( x != NULL ) (*x) = xji ;
        if ( y != NULL ) (*y) = yji ;
        if ( z != NULL ) (*z) = zji ;

        /* . Only continue if there are p orbitals or higher and the distance is big enough. */
        if ( ( norbitals > 1 ) && ( r2 > SMALL_RIJ2 ) )
        {
# ifdef MNDODORBITALS
            auto Real  b, b2, ca, cb, c2a, c2b, sa, sb, s2a, s2b ;
# else
            auto Real  b, b2, ca, cb, sa, sb ;
# endif
            auto Real  caX = 0.0e+00, caY = 0.0e+00, cbX = 0.0e+00, cbY = 0.0e+00, cbZ = 0.0e+00,
                         saX = 0.0e+00, saY = 0.0e+00, sbX = 0.0e+00, sbY = 0.0e+00, sbZ = 0.0e+00 ;
            auto Integer     i, ij, j, m, mn, n ;
            auto Real2DArray *clesser         = NULL, *clesserX         = NULL, *clesserY         = NULL, *clesserZ         = NULL,
                        *ctransformation = NULL, *ctransformationX = NULL, *ctransformationY = NULL, *ctransformationZ = NULL,
                        *otransformation = NULL, *otransformationX = NULL, *otransformationY = NULL, *otransformationZ = NULL,
                        *stransformation = NULL ;

            /* . Allocate the orbital transformation matrix and initialize to the ss case. */
            otransformation = Real2DArray_Allocate ( norbitals, norbitals, NULL ) ;
            Real2DArray_Set ( otransformation, 0.0e+00 ) ;
            Real2DArray_Item ( otransformation, S, S ) = 1.0e+00 ;
            if ( doGradients )
            {
                otransformationX = Real2DArray_Allocate ( norbitals, norbitals, NULL ) ; Real2DArray_Set ( otransformationX, 0.0e+00 ) ;
                otransformationY = Real2DArray_Allocate ( norbitals, norbitals, NULL ) ; Real2DArray_Set ( otransformationY, 0.0e+00 ) ;
                otransformationZ = Real2DArray_Allocate ( norbitals, norbitals, NULL ) ; Real2DArray_Set ( otransformationZ, 0.0e+00 ) ;
            }

            /* . Check for z-axis alignment. */
            axesSwapped = ( fabs ( zji / r ) > ALIGNMENT_TOLERANCE ) ;
/*axesSwapped = True ;*/
# ifdef PRINTMOPACTEIS
{
/* . Printing. */
printf ( "\n\nSWAP DATA = %25.15f %25.15f %25.15f\n\n", zji, r, ALIGNMENT_TOLERANCE ) ;
if ( axesSwapped ) printf ( "\n\nTWO-ELECTRON INTEGRALS - AXES SWAPPED:\n%d %d\n"    , idata->atomicNumber, jdata->atomicNumber ) ;
else               printf ( "\n\nTWO-ELECTRON INTEGRALS - AXES NOT SWAPPED:\n%d %d\n", idata->atomicNumber, jdata->atomicNumber ) ;
}
# endif

            if ( axesSwapped )
            {
                b = xji ; xji = zji ; zji = b ;
                stransformation = Real2DArray_Allocate ( norbitals, norbitals, NULL ) ; Real2DArray_Set ( stransformation, 0.0e+00 ) ;
            }

            /* . p-orbital transformation. */
            b2 = xji * xji + yji * yji ;
            b  = sqrt ( b2 ) ;
            sb = b / r ;
            ca = xji / b ;
            sa = yji / b ;
            cb = zji / r   ;

            /* . Rotation matrix. */
            Real2DArray_Item ( otransformation,  PX, PSIGMA   ) = ca*sb ;
            Real2DArray_Item ( otransformation,  PX, PPIPLUS  ) = ca*cb ;
            Real2DArray_Item ( otransformation,  PX, PPIMINUS ) = -sa   ;
            Real2DArray_Item ( otransformation,  PY, PSIGMA   ) = sa*sb ;
            Real2DArray_Item ( otransformation,  PY, PPIPLUS  ) = sa*cb ;
            Real2DArray_Item ( otransformation,  PY, PPIMINUS ) = ca    ;
            Real2DArray_Item ( otransformation,  PZ, PSIGMA   ) = cb    ;
            Real2DArray_Item ( otransformation,  PZ, PPIPLUS  ) = -sb   ;
/*            Real2DArray_Item ( otransformation,  PZ, PPIMINUS ) = 0.0 ; */

            /* . Derivatives. */
            if ( doGradients )
            {
                /* . Various factors. */
                caX =   yji * yji / ( b2 * b ) ;
                caY = - xji * yji / ( b2 * b ) ;
                saX =   caY ;
                saY =   xji * xji / ( b2 * b ) ;
                cbX = - xji * zji       / (     r2 * r ) ;
                cbY = - yji * zji       / (     r2 * r ) ;
                cbZ =   b2              / (     r2 * r ) ;
                sbX =   xji * zji * zji / ( b * r2 * r ) ;
                sbY =   yji * zji * zji / ( b * r2 * r ) ;
                sbZ = - b * zji         / (     r2 * r ) ;

                /* . X. */
                Real2DArray_Item ( otransformationX,  PX, PSIGMA   ) = caX*sb + ca*sbX ;
                Real2DArray_Item ( otransformationX,  PX, PPIPLUS  ) = caX*cb + ca*cbX ;
                Real2DArray_Item ( otransformationX,  PX, PPIMINUS ) = -saX ;
                Real2DArray_Item ( otransformationX,  PY, PSIGMA   ) = saX*sb + sa*sbX ;
                Real2DArray_Item ( otransformationX,  PY, PPIPLUS  ) = saX*cb + sa*cbX ;
                Real2DArray_Item ( otransformationX,  PY, PPIMINUS ) = caX  ;
                Real2DArray_Item ( otransformationX,  PZ, PSIGMA   ) = cbX  ;
                Real2DArray_Item ( otransformationX,  PZ, PPIPLUS  ) = -sbX ;

                /* . Y. */
                Real2DArray_Item ( otransformationY,  PX, PSIGMA   ) = caY*sb + ca*sbY ;
                Real2DArray_Item ( otransformationY,  PX, PPIPLUS  ) = caY*cb + ca*cbY ;
                Real2DArray_Item ( otransformationY,  PX, PPIMINUS ) = -saY ;
                Real2DArray_Item ( otransformationY,  PY, PSIGMA   ) = saY*sb + sa*sbY ;
                Real2DArray_Item ( otransformationY,  PY, PPIPLUS  ) = saY*cb + sa*cbY ;
                Real2DArray_Item ( otransformationY,  PY, PPIMINUS ) = caY  ;
                Real2DArray_Item ( otransformationY,  PZ, PSIGMA   ) = cbY  ;
                Real2DArray_Item ( otransformationY,  PZ, PPIPLUS  ) = -sbY ;

                /* . Z. */
                Real2DArray_Item ( otransformationZ,  PX, PSIGMA   ) = ca*sbZ ;
                Real2DArray_Item ( otransformationZ,  PX, PPIPLUS  ) = ca*cbZ ;
                Real2DArray_Item ( otransformationZ,  PY, PSIGMA   ) = sa*sbZ ;
                Real2DArray_Item ( otransformationZ,  PY, PPIPLUS  ) = sa*cbZ ;
                Real2DArray_Item ( otransformationZ,  PZ, PSIGMA   ) = cbZ  ;
                Real2DArray_Item ( otransformationZ,  PZ, PPIPLUS  ) = -sbZ ;
            }

            /* . Swap transformation. */
            if ( axesSwapped )
            {
                Real2DArray_Item ( stransformation,  S,  S ) = 1.0e+00 ;
                Real2DArray_Item ( stransformation, PX, PZ ) = 1.0e+00 ;
                Real2DArray_Item ( stransformation, PY, PY ) = 1.0e+00 ;
                Real2DArray_Item ( stransformation, PZ, PX ) = 1.0e+00 ;
            }

# ifdef MNDODORBITALS
            /* . d-orbital transformation. */
            if ( norbitals == 9 )
            {
                auto Real pt5sq3 ;

                /* . Initialization. */
                pt5sq3 = 0.5e+00 * sqrt ( 3.0e+00 ) ;
                c2a = 2.0e+00 * ca * ca - 1.0e+00 ;
                c2b = 2.0e+00 * cb * cb - 1.0e+00 ;
                s2a = 2.0e+00 * sa * ca ;
                s2b = 2.0e+00 * sb * cb ;

                /* . Rotation matrix. */
                Real2DArray_Item ( otransformation, DX2Y2, DSIGMA      ) = pt5sq3 * c2a * sb * sb ;
                Real2DArray_Item ( otransformation, DX2Y2, DPIPLUS     ) = 0.5e+00 * c2a * s2b ;
                Real2DArray_Item ( otransformation, DX2Y2, DPIMINUS    ) = -s2a * sb ;
                Real2DArray_Item ( otransformation, DX2Y2, DDELTAPLUS  ) = c2a * ( cb * cb + 0.5e+00 * sb * sb ) ;
                Real2DArray_Item ( otransformation, DX2Y2, DDELTAMINUS ) = -s2a * cb ;
                Real2DArray_Item ( otransformation, DXZ  , DSIGMA      ) = pt5sq3 * ca * s2b ;
                Real2DArray_Item ( otransformation, DXZ  , DPIPLUS     ) = ca * c2b ;
                Real2DArray_Item ( otransformation, DXZ  , DPIMINUS    ) = -sa * cb ;
                Real2DArray_Item ( otransformation, DXZ  , DDELTAPLUS  ) = -0.5e+00 * ca * s2b ;
                Real2DArray_Item ( otransformation, DXZ  , DDELTAMINUS ) = sa * sb ;
                Real2DArray_Item ( otransformation, DZ2  , DSIGMA      ) = cb * cb - 0.5e+00 * sb * sb ;
                Real2DArray_Item ( otransformation, DZ2  , DPIPLUS     ) = -pt5sq3 * s2b ;
/*                Real2DArray_Item ( otransformation, DZ2  , DPIMINUS    ) = 0.0e+00 ; */
                Real2DArray_Item ( otransformation, DZ2  , DDELTAPLUS  ) = pt5sq3 * sb * sb ;
/*                Real2DArray_Item ( otransformation, DZ2  , DDELTAMINUS ) = 0.0e+00 ; */
                Real2DArray_Item ( otransformation, DYZ  , DSIGMA      ) = pt5sq3 * sa * s2b ;
                Real2DArray_Item ( otransformation, DYZ  , DPIPLUS     ) = sa * c2b ;
                Real2DArray_Item ( otransformation, DYZ  , DPIMINUS    ) = ca * cb ;
                Real2DArray_Item ( otransformation, DYZ  , DDELTAPLUS  ) = -0.5e+00 * sa * s2b ;
                Real2DArray_Item ( otransformation, DYZ  , DDELTAMINUS ) = -ca * sb ;
                Real2DArray_Item ( otransformation, DXY  , DSIGMA      ) = pt5sq3 * s2a * sb * sb ;
                Real2DArray_Item ( otransformation, DXY  , DPIPLUS     ) = 0.5e+00 * s2a * s2b ;
                Real2DArray_Item ( otransformation, DXY  , DPIMINUS    ) = c2a * sb ;
                Real2DArray_Item ( otransformation, DXY  , DDELTAPLUS  ) = s2a * ( cb * cb + 0.5e+00 * sb * sb ) ;
                Real2DArray_Item ( otransformation, DXY  , DDELTAMINUS ) = c2a * cb ;

                /* . Derivatives. */
                if ( doGradients )
                {
                    auto Real c2aX, c2aY, c2bX, c2bY, c2bZ, s2aX, s2aY, s2bX, s2bY, s2bZ ;

                    /* . Various factors. */
                    c2aX = 4.0e+00 * ca * caX ;
                    c2aY = 4.0e+00 * ca * caY ;
                    s2aX = 2.0e+00 * ( saX * ca + sa * caX ) ;
                    s2aY = 2.0e+00 * ( saY * ca + sa * caY ) ;
                    c2bX = 4.0e+00 * cb * cbX ;
                    c2bY = 4.0e+00 * cb * cbY ;
                    c2bZ = 4.0e+00 * cb * cbZ ;
                    s2bX = 2.0e+00 * ( sbX * cb + sb * cbX ) ;
                    s2bY = 2.0e+00 * ( sbY * cb + sb * cbY ) ;
                    s2bZ = 2.0e+00 * ( sbZ * cb + sb * cbZ ) ;

                    /* . X. */
                    Real2DArray_Item ( otransformationX, DX2Y2, DSIGMA      ) = pt5sq3 * ( c2aX * sb * sb + 2.0e+00 * c2a * sb * sbX ) ;
                    Real2DArray_Item ( otransformationX, DX2Y2, DPIPLUS     ) = 0.5e+00 * ( c2aX * s2b + c2a * s2bX ) ;
                    Real2DArray_Item ( otransformationX, DX2Y2, DPIMINUS    ) = - s2aX * sb - s2a * sbX ;
                    Real2DArray_Item ( otransformationX, DX2Y2, DDELTAPLUS  ) = c2aX * ( cb * cb + 0.5e+00 * sb * sb ) + c2a * ( 2.0e+00 * cb * cbX + sb * sbX ) ;
                    Real2DArray_Item ( otransformationX, DX2Y2, DDELTAMINUS ) = - s2aX * cb - s2a * cbX ;
                    Real2DArray_Item ( otransformationX, DXZ  , DSIGMA      ) = pt5sq3 * ( caX * s2b + ca * s2bX ) ;
                    Real2DArray_Item ( otransformationX, DXZ  , DPIPLUS     ) = caX * c2b + ca * c2bX ;
                    Real2DArray_Item ( otransformationX, DXZ  , DPIMINUS    ) = - saX * cb - sa * cbX ;
                    Real2DArray_Item ( otransformationX, DXZ  , DDELTAPLUS  ) = - 0.5e+00 * ( caX * s2b + ca * s2bX ) ;
                    Real2DArray_Item ( otransformationX, DXZ  , DDELTAMINUS ) = saX * sb + sa * sbX ;
                    Real2DArray_Item ( otransformationX, DZ2  , DSIGMA      ) = 2.0e+00 * cb * cbX - sb * sbX ;
                    Real2DArray_Item ( otransformationX, DZ2  , DPIPLUS     ) = - pt5sq3 * s2bX ;
                    Real2DArray_Item ( otransformationX, DZ2  , DDELTAPLUS  ) = pt5sq3 * 2.0e+00 * sb * sbX ;
                    Real2DArray_Item ( otransformationX, DYZ  , DSIGMA      ) = pt5sq3 * ( saX * s2b + sa * s2bX ) ;
                    Real2DArray_Item ( otransformationX, DYZ  , DPIPLUS     ) = saX * c2b + sa * c2bX ;
                    Real2DArray_Item ( otransformationX, DYZ  , DPIMINUS    ) = caX * cb + ca * cbX ;
                    Real2DArray_Item ( otransformationX, DYZ  , DDELTAPLUS  ) = - 0.5e+00 * ( saX * s2b + sa * s2bX ) ;
                    Real2DArray_Item ( otransformationX, DYZ  , DDELTAMINUS ) = - caX * sb - ca * sbX ;
                    Real2DArray_Item ( otransformationX, DXY  , DSIGMA      ) = pt5sq3 * ( s2aX * sb * sb + 2.0e+00 * s2a * sb * sbX ) ;
                    Real2DArray_Item ( otransformationX, DXY  , DPIPLUS     ) = 0.5e+00 * ( s2aX * s2b + s2a * s2bX ) ;
                    Real2DArray_Item ( otransformationX, DXY  , DPIMINUS    ) = c2aX * sb + c2a * sbX ;
                    Real2DArray_Item ( otransformationX, DXY  , DDELTAPLUS  ) = s2aX * ( cb * cb + 0.5e+00 * sb * sb ) + s2a * ( 2.0e+00 * cb * cbX + sb * sbX ) ;
                    Real2DArray_Item ( otransformationX, DXY  , DDELTAMINUS ) = c2aX * cb + c2a * cbX ;

                    /* . Y. */
                    Real2DArray_Item ( otransformationY, DX2Y2, DSIGMA      ) = pt5sq3 * ( c2aY * sb * sb + 2.0e+00 * c2a * sb * sbY ) ;
                    Real2DArray_Item ( otransformationY, DX2Y2, DPIPLUS     ) = 0.5e+00 * ( c2aY * s2b + c2a * s2bY ) ;
                    Real2DArray_Item ( otransformationY, DX2Y2, DPIMINUS    ) = - s2aY * sb - s2a * sbY ;
                    Real2DArray_Item ( otransformationY, DX2Y2, DDELTAPLUS  ) = c2aY * ( cb * cb + 0.5e+00 * sb * sb ) + c2a * ( 2.0e+00 * cb * cbY + sb * sbY ) ;
                    Real2DArray_Item ( otransformationY, DX2Y2, DDELTAMINUS ) = - s2aY * cb - s2a * cbY ;
                    Real2DArray_Item ( otransformationY, DXZ  , DSIGMA      ) = pt5sq3 * ( caY * s2b + ca * s2bY ) ;
                    Real2DArray_Item ( otransformationY, DXZ  , DPIPLUS     ) = caY * c2b + ca * c2bY ;
                    Real2DArray_Item ( otransformationY, DXZ  , DPIMINUS    ) = - saY * cb - sa * cbY ;
                    Real2DArray_Item ( otransformationY, DXZ  , DDELTAPLUS  ) = - 0.5e+00 * ( caY * s2b + ca * s2bY ) ;
                    Real2DArray_Item ( otransformationY, DXZ  , DDELTAMINUS ) = saY * sb + sa * sbY ;
                    Real2DArray_Item ( otransformationY, DZ2  , DSIGMA      ) = 2.0e+00 * cb * cbY - sb * sbY ;
                    Real2DArray_Item ( otransformationY, DZ2  , DPIPLUS     ) = - pt5sq3 * s2bY ;
                    Real2DArray_Item ( otransformationY, DZ2  , DDELTAPLUS  ) = pt5sq3 * 2.0e+00 * sb * sbY ;
                    Real2DArray_Item ( otransformationY, DYZ  , DSIGMA      ) = pt5sq3 * ( saY * s2b + sa * s2bY ) ;
                    Real2DArray_Item ( otransformationY, DYZ  , DPIPLUS     ) = saY * c2b + sa * c2bY ;
                    Real2DArray_Item ( otransformationY, DYZ  , DPIMINUS    ) = caY * cb + ca * cbY ;
                    Real2DArray_Item ( otransformationY, DYZ  , DDELTAPLUS  ) = - 0.5e+00 * ( saY * s2b + sa * s2bY ) ;
                    Real2DArray_Item ( otransformationY, DYZ  , DDELTAMINUS ) = - caY * sb - ca * sbY ;
                    Real2DArray_Item ( otransformationY, DXY  , DSIGMA      ) = pt5sq3 * ( s2aY * sb * sb + 2.0e+00 * s2a * sb * sbY ) ;
                    Real2DArray_Item ( otransformationY, DXY  , DPIPLUS     ) = 0.5e+00 * ( s2aY * s2b + s2a * s2bY ) ;
                    Real2DArray_Item ( otransformationY, DXY  , DPIMINUS    ) = c2aY * sb + c2a * sbY ;
                    Real2DArray_Item ( otransformationY, DXY  , DDELTAPLUS  ) = s2aY * ( cb * cb + 0.5e+00 * sb * sb ) + s2a * ( 2.0e+00 * cb * cbY + sb * sbY ) ;
                    Real2DArray_Item ( otransformationY, DXY  , DDELTAMINUS ) = c2aY * cb + c2a * cbY ;

                    /* . Z. */
                    Real2DArray_Item ( otransformationZ, DX2Y2, DSIGMA      ) = pt5sq3 * 2.0e+00 * c2a * sb * sbZ ;
                    Real2DArray_Item ( otransformationZ, DX2Y2, DPIPLUS     ) = 0.5e+00 * c2a * s2bZ ;
                    Real2DArray_Item ( otransformationZ, DX2Y2, DPIMINUS    ) = - s2a * sbZ ;
                    Real2DArray_Item ( otransformationZ, DX2Y2, DDELTAPLUS  ) = c2a * ( 2.0e+00 * cb * cbZ + sb * sbZ ) ;
                    Real2DArray_Item ( otransformationZ, DX2Y2, DDELTAMINUS ) = - s2a * cbZ ;
                    Real2DArray_Item ( otransformationZ, DXZ  , DSIGMA      ) = pt5sq3 * ca * s2bZ ;
                    Real2DArray_Item ( otransformationZ, DXZ  , DPIPLUS     ) = ca * c2bZ ;
                    Real2DArray_Item ( otransformationZ, DXZ  , DPIMINUS    ) = - sa * cbZ ;
                    Real2DArray_Item ( otransformationZ, DXZ  , DDELTAPLUS  ) = - 0.5e+00 * ca * s2bZ ;
                    Real2DArray_Item ( otransformationZ, DXZ  , DDELTAMINUS ) = sa * sbZ ;
                    Real2DArray_Item ( otransformationZ, DZ2  , DSIGMA      ) = 2.0e+00 * cb * cbZ - sb * sbZ ;
                    Real2DArray_Item ( otransformationZ, DZ2  , DPIPLUS     ) = - pt5sq3 * s2bZ ;
                    Real2DArray_Item ( otransformationZ, DZ2  , DDELTAPLUS  ) = pt5sq3 * 2.0e+00 * sb * sbZ ;
                    Real2DArray_Item ( otransformationZ, DYZ  , DSIGMA      ) = pt5sq3 * sa * s2bZ ;
                    Real2DArray_Item ( otransformationZ, DYZ  , DPIPLUS     ) = sa * c2bZ ;
                    Real2DArray_Item ( otransformationZ, DYZ  , DPIMINUS    ) = ca * cbZ ;
                    Real2DArray_Item ( otransformationZ, DYZ  , DDELTAPLUS  ) = - 0.5e+00 * sa * s2bZ ;
                    Real2DArray_Item ( otransformationZ, DYZ  , DDELTAMINUS ) = - ca * sbZ ;
                    Real2DArray_Item ( otransformationZ, DXY  , DSIGMA      ) = pt5sq3 * 2.0e+00 * s2a * sb * sbZ ;
                    Real2DArray_Item ( otransformationZ, DXY  , DPIPLUS     ) = 0.5e+00 * s2a * s2bZ ;
                    Real2DArray_Item ( otransformationZ, DXY  , DPIMINUS    ) = c2a * sbZ ;
                    Real2DArray_Item ( otransformationZ, DXY  , DDELTAPLUS  ) = s2a * ( 2.0e+00 * cb * cbZ + sb * sbZ ) ;
                    Real2DArray_Item ( otransformationZ, DXY  , DDELTAMINUS ) = c2a * cbZ ;
                }

                /* . Swap transformation. */
                if ( axesSwapped )
                {
                    Real2DArray_Item ( stransformation, DXZ,   DXZ   ) =  1.0e+00 ;
                    Real2DArray_Item ( stransformation, DXY,   DYZ   ) =  1.0e+00 ;
                    Real2DArray_Item ( stransformation, DYZ,   DXY   ) =  1.0e+00 ;
                    Real2DArray_Item ( stransformation, DX2Y2, DX2Y2 ) =  0.5e+00 ;
                    Real2DArray_Item ( stransformation, DX2Y2, DZ2   ) =  pt5sq3  ;
                    Real2DArray_Item ( stransformation, DZ2,   DX2Y2 ) =  pt5sq3  ;
                    Real2DArray_Item ( stransformation, DZ2,   DZ2   ) = -0.5e+00 ;
                }
            }
# endif /*MNDODORBITALS*/

            /* . Swap transformation. */
            if ( axesSwapped )
            {
                auto Real2DArray *new = NULL, *swap = NULL ;

                /* . Premultiply the existing transformation by the swap transformation. */
                new  = Real2DArray_Allocate ( norbitals, norbitals, NULL ) ;
                Real2DArray_MatrixMultiply ( False, False, 1.0e+00, stransformation, otransformation, 0.0e+00, new, NULL ) ;
                swap = otransformation ; otransformation = new ; new = swap ;

                /* . Same for gradients but also swap X and Z matrices. */
                if ( doGradients )
                {
                    Real2DArray_MatrixMultiply ( False, False, 1.0e+00, stransformation, otransformationX, 0.0e+00, new, NULL ) ; swap = otransformationX ; otransformationX = new ; new = swap ;
                    Real2DArray_MatrixMultiply ( False, False, 1.0e+00, stransformation, otransformationY, 0.0e+00, new, NULL ) ; swap = otransformationY ; otransformationY = new ; new = swap ;
                    Real2DArray_MatrixMultiply ( False, False, 1.0e+00, stransformation, otransformationZ, 0.0e+00, new, NULL ) ; swap = otransformationZ ; otransformationZ = new ; new = swap ;
                    swap = otransformationZ ; otransformationZ = otransformationX ; otransformationX = swap ;
                }

                /* . Clear up. */
                Real2DArray_Deallocate ( &new ) ;
            }

            /* . Allocate the largest transformation matrix. */
            n = ( norbitals * ( norbitals + 1 ) ) / 2 ;
            ctransformation = Real2DArray_Allocate ( n, n, NULL ) ;
            Real2DArray_Set ( ctransformation, 0.0e+00 ) ;

            /* . Form the matrix. */
            for ( i = ij = 0 ; i < norbitals ; i++ )
            {
                for ( j = 0 ; j <= i ; j++, ij++ )
                {
                    for ( m = mn = 0 ; m < norbitals ; m++, mn++ )
                    {
                        for ( n = 0 ; n < m ; n++, mn++ ) Real2DArray_Item ( ctransformation, ij, mn ) = Real2DArray_Item ( otransformation, i, m ) * Real2DArray_Item ( otransformation, j, n ) +
                                                                                                         Real2DArray_Item ( otransformation, i, n ) * Real2DArray_Item ( otransformation, j, m ) ;
                        Real2DArray_Item ( ctransformation, ij, mn ) = Real2DArray_Item ( otransformation, i, m ) * Real2DArray_Item ( otransformation, j, m ) ;
                    }
                }
            }

            /* . Derivatives. */
            if ( doGradients )
            {
                /* . Allocation. */
                n = ( norbitals * ( norbitals + 1 ) ) / 2 ;
                ctransformationX = Real2DArray_Allocate ( n, n, NULL ) ; Real2DArray_Set ( ctransformationX, 0.0e+00 ) ;
                ctransformationY = Real2DArray_Allocate ( n, n, NULL ) ; Real2DArray_Set ( ctransformationY, 0.0e+00 ) ;
                ctransformationZ = Real2DArray_Allocate ( n, n, NULL ) ; Real2DArray_Set ( ctransformationZ, 0.0e+00 ) ;

                /* . Form the matrices. */
                for ( i = ij = 0 ; i < norbitals ; i++ )
                {
                    for ( j = 0 ; j <= i ; j++, ij++ )
                    {
                        for ( m = mn = 0 ; m < norbitals ; m++, mn++ )
                        {
                            for ( n = 0 ; n < m ; n++, mn++ )
                            {
                                Real2DArray_Item ( ctransformationX, ij, mn ) = Real2DArray_Item ( otransformationX, i, m ) * Real2DArray_Item ( otransformation,  j, n ) +
                                                                                Real2DArray_Item ( otransformation,  i, m ) * Real2DArray_Item ( otransformationX, j, n ) +
                                                                                Real2DArray_Item ( otransformationX, i, n ) * Real2DArray_Item ( otransformation,  j, m ) +
                                                                                Real2DArray_Item ( otransformation,  i, n ) * Real2DArray_Item ( otransformationX, j, m ) ;
                                Real2DArray_Item ( ctransformationY, ij, mn ) = Real2DArray_Item ( otransformationY, i, m ) * Real2DArray_Item ( otransformation,  j, n ) +
                                                                                Real2DArray_Item ( otransformation,  i, m ) * Real2DArray_Item ( otransformationY, j, n ) +
                                                                                Real2DArray_Item ( otransformationY, i, n ) * Real2DArray_Item ( otransformation,  j, m ) +
                                                                                Real2DArray_Item ( otransformation,  i, n ) * Real2DArray_Item ( otransformationY, j, m ) ;
                                Real2DArray_Item ( ctransformationZ, ij, mn ) = Real2DArray_Item ( otransformationZ, i, m ) * Real2DArray_Item ( otransformation,  j, n ) +
                                                                                Real2DArray_Item ( otransformation,  i, m ) * Real2DArray_Item ( otransformationZ, j, n ) +
                                                                                Real2DArray_Item ( otransformationZ, i, n ) * Real2DArray_Item ( otransformation,  j, m ) +
                                                                                Real2DArray_Item ( otransformation,  i, n ) * Real2DArray_Item ( otransformationZ, j, m ) ;
                            }
                            Real2DArray_Item ( ctransformationX, ij, mn ) = Real2DArray_Item ( otransformationX, i, m ) * Real2DArray_Item ( otransformation,  j, m ) +
                                                                            Real2DArray_Item ( otransformation,  i, m ) * Real2DArray_Item ( otransformationX, j, m ) ;
                            Real2DArray_Item ( ctransformationY, ij, mn ) = Real2DArray_Item ( otransformationY, i, m ) * Real2DArray_Item ( otransformation,  j, m ) +
                                                                            Real2DArray_Item ( otransformation,  i, m ) * Real2DArray_Item ( otransformationY, j, m ) ;
                            Real2DArray_Item ( ctransformationZ, ij, mn ) = Real2DArray_Item ( otransformationZ, i, m ) * Real2DArray_Item ( otransformation,  j, m ) +
                                                                            Real2DArray_Item ( otransformation,  i, m ) * Real2DArray_Item ( otransformationZ, j, m ) ;
                        }
                    }
                }
            }

            /* . Define the smaller matrix as a sub-block of the larger one. */
            /* . Should use slices or views here. */
            m = Minimum ( ni, nj ) ;
            if ( ( m > 1 ) && ( m < norbitals ) )
            {
                n = ( m * ( m + 1 ) ) / 2 ;
                clesser = Real2DArray_Allocate ( n, n, NULL ) ;
                for ( i = 0 ; i < n ; i++ )
                {
                    for ( j = 0 ; j < n ; j++ ) Real2DArray_Item ( clesser, i, j ) = Real2DArray_Item ( ctransformation, i, j ) ;
                }
                if ( doGradients )
                {
                    clesserX = Real2DArray_Allocate ( n, n, NULL ) ;
                    clesserY = Real2DArray_Allocate ( n, n, NULL ) ;
                    clesserZ = Real2DArray_Allocate ( n, n, NULL ) ;
                    for ( i = 0 ; i < n ; i++ )
                    {
                        for ( j = 0 ; j < n ; j++ )
                        {
                            Real2DArray_Item ( clesserX, i, j ) = Real2DArray_Item ( ctransformationX, i, j ) ;
                            Real2DArray_Item ( clesserY, i, j ) = Real2DArray_Item ( ctransformationY, i, j ) ;
                            Real2DArray_Item ( clesserZ, i, j ) = Real2DArray_Item ( ctransformationZ, i, j ) ;
                        }
                    }
                }
            }

            /* . Assign the matrices to i and j.*/
            if ( ni == norbitals ) it = ctransformation ;
            else                   it = clesser         ;
            if ( nj == norbitals ) jt = ctransformation ;
            else                   jt = clesser         ;
            if ( doGradients )
            {
                if ( ni == norbitals ) { itX = ctransformationX ; itY = ctransformationY ; itZ = ctransformationZ ; }
                else                   { itX = clesserX         ; itY = clesserY         ; itZ = clesserZ         ; }
                if ( nj == norbitals ) { jtX = ctransformationX ; jtY = ctransformationY ; jtZ = ctransformationZ ; }
                else                   { jtX = clesserX         ; jtY = clesserY         ; jtZ = clesserZ         ; }
            }

            /* . Clear up. */
            Real2DArray_Deallocate ( &otransformation ) ;
            Real2DArray_Deallocate ( &stransformation ) ;
            if ( doGradients )
            {
                Real2DArray_Deallocate ( &otransformationX ) ;
                Real2DArray_Deallocate ( &otransformationY ) ;
                Real2DArray_Deallocate ( &otransformationZ ) ;
            }
        }

        /* . Finish up. */
        (*itransformation) = it ;
        (*jtransformation) = jt ;
        if ( doGradients )
        {
            (*itransformationX) = itX ;
            (*itransformationY) = itY ;
            (*itransformationZ) = itZ ;
            (*jtransformationX) = jtX ;
            (*jtransformationY) = jtY ;
            (*jtransformationZ) = jtZ ;
        }
    }
    return r ;
}
# undef ALIGNMENT_TOLERANCE

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate all two-center TEIs and optionally their derivatives in the local frame.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MNDOIntegralUtilities_LocalFrame2CTEIs ( const MNDOParameters *idata, const MNDOParameters *jdata, const Real r, Real2DArray  *lfteis, Real1DArray  *core1b, Real1DArray  *core2a,
                                                                                                                                 Real2DArray *dlfteis, Real1DArray *dcore1b, Real1DArray *dcore2a )
{
    Boolean doGradients ;
# ifdef MNDODORBITALS
    const Integer *negative = NULL, *positive = NULL, *unique = NULL ;
          Integer c, i, iam = 0, ij, j, jam = 0, k, kl, l, nnegative, npositive, nunique, t ;
# else
    const Integer *positive = NULL ;
          Integer c, iam = 0, jam = 0, npositive, t ;
# endif

    /* . Check for gradients. */
    doGradients = ( dlfteis != NULL ) && ( dcore1b != NULL ) && ( dcore2a != NULL ) ;

    /* . Initialization. */
    Real2DArray_Set ( lfteis, 0.0e+00 ) ;
    Real1DArray_Set ( core1b, 0.0e+00 ) ;
    Real1DArray_Set ( core2a, 0.0e+00 ) ;
    if ( doGradients )
    {
        Real2DArray_Set ( dlfteis, 0.0e+00 ) ;
        Real1DArray_Set ( dcore1b, 0.0e+00 ) ;
        Real1DArray_Set ( dcore2a, 0.0e+00 ) ;
    }

    /* . Get the highest AM for each atom. */
    switch ( idata->norbitals )
    {
        case 1: iam = 0 ; break ;
        case 4: iam = 1 ; break ;
# ifdef MNDODORBITALS
        case 9: iam = 2 ; break ;
# endif
    }
    switch ( jdata->norbitals )
    {
        case 1: jam = 0 ; break ;
        case 4: jam = 1 ; break ;
# ifdef MNDODORBITALS
        case 9: jam = 2 ; break ;
# endif
    }

    /* . Unique non-zero integrals in the local frame and then those related by symmetry. */
    /* . sp integrals. */
    MNDOIntegralUtilities_LocalFrame2CTEIsSP ( idata, jdata, r,        lfteis, dlfteis ) ;
    MNDOIntegralUtilities_LocalFrame2COEIsSP ( idata, jdata, r, False, core1b, dcore1b ) ;
    MNDOIntegralUtilities_LocalFrame2COEIsSP ( jdata, idata, r, True,  core2a, dcore2a ) ;

# ifdef MNDODORBITALS
    /* . Integrals involving d-orbitals. */
    /* . TEIs. */
    /* . Define the sets of integrals to evaluate. */
    nnegative = NNEGATIVE[iam][jam] ;
    npositive = NPOSITIVE[iam][jam] ;
    nunique   = NUNIQUE  [iam][jam] ;
    switch ( jam )
    {
        case 0:
            negative = SNEGATIVE ;
            positive = SPOSITIVE ;
            unique   = SUNIQUE   ;
            break ;
        case 1:
            negative = SPNEGATIVE ;
            positive = SPPOSITIVE ;
            unique   = SPUNIQUE   ;
            break ;
        case 2:
            negative = SPDNEGATIVE ;
            positive = SPDPOSITIVE ;
            unique   = SPDUNIQUE   ;
            break ;
    }

    /* . Integrals. */
    for ( c = t = 0 ; c < nunique ; c++, t += 4 )
    {
        i  = unique[t  ] ;
        j  = unique[t+1] ;
        k  = unique[t+2] ;
        l  = unique[t+3] ;
        ij = ( i * ( i + 1 ) ) / 2 + j ;
        kl = ( k * ( k + 1 ) ) / 2 + l ;
        Real2DArray_Item ( lfteis, ij, kl ) = MNDOIntegralUtilities_LocalFrame2CTEI ( &MNDOIntegralUtilities_2CChargeInteraction, idata, jdata, ij, kl, ORBITALAM[i], ORBITALAM[j], ORBITALAM[k], ORBITALAM[l], 0, r ) ;
    }
    for ( c = t = 0 ; c < npositive ; c++, t+= 4 ) Real2DArray_Item ( lfteis, positive[t], positive[t+1] ) =   Real2DArray_Item ( lfteis, positive[t+2], positive[t+3] ) ;
    for ( c = t = 0 ; c < nnegative ; c++, t+= 4 ) Real2DArray_Item ( lfteis, negative[t], negative[t+1] ) = - Real2DArray_Item ( lfteis, negative[t+2], negative[t+3] ) ;

    /* . First atom electron - second atom core terms. */
    for ( c = t = 0 ; c < NCUNIQUE[iam] ; c++, t += 2 )
    {
        i  = CUNIQUE[t  ] ;
        j  = CUNIQUE[t+1] ;
        ij = ( i * ( i + 1 ) ) / 2 + j ;
        Real1DArray_Item ( core1b, ij ) = MNDOIntegralUtilities_LocalFrame2CTEI ( &MNDOIntegralUtilities_2CChargeInteraction, idata, jdata, ij, SS, ORBITALAM[i], ORBITALAM[j], 0, 0, 2, r ) ;
    }
    for ( c = t = 0 ; c < NCPOSITIVE[iam] ; c++, t+= 2 ) Real1DArray_Item ( core1b, CPOSITIVE[t] ) =   Real1DArray_Item ( core1b, CPOSITIVE[t+1] ) ;
    for ( c = t = 0 ; c < NCNEGATIVE[iam] ; c++, t+= 2 ) Real1DArray_Item ( core1b, CNEGATIVE[t] ) = - Real1DArray_Item ( core1b, CNEGATIVE[t+1] ) ;
    Real1DArray_Scale ( core1b, - jdata->zcore ) ;

    /* . Second atom electron - first atom core terms. */
    for ( c = t = 0 ; c < NCUNIQUE[jam] ; c++, t += 2 )
    {
        k  = CUNIQUE[t  ] ;
        l  = CUNIQUE[t+1] ;
        kl = ( k * ( k + 1 ) ) / 2 + l ;
        Real1DArray_Item ( core2a, kl ) = MNDOIntegralUtilities_LocalFrame2CTEI ( &MNDOIntegralUtilities_2CChargeInteraction, idata, jdata, SS, kl, 0, 0, ORBITALAM[k], ORBITALAM[l], 1, r ) ;
    }
    for ( c = t = 0 ; c < NCPOSITIVE[jam] ; c++, t+= 2 ) Real1DArray_Item ( core2a, CPOSITIVE[t] ) =   Real1DArray_Item ( core2a, CPOSITIVE[t+1] ) ;
    for ( c = t = 0 ; c < NCNEGATIVE[jam] ; c++, t+= 2 ) Real1DArray_Item ( core2a, CNEGATIVE[t] ) = - Real1DArray_Item ( core2a, CNEGATIVE[t+1] ) ;
    Real1DArray_Scale ( core2a, - idata->zcore ) ;

    /* . Derivatives. */
    if ( doGradients )
    {
        for ( c = t = 0 ; c < nunique ; c++, t += 4 )
        {
            i  = unique[t  ] ;
            j  = unique[t+1] ;
            k  = unique[t+2] ;
            l  = unique[t+3] ;
            ij = ( i * ( i + 1 ) ) / 2 + j ;
            kl = ( k * ( k + 1 ) ) / 2 + l ;
            Real2DArray_Item ( dlfteis, ij, kl ) = MNDOIntegralUtilities_LocalFrame2CTEI ( &MNDOIntegralUtilities_2CChargeInteractionD, idata, jdata, ij, kl, ORBITALAM[i], ORBITALAM[j], ORBITALAM[k], ORBITALAM[l], 0, r ) ;
        }
        for ( c = t = 0 ; c < npositive ; c++, t+= 4 ) Real2DArray_Item ( dlfteis, positive[t], positive[t+1] ) =   Real2DArray_Item ( dlfteis, positive[t+2], positive[t+3] ) ;
        for ( c = t = 0 ; c < nnegative ; c++, t+= 4 ) Real2DArray_Item ( dlfteis, negative[t], negative[t+1] ) = - Real2DArray_Item ( dlfteis, negative[t+2], negative[t+3] ) ;

        /* . First atom electron - second atom core terms. */
        for ( c = t = 0 ; c < NCUNIQUE[iam] ; c++, t += 2 )
        {
            i  = CUNIQUE[t  ] ;
            j  = CUNIQUE[t+1] ;
            ij = ( i * ( i + 1 ) ) / 2 + j ;
            Real1DArray_Item ( dcore1b, ij ) = MNDOIntegralUtilities_LocalFrame2CTEI ( &MNDOIntegralUtilities_2CChargeInteractionD, idata, jdata, ij, SS, ORBITALAM[i], ORBITALAM[j], 0, 0, 2, r ) ;
        }
        for ( c = t = 0 ; c < NCPOSITIVE[iam] ; c++, t+= 2 ) Real1DArray_Item ( dcore1b, CPOSITIVE[t] ) =   Real1DArray_Item ( dcore1b, CPOSITIVE[t+1] ) ;
        for ( c = t = 0 ; c < NCNEGATIVE[iam] ; c++, t+= 2 ) Real1DArray_Item ( dcore1b, CNEGATIVE[t] ) = - Real1DArray_Item ( dcore1b, CNEGATIVE[t+1] ) ;
        Real1DArray_Scale ( dcore1b, - jdata->zcore ) ;

        /* . Second atom electron - first atom core terms. */
        for ( c = t = 0 ; c < NCUNIQUE[jam] ; c++, t += 2 )
        {
            k  = CUNIQUE[t  ] ;
            l  = CUNIQUE[t+1] ;
            kl = ( k * ( k + 1 ) ) / 2 + l ;
            Real1DArray_Item ( dcore2a, kl ) = MNDOIntegralUtilities_LocalFrame2CTEI ( &MNDOIntegralUtilities_2CChargeInteractionD, idata, jdata, SS, kl, 0, 0, ORBITALAM[k], ORBITALAM[l], 1, r ) ;
        }
        for ( c = t = 0 ; c < NCPOSITIVE[jam] ; c++, t+= 2 ) Real1DArray_Item ( dcore2a, CPOSITIVE[t] ) =   Real1DArray_Item ( dcore2a, CPOSITIVE[t+1] ) ;
        for ( c = t = 0 ; c < NCNEGATIVE[jam] ; c++, t+= 2 ) Real1DArray_Item ( dcore2a, CNEGATIVE[t] ) = - Real1DArray_Item ( dcore2a, CNEGATIVE[t+1] ) ;
        Real1DArray_Scale ( dcore2a, - idata->zcore ) ;
    }
# else
    /* . TEIs. */
    /* . Define the sets of integrals to evaluate. */
    npositive = NPOSITIVE[iam][jam] ;
    switch ( jam )
    {
        case 0:
            positive = SPOSITIVE ;
            break ;
        case 1:
            positive = SPPOSITIVE ;
            break ;
    }

    /* . Integrals. */
    for ( c = t = 0 ; c < npositive ; c++, t+= 4 ) Real2DArray_Item ( lfteis, positive[t], positive[t+1] ) =   Real2DArray_Item ( lfteis, positive[t+2], positive[t+3] ) ;

    /* . First atom electron - second atom core terms. */
    for ( c = t = 0 ; c < NCPOSITIVE[iam] ; c++, t+= 2 ) Real1DArray_Item ( core1b, CPOSITIVE[t] ) =   Real1DArray_Item ( core1b, CPOSITIVE[t+1] ) ;
    Real1DArray_Scale ( core1b, - jdata->zcore ) ;

    /* . Second atom electron - first atom core terms. */
    for ( c = t = 0 ; c < NCPOSITIVE[jam] ; c++, t+= 2 ) Real1DArray_Item ( core2a, CPOSITIVE[t] ) =   Real1DArray_Item ( core2a, CPOSITIVE[t+1] ) ;
    Real1DArray_Scale ( core2a, - idata->zcore ) ;

    /* . Derivatives. */
    if ( doGradients )
    {
        for ( c = t = 0 ; c < npositive ; c++, t+= 4 ) Real2DArray_Item ( dlfteis, positive[t], positive[t+1] ) =   Real2DArray_Item ( dlfteis, positive[t+2], positive[t+3] ) ;

        /* . First atom electron - second atom core terms. */
        for ( c = t = 0 ; c < NCPOSITIVE[iam] ; c++, t+= 2 ) Real1DArray_Item ( dcore1b, CPOSITIVE[t] ) =   Real1DArray_Item ( dcore1b, CPOSITIVE[t+1] ) ;
        Real1DArray_Scale ( dcore1b, - jdata->zcore ) ;

        /* . Second atom electron - first atom core terms. */
        for ( c = t = 0 ; c < NCPOSITIVE[jam] ; c++, t+= 2 ) Real1DArray_Item ( dcore2a, CPOSITIVE[t] ) =   Real1DArray_Item ( dcore2a, CPOSITIVE[t+1] ) ;
        Real1DArray_Scale ( dcore2a, - idata->zcore ) ;
    }
# endif
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the unique OEIs and TEIs involving sp orbitals in the local frame.
! . Derivatives (wrt r) can also be optionally calculated.
! . Unfortunately this code is incompatible with the d-orbital code (as the formulae are different).
! . This is identical to REPP with a few sign changes to cope with the change in transformation (1,4,7,8,9,12,13,14).
! . The OEIs and TEIs are treated separately.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define EV1 0.5e+00
# define EV2 0.25e+00
# define EV3 0.125e+00
# define EV4 0.0625e+00
# define PP  0.5e+00
# define TD  2.0e+00
void MNDOIntegralUtilities_LocalFrame2COEIsSP ( const MNDOParameters *idata, const MNDOParameters *jdata, const Real r, const Boolean swapped, Real1DArray *core, Real1DArray *dcore )
{
    Boolean    doGradients ;
    Real  da, qa ;
    Real  ade, aee, aqe, dze, gdze, gqxxe, gqzze, gri0, gri2, gri3, qxxe, qzze, ri0, ri2, ri3, rsq, sqr1, sqr2, sqr3, sqr4, sqr5, sqr6, xxx ;

    /* . idata is the electronic atom, jdata is the core. */

    /* . Check for gradients. */
    doGradients = ( dcore != NULL ) ;

    /* . s/s - always done. */
    aee = pow ( ( idata->po[0] + jdata->po[8] ), 2 ) ;
    rsq = r * r ;
    ri0 = 1.0e+00 / sqrt ( rsq + aee ) ;
    Real1DArray_Item ( core, SS ) = ri0 ;
    if ( doGradients )
    {
        gri0 = - r * ri0 * ri0 * ri0 ;
        Real1DArray_Item ( dcore, SS ) = gri0 ;
    }

    /* . sp/s. */
    /* . Redo ri0 with po[6]. */
    if ( idata->norbitals > 3 )
    {
        da   = idata->dd ;
        qa   = idata->qq * TD ;
        if ( swapped ) { da *= -1.0e+00 ; qa *= -1.0e+00 ; }
        aee  = pow ( ( idata->po[6] + jdata->po[8] ), 2 ) ;
        ade  = pow ( ( idata->po[1] + jdata->po[8] ), 2 ) ;
        aqe  = pow ( ( idata->po[2] + jdata->po[8] ), 2 ) ;
        xxx  = r+da      ; sqr1 = 1.0e+00 / sqrt ( xxx*xxx + ade ) ;
        xxx  = r-da      ; sqr2 = 1.0e+00 / sqrt ( xxx*xxx + ade ) ;
        xxx  = r+qa      ; sqr3 = 1.0e+00 / sqrt ( xxx*xxx + aqe ) ;
        xxx  = r-qa      ; sqr4 = 1.0e+00 / sqrt ( xxx*xxx + aqe ) ;
        xxx  = rsq + aqe ; sqr5 = 1.0e+00 / sqrt ( xxx ) ;
                           sqr6 = 1.0e+00 / sqrt ( xxx + qa*qa ) ;
        dze  = - EV1 * ( sqr1 - sqr2 ) ;
        qzze = EV2 * ( sqr3 + sqr4 ) - EV1 * sqr5 ;
        qxxe = EV1 * ( sqr6 - sqr5 ) ;
        ri0  = 1.0e+00 / sqrt ( rsq + aee ) ;
        ri2  = ri0 + qzze ;
        ri3  = ri0 + qxxe ;
        Real1DArray_Item ( core, PZS  ) = dze ;
        Real1DArray_Item ( core, PZPZ ) = ri2 ;
        Real1DArray_Item ( core, PXPX ) = ri3 ;
        if ( doGradients )
        {
            xxx   = r * sqr5*sqr5*sqr5 ;
            gdze  = EV1 * ( (r+da) * sqr1*sqr1*sqr1 - (r-da) * sqr2*sqr2*sqr2 ) ;
            gqzze = - EV2 * ( (r+qa) * sqr3*sqr3*sqr3 + (r-qa) * sqr4*sqr4*sqr4 ) + EV1 * xxx ;
            gqxxe = - EV1 * ( r * sqr6*sqr6*sqr6 - xxx ) ;
            gri0  = - r * ri0 * ri0 * ri0 ;
            gri2  = gri0 + gqzze ;
            gri3  = gri0 + gqxxe ;
            Real1DArray_Item ( dcore, PZS  ) = gdze ;
            Real1DArray_Item ( dcore, PZPZ ) = gri2 ;
            Real1DArray_Item ( dcore, PXPX ) = gri3 ;
        }
    }
}

void MNDOIntegralUtilities_LocalFrame2CTEIsSP ( const MNDOParameters *idata, const MNDOParameters *jdata, const Real r, Real2DArray *lfteis, Real2DArray *dlfteis )
{
    Boolean    doGradients ;
    Real  da = 0.0e+00, db = 0.0e+00, qa = 0.0e+00, qa0, qb = 0.0e+00, qb0 ;
    Real  ade, aee, aed, adq, aqe, aeq, aqd, aqq, arg35, arg38, arg39, axx, rsq, www, xxx, yyy, zzz,
            dze = 0.0e+00, edz = 0.0e+00, dxdx, dzdz, qxxe = 0.0e+00, eqxx = 0.0e+00, qzze = 0.0e+00, eqzz = 0.0e+00, dzqxx, qxxdz, dxqxz, qxzdx, dzqzz, qzzdz, qxxqxx, qxxqyy, qxxqzz, qzzqxx, qxzqxz, qzzqzz,
            gdze = 0.0e+00, gedz = 0.0e+00, gdxdx, gdzdz, gqxxe = 0.0e+00, geqxx = 0.0e+00, gqzze = 0.0e+00, geqzz = 0.0e+00, gdzqxx, gqxxdz, gdxqxz, gqxzdx, gdzqzz, gqzzdz, gqxxqxx, gqxxqyy, gqxxqzz, gqzzqxx, gqxzqxz, gqzzqzz,
            gri0 = 0.0e+00,  gri2, gri3, gri10, gri11, ri0, ri2, ri3, ri10, ri11,
            sqr1, sqr2, sqr3, sqr4, sqr5, sqr6, sqr7, sqr8, sqr9, sqr10, sqr11, sqr12, sqr13, sqr14, sqr15, sqr16, sqr17, sqr18, sqr19, sqr20,
            sqr21, sqr22, sqr23, sqr24, sqr25, sqr26, sqr27, sqr28, sqr29, sqr30, sqr31, sqr32, sqr33, sqr34, sqr35, sqr36, sqr37, sqr38, sqr39, sqr40,
            sqr41, sqr42, sqr43, sqr44, sqr45, sqr46, sqr47, sqr48, sqr49, sqr50, sqr51, sqr52, sqr53, sqr54, sqr55, sqr56, sqr57, sqr58, sqr59, sqr60,
/*            sqr61, sqr62, sqr63, */
            sqr64, sqr65, sqr66, sqr67, sqr68, sqr69, sqr70, sqr71 ;

    /* . Check for gradients. */
    doGradients = ( dlfteis != NULL ) ;

    /* . s/s - always done. */
    aee = pow ( ( idata->po[0] + jdata->po[0] ), 2 ) ;
    rsq = r * r ;
    ri0 = 1.0e+00 / sqrt ( rsq + aee ) ;
    Real2DArray_Item ( lfteis, SS, SS ) = ri0 ;
    if ( doGradients )
    {
        gri0 = - r * ri0 * ri0 * ri0 ;
        Real2DArray_Item ( dlfteis, SS, SS ) = gri0 ;
    }

    /* . sp/s. */
    if ( idata->norbitals > 3 )
    {
        da   = idata->dd ;
        qa   = idata->qq * TD ;
        ade  = pow ( ( idata->po[1] + jdata->po[0] ), 2 ) ;
        aqe  = pow ( ( idata->po[2] + jdata->po[0] ), 2 ) ;
        xxx  = r+da      ; sqr1 = 1.0e+00 / sqrt ( xxx*xxx + ade ) ;
        xxx  = r-da      ; sqr2 = 1.0e+00 / sqrt ( xxx*xxx + ade ) ;
        xxx  = r+qa      ; sqr3 = 1.0e+00 / sqrt ( xxx*xxx + aqe ) ;
        xxx  = r-qa      ; sqr4 = 1.0e+00 / sqrt ( xxx*xxx + aqe ) ;
        xxx  = rsq + aqe ; sqr5 = 1.0e+00 / sqrt ( xxx ) ;
                           sqr6 = 1.0e+00 / sqrt ( xxx + qa*qa ) ;
        dze  = - EV1 * ( sqr1 - sqr2 ) ;
        qzze = EV2 * ( sqr3 + sqr4 ) - EV1 * sqr5 ;
        qxxe = EV1 * ( sqr6 - sqr5 ) ;
        ri2  = ri0 + qzze ;
        ri3  = ri0 + qxxe ;
        Real2DArray_Item ( lfteis, PZS , SS ) = dze ;
        Real2DArray_Item ( lfteis, PZPZ, SS ) = ri2 ;
        Real2DArray_Item ( lfteis, PXPX, SS ) = ri3 ;
        if ( doGradients )
        {
            xxx   = r * sqr5*sqr5*sqr5 ;
            gdze  = EV1 * ( (r+da) * sqr1*sqr1*sqr1 - (r-da) * sqr2*sqr2*sqr2 ) ;
            gqzze = - EV2 * ( (r+qa) * sqr3*sqr3*sqr3 + (r-qa) * sqr4*sqr4*sqr4 ) + EV1 * xxx ;
            gqxxe = - EV1 * ( r * sqr6*sqr6*sqr6 - xxx ) ;
            gri2  = gri0 + gqzze ;
            gri3  = gri0 + gqxxe ;
            Real2DArray_Item ( dlfteis, PZS , SS ) = gdze ;
            Real2DArray_Item ( dlfteis, PZPZ, SS ) = gri2 ;
            Real2DArray_Item ( dlfteis, PXPX, SS ) = gri3 ;
        }
    }

    /* . s/sp. */
    if ( jdata->norbitals > 3 )
    {
        db    = jdata->dd ;
        qb    = jdata->qq * TD ;
        aed   = pow ( ( idata->po[0] + jdata->po[1] ), 2 ) ;
        aeq   = pow ( ( idata->po[0] + jdata->po[2] ), 2 ) ;
        xxx   = r-db      ; sqr7  = 1.0e+00 / sqrt ( xxx*xxx + aed ) ;
        xxx   = r+db      ; sqr8  = 1.0e+00 / sqrt ( xxx*xxx + aed ) ;
        xxx   = r-qb      ; sqr9  = 1.0e+00 / sqrt ( xxx*xxx + aeq ) ;
        xxx   = r+qb      ; sqr10 = 1.0e+00 / sqrt ( xxx*xxx + aeq ) ;
        xxx   = rsq + aeq ; sqr11 = 1.0e+00 / sqrt ( xxx ) ;
                            sqr12 = 1.0e+00 / sqrt ( xxx + qb*qb ) ;
        edz   = - EV1 * ( sqr7 - sqr8 ) ;
        eqzz  = EV2 * ( sqr9 + sqr10 ) - EV1 * sqr11 ;
        eqxx  = EV1 * ( sqr12 - sqr11 ) ;
        ri10  = ri0 + eqzz ;
        ri11  = ri0 + eqxx ;
        Real2DArray_Item ( lfteis, SS, PZS  ) = edz  ;
        Real2DArray_Item ( lfteis, SS, PZPZ ) = ri10 ;
        Real2DArray_Item ( lfteis, SS, PXPX ) = ri11 ;
        if ( doGradients )
        {
            xxx   = r * sqr11*sqr11*sqr11 ;
            gedz  = EV1 * ( (r-db) * sqr7*sqr7*sqr7 - (r+db) * sqr8*sqr8*sqr8 ) ;
            geqzz = - EV2 * ( (r-qb) * sqr9*sqr9*sqr9 + (r+qb) * sqr10*sqr10*sqr10 ) + EV1 * xxx ;
            geqxx = - EV1 * ( r * sqr12*sqr12*sqr12 - xxx ) ;
            gri10 = gri0 + geqzz ;
            gri11 = gri0 + geqxx ;
            Real2DArray_Item ( dlfteis, SS, PZS  ) = gedz  ;
            Real2DArray_Item ( dlfteis, SS, PZPZ ) = gri10 ;
            Real2DArray_Item ( dlfteis, SS, PXPX ) = gri11 ;
        }
    }

    /* . sp/sp. */
    if ( ( idata->norbitals > 3 ) && ( jdata->norbitals > 3 ) )
    {
        axx   = pow ( ( idata->po[1] + jdata->po[1] ), 2 ) ;
        adq   = pow ( ( idata->po[1] + jdata->po[2] ), 2 ) ;
        aqd   = pow ( ( idata->po[2] + jdata->po[1] ), 2 ) ;
        aqq   = pow ( ( idata->po[2] + jdata->po[2] ), 2 ) ;
        xxx   = da - db         ; sqr13 = 1.0e+00 / sqrt ( rsq + axx + xxx * xxx ) ;
        xxx   = da + db         ; sqr14 = 1.0e+00 / sqrt ( rsq + axx + xxx * xxx ) ;
        xxx   = r + da - db     ; sqr15 = 1.0e+00 / sqrt ( xxx * xxx + axx ) ;
        xxx   = r - da + db     ; sqr16 = 1.0e+00 / sqrt ( xxx * xxx + axx ) ;
        xxx   = r - da - db     ; sqr17 = 1.0e+00 / sqrt ( xxx * xxx + axx ) ;
        xxx   = r + da + db     ; sqr18 = 1.0e+00 / sqrt ( xxx * xxx + axx ) ;
        xxx   = r + da          ;
        yyy   = xxx * xxx + adq ; sqr19 = 1.0e+00 / sqrt ( yyy ) ;
                                  sqr20 = 1.0e+00 / sqrt ( yyy + qb * qb ) ;
        xxx   = r - da          ;
        yyy   = xxx * xxx + adq ; sqr21 = 1.0e+00 / sqrt ( yyy ) ;
                                  sqr22 = 1.0e+00 / sqrt ( yyy + qb * qb ) ;
        xxx   = r - db          ;
        yyy   = xxx * xxx + aqd ; sqr23 = 1.0e+00 / sqrt ( yyy ) ;
                                  sqr24 = 1.0e+00 / sqrt ( yyy + qa * qa ) ;
        xxx   = r + db          ;
        yyy   = xxx * xxx + aqd ; sqr25 = 1.0e+00 / sqrt ( yyy ) ;
                                  sqr26 = 1.0e+00 / sqrt ( yyy + qa * qa ) ;
        xxx   = r + da - qb     ; sqr27 = 1.0e+00 / sqrt ( xxx * xxx + adq ) ;
        xxx   = r - da - qb     ; sqr28 = 1.0e+00 / sqrt ( xxx * xxx + adq ) ;
        xxx   = r + da + qb     ; sqr29 = 1.0e+00 / sqrt ( xxx * xxx + adq ) ;
        xxx   = r - da + qb     ; sqr30 = 1.0e+00 / sqrt ( xxx * xxx + adq ) ;
        xxx   = r + qa - db     ; sqr31 = 1.0e+00 / sqrt ( xxx * xxx + aqd ) ;
        xxx   = r + qa + db     ; sqr32 = 1.0e+00 / sqrt ( xxx * xxx + aqd ) ;
        xxx   = r - qa - db     ; sqr33 = 1.0e+00 / sqrt ( xxx * xxx + aqd ) ;
        xxx   = r - qa + db     ; sqr34 = 1.0e+00 / sqrt ( xxx * xxx + aqd ) ;
        arg35 = rsq + aqq       ; sqr35 = 1.0e+00 / sqrt ( arg35 ) ;
        xxx   = qa - qb         ; sqr36 = 1.0e+00 / sqrt ( arg35 + xxx * xxx ) ;
        xxx   = qa + qb         ; sqr37 = 1.0e+00 / sqrt ( arg35 + xxx * xxx ) ;
        xxx   = arg35 + qa * qa ; sqr38 = 1.0e+00 / sqrt ( xxx ) ;
                                  sqr39 = 1.0e+00 / sqrt ( arg35 + qb * qb ) ;
                                  sqr40 = 1.0e+00 / sqrt ( xxx + qb * qb ) ;
        xxx   = r - qb          ;
        yyy   = xxx * xxx + aqq ; sqr41 = 1.0e+00 / sqrt ( yyy ) ;
                                  sqr42 = 1.0e+00 / sqrt ( yyy + qa * qa ) ;
        xxx   = r + qb          ;
        yyy   = xxx * xxx + aqq ; sqr43 = 1.0e+00 / sqrt ( yyy ) ;
                                  sqr44 = 1.0e+00 / sqrt ( yyy + qa * qa ) ;
        xxx   = r + qa          ;
        yyy   = xxx * xxx + aqq ; sqr45 = 1.0e+00 / sqrt ( yyy ) ;
                                  sqr46 = 1.0e+00 / sqrt ( yyy + qb * qb ) ;
        xxx   = r - qa          ;
        yyy   = xxx * xxx + aqq ; sqr47 = 1.0e+00 / sqrt ( yyy ) ;
                                  sqr48 = 1.0e+00 / sqrt ( yyy + qb * qb ) ;
        xxx   = r + qa - qb     ; sqr49 = 1.0e+00 / sqrt ( xxx * xxx + aqq ) ;
        xxx   = r + qa + qb     ; sqr50 = 1.0e+00 / sqrt ( xxx * xxx + aqq ) ;
        xxx   = r - qa - qb     ; sqr51 = 1.0e+00 / sqrt ( xxx * xxx + aqq ) ;
        xxx   = r - qa + qb     ; sqr52 = 1.0e+00 / sqrt ( xxx * xxx + aqq ) ;
        qa0   = idata->qq ;
        qb0   = jdata->qq ;
        xxx   = pow ( ( da - qb0 ), 2 ) ;
        yyy   = pow ( ( r  - qb0 ), 2 ) ;
        zzz   = pow ( ( da + qb0 ), 2 ) ;
        www   = pow ( ( r + qb0  ), 2 ) ;
        sqr53 = 1.0e+00 / sqrt ( xxx + yyy + adq ) ;
        sqr54 = 1.0e+00 / sqrt ( xxx + www + adq ) ;
        sqr55 = 1.0e+00 / sqrt ( zzz + yyy + adq ) ;
        sqr56 = 1.0e+00 / sqrt ( zzz + www + adq ) ;
        xxx   = pow ( ( qa0 - db ), 2 ) ;
        yyy   = pow ( ( qa0 + db ), 2 ) ;
        zzz   = pow ( ( r + qa0  ), 2 ) ;
        www   = pow ( ( r - qa0  ), 2 ) ;
        sqr57 = 1.0e+00 / sqrt ( zzz + xxx + aqd ) ;
        sqr58 = 1.0e+00 / sqrt ( www + xxx + aqd ) ;
        sqr59 = 1.0e+00 / sqrt ( zzz + yyy + aqd ) ;
        sqr60 = 1.0e+00 / sqrt ( www + yyy + aqd ) ;
        xxx   = pow ( ( qa0 - qb0 ), 2 ) ;
        yyy   = pow ( ( qa0 + qb0 ), 2 ) ;
/*
        sqr61 = 1.0e+00 / sqrt ( arg35 + TD * xxx ) ;
        sqr62 = 1.0e+00 / sqrt ( arg35 + TD * yyy ) ;
        sqr63 = 1.0e+00 / sqrt ( arg35 + TD * (qa0 * qa0 + qb0 * qb0) ) ;
*/
        zzz   = pow ( ( r + qa0 - qb0 ), 2 ) ;
        sqr64 = 1.0e+00 / sqrt ( zzz + xxx + aqq ) ;
        sqr65 = 1.0e+00 / sqrt ( zzz + yyy + aqq ) ;
        zzz   = pow ( ( r + qa0 + qb0 ), 2 ) ;
        sqr66 = 1.0e+00 / sqrt ( zzz + xxx + aqq ) ;
        sqr67 = 1.0e+00 / sqrt ( zzz + yyy + aqq ) ;
        zzz   = pow ( ( r - qa0 - qb0 ), 2 ) ;
        sqr68 = 1.0e+00 / sqrt ( zzz + xxx + aqq ) ;
        sqr69 = 1.0e+00 / sqrt ( zzz + yyy + aqq ) ;
        zzz   = pow ( ( r - qa0 + qb0 ), 2 ) ;
        sqr70 = 1.0e+00 / sqrt ( zzz + xxx + aqq ) ;
        sqr71 = 1.0e+00 / sqrt ( zzz + yyy + aqq ) ;
        dxdx   =   EV1 * sqr13 - EV1 * sqr14 ;
        dzdz   =   EV2 * sqr15 + EV2 * sqr16 - EV2 * sqr17 - EV2 * sqr18 ;
        dzqxx  =   EV2 * sqr19 - EV2 * sqr20 - EV2 * sqr21 + EV2 * sqr22 ;
        qxxdz  =   EV2 * sqr23 - EV2 * sqr24 - EV2 * sqr25 + EV2 * sqr26 ;
        dzqzz  = - EV3 * sqr27 + EV3 * sqr28 - EV3 * sqr29 + EV3 * sqr30 - EV2 * sqr21 + EV2 * sqr19 ;
        qzzdz  = - EV3 * sqr31 + EV3 * sqr32 - EV3 * sqr33 + EV3 * sqr34 + EV2 * sqr23 - EV2 * sqr25 ;
        qxxqxx =   EV3 * sqr36 + EV3 * sqr37 - EV2 * sqr38 - EV2 * sqr39 + EV2 * sqr35 ;
        qxxqyy =   EV2 * sqr40 - EV2 * sqr38 - EV2 * sqr39 + EV2 * sqr35 ;
        qxxqzz =   EV3 * sqr42 + EV3 * sqr44 - EV3 * sqr41 - EV3 * sqr43 - EV2 * sqr38 + EV2 * sqr35 ;
        qzzqxx =   EV3 * sqr46 + EV3 * sqr48 - EV3 * sqr45 - EV3 * sqr47 - EV2 * sqr39 + EV2 * sqr35 ;
        qzzqzz =   EV4 * sqr49 + EV4 * sqr50 + EV4 * sqr51 + EV4 * sqr52 - EV3 * sqr47 - EV3 * sqr45 - EV3 * sqr41 - EV3 * sqr43 + EV2 * sqr35 ;
        dxqxz  = - EV2 * sqr53 + EV2 * sqr54 + EV2 * sqr55 - EV2 * sqr56 ;
        qxzdx  = - EV2 * sqr57 + EV2 * sqr58 + EV2 * sqr59 - EV2 * sqr60 ;
        qxzqxz =   EV3 * sqr64 - EV3 * sqr66 - EV3 * sqr68 + EV3 * sqr70 - EV3 * sqr65 + EV3 * sqr67 + EV3 * sqr69 - EV3 * sqr71 ;
        Real2DArray_Item ( lfteis, PZS , PZS  ) = dzdz ;
        Real2DArray_Item ( lfteis, PXS , PXS  ) = dxdx ;
        Real2DArray_Item ( lfteis, PZPZ, PZS  ) = edz + qzzdz ;
        Real2DArray_Item ( lfteis, PXPX, PZS  ) = edz + qxxdz ;
        Real2DArray_Item ( lfteis, PXPZ, PXS  ) = qxzdx ;
        Real2DArray_Item ( lfteis, PZS , PZPZ ) = dze + dzqzz ;
        Real2DArray_Item ( lfteis, PZS , PXPX ) = dze + dzqxx ;
        Real2DArray_Item ( lfteis, PXS , PXPZ ) = dxqxz ;
        Real2DArray_Item ( lfteis, PZPZ, PZPZ ) = ri0 + eqzz + qzze + qzzqzz ;
        Real2DArray_Item ( lfteis, PXPX, PZPZ ) = ri0 + eqzz + qxxe + qxxqzz ;
        Real2DArray_Item ( lfteis, PZPZ, PXPX ) = ri0 + eqxx + qzze + qzzqxx ;
        Real2DArray_Item ( lfteis, PXPX, PXPX ) = ri0 + eqxx + qxxe + qxxqxx ;
        Real2DArray_Item ( lfteis, PXPZ, PXPZ ) = qxzqxz ;
        Real2DArray_Item ( lfteis, PXPX, PYPY ) = ri0 + eqxx + qxxe + qxxqyy ;
        Real2DArray_Item ( lfteis, PYPX, PYPX ) = PP * ( qxxqxx - qxxqyy ) ;
        if ( doGradients )
        {
            gdxdx   = - EV1 * r * ( sqr13*sqr13*sqr13 - sqr14*sqr14*sqr14 ) ;
            gdzdz   = - EV2 * ( (r+da-db)*sqr15*sqr15*sqr15 + (r-da+db)*sqr16*sqr16*sqr16 - (r-da-db)*sqr17*sqr17*sqr17 - (r+da+db)*sqr18*sqr18*sqr18 ) ;
            www     = (r+da)*sqr19*sqr19*sqr19 ;
            xxx     = (r-da)*sqr21*sqr21*sqr21 ;
            gdzqxx  = - EV2 * ( www - (r+da)*sqr20*sqr20*sqr20 - xxx + (r-da)*sqr22*sqr22*sqr22 ) ;
            yyy     = (r-db)*sqr23*sqr23*sqr23 ;
            zzz     = (r+db)*sqr25*sqr25*sqr25 ;
            gqxxdz  = - EV2 * ( yyy - (r-db)*sqr24*sqr24*sqr24 - zzz + (r+db)*sqr26*sqr26*sqr26 ) ;
            gdzqzz  = - EV3 * ( - (r+da-qb)*sqr27*sqr27*sqr27 + (r-da-qb)*sqr28*sqr28*sqr28 - (r+da+qb)*sqr29*sqr29*sqr29 + (r-da+qb)*sqr30*sqr30*sqr30 ) + EV2 * ( xxx - www ) ;
            gqzzdz  = - EV3 * ( - (r+qa-db)*sqr31*sqr31*sqr31 + (r+qa+db)*sqr32*sqr32*sqr32 - (r-qa-db)*sqr33*sqr33*sqr33 + (r-qa+db)*sqr34*sqr34*sqr34 ) - EV2 * ( yyy - zzz ) ;
            arg35   = r*sqr35*sqr35*sqr35 ;
            arg38   = r*sqr38*sqr38*sqr38 ;
            arg39   = r*sqr39*sqr39*sqr39 ;
            gqxxqxx = - EV3 * ( r*sqr36*sqr36*sqr36 + r*sqr37*sqr37*sqr37 ) + EV2 * ( arg38 + arg39 - arg35 ) ;
            gqxxqyy = - EV2 * ( r*sqr40*sqr40*sqr40 - arg38 - arg39 + arg35 ) ;
            www     = (r-qb)*sqr41*sqr41*sqr41 ;
            xxx     = (r+qb)*sqr43*sqr43*sqr43 ;
            yyy     = (r+qa)*sqr45*sqr45*sqr45 ;
            zzz     = (r-qa)*sqr47*sqr47*sqr47 ;
            gqxxqzz = - EV3 * ( (r-qb)*sqr42*sqr42*sqr42 + (r+qb)*sqr44*sqr44*sqr44 - www - xxx ) + EV2 * ( arg38 - arg35 ) ;
            gqzzqxx = - EV3 * ( (r+qa)*sqr46*sqr46*sqr46 + (r-qa)*sqr48*sqr48*sqr48 - yyy - zzz ) + EV2 * ( arg39 - arg35 ) ;
            gqzzqzz = - EV4 * ( (r+qa-qb)*sqr49*sqr49*sqr49 + (r+qa+qb)*sqr50*sqr50*sqr50 + (r-qa-qb)*sqr51*sqr51*sqr51 + (r-qa+qb)*sqr52*sqr52*sqr52 ) + EV3 * ( zzz + yyy + www + xxx ) - EV2 * arg35 ;
            gdxqxz  = - EV2 * ( - (r-qb0)*sqr53*sqr53*sqr53 + (r+qb0)*sqr54*sqr54*sqr54 + (r-qb0)*sqr55*sqr55*sqr55 - (r+qb0)*sqr56*sqr56*sqr56 ) ;
            gqxzdx  = - EV2 * ( - (r+qa0)*sqr57*sqr57*sqr57 + (r-qa0)*sqr58*sqr58*sqr58 + (r+qa0)*sqr59*sqr59*sqr59 - (r-qa0)*sqr60*sqr60*sqr60 ) ;
            gqxzqxz = - EV3 * ( (r+qa0-qb0)*sqr64*sqr64*sqr64 - (r+qa0+qb0)*sqr66*sqr66*sqr66 - (r-qa0-qb0)*sqr68*sqr68*sqr68 + (r-qa0+qb0)*sqr70*sqr70*sqr70 -
                                (r+qa0-qb0)*sqr65*sqr65*sqr65 + (r+qa0+qb0)*sqr67*sqr67*sqr67 + (r-qa0-qb0)*sqr69*sqr69*sqr69 - (r-qa0+qb0)*sqr71*sqr71*sqr71 ) ;
            Real2DArray_Item ( dlfteis, PZS , PZS  ) = gdzdz ;
            Real2DArray_Item ( dlfteis, PXS , PXS  ) = gdxdx ;
            Real2DArray_Item ( dlfteis, PZPZ, PZS  ) = gedz + gqzzdz ;
            Real2DArray_Item ( dlfteis, PXPX, PZS  ) = gedz + gqxxdz ;
            Real2DArray_Item ( dlfteis, PXPZ, PXS  ) = gqxzdx ;
            Real2DArray_Item ( dlfteis, PZS , PZPZ ) = gdze + gdzqzz ;
            Real2DArray_Item ( dlfteis, PZS , PXPX ) = gdze + gdzqxx ;
            Real2DArray_Item ( dlfteis, PXS , PXPZ ) = gdxqxz ;
            Real2DArray_Item ( dlfteis, PZPZ, PZPZ ) = gri0 + geqzz + gqzze + gqzzqzz ;
            Real2DArray_Item ( dlfteis, PXPX, PZPZ ) = gri0 + geqzz + gqxxe + gqxxqzz ;
            Real2DArray_Item ( dlfteis, PZPZ, PXPX ) = gri0 + geqxx + gqzze + gqzzqxx ;
            Real2DArray_Item ( dlfteis, PXPX, PXPX ) = gri0 + geqxx + gqxxe + gqxxqxx ;
            Real2DArray_Item ( dlfteis, PXPZ, PXPZ ) = gqxzqxz ;
            Real2DArray_Item ( dlfteis, PXPX, PYPY ) = gri0 + geqxx + gqxxe + gqxxqyy ;
            Real2DArray_Item ( dlfteis, PYPX, PYPX ) = PP * ( gqxxqxx - gqxxqyy ) ;
        }
    }
}
# undef EV1
# undef EV2
# undef EV3
# undef EV4
# undef PP
# undef TD

# ifdef MNDODORBITALS
/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate a two-center TEI or its derivative in the local frame.
!---------------------------------------------------------------------------------------------------------------------------------*/
/*
! . i, j, k, l - either 0, 1, or 2 with i >= j and k >= l.
! . c - either 0, 1 or 2.
*/
Real MNDOIntegralUtilities_LocalFrame2CTEI ( const ChargeInteractionFunction Evaluate, const MNDOParameters *idata, const MNDOParameters *jdata, const Integer ij, const Integer kl,
                                                                                                      const Integer i, const Integer j, const Integer k, const Integer l, const Integer c, const Real r )
{
    Real integral = 0.0e+00 ;
    if ( ( idata != NULL ) && ( jdata != NULL ) && ( NCHTERMS[ij] > 0 ) && ( NCHTERMS[kl] > 0 ) )
    {
        auto Real add, chijkl[MAXCHTERMS], dij = 0.0e+00, dkl = 0.0e+00, pij = 0.0e+00, pkl = 0.0e+00 ;
        auto Integer    lij, lkl, lmin, lm1, lm2, l1, l1max, l1min, l1offset, l1terms[MAXCHTERMS], l2, l2max, l2min, l2offset, l2terms[MAXCHTERMS], m, mterms[MAXCHTERMS], nterms, t ;
        /* . Get loop indices. */
        /* . Possibilites for ( i, j, l1min, l1max ) are ( 0, 0, 0, 0 ), ( 1, 0, 1, 1 ), ( 1, 1, 0, 2 ), ( 2, 0, 2, 2 ), ( 2, 1, 1, 2 ) and ( 2, 2, 0, 2 ). */
        /* . Special cases (l1 or l2 == 0) are the diagonal ones; i.e. i == j with i = 0, 1 or 2. */
        l1min = i - j ;
        l1max = Minimum ( i + j, 2 ) ;
        lij   = ( i * ( i + 1 ) ) / 2 + j ;
        l2min = k - l ;
        l2max = Minimum ( k + l, 2 ) ;
        lkl   = ( k * ( k + 1 ) ) / 2 + l ;
        /* . Preprocessing for number of terms. */
        for ( l1 = l1min, nterms = 0 ; l1 <= l1max ; l1++ )
        {
            l1offset = ij * CHINCREMENT1 + l1 * CHINCREMENT2 + CHINCREMENT3 ;
            for ( l2 = l2min ; l2 <= l2max ; l2++ )
            {
                l2offset = kl * CHINCREMENT1 + l2 * CHINCREMENT2 + CHINCREMENT3 ;
                lmin     = Minimum ( l1, l2 ) ;
                for ( m = -lmin ; m <= lmin ; m++ )
                {
                    lm1 = CHINDICES[l1offset+m] ;
                    lm2 = CHINDICES[l2offset+m] ;
                    if ( ( lm1 > 0 ) && ( lm2 > 0 ) )
                    {
                        l1terms[nterms] = l1 ;
                        l2terms[nterms] = l2 ;
                        mterms [nterms] = abs ( m ) ;
                        chijkl [nterms] = CHTERMS[lm1-1] * CHTERMS[lm2-1] ;
                        nterms ++ ;
                    }
                }
            }
        }
        /* . Calculate the terms. */
        for ( t = 0, integral = 0.0e+00 ; t < nterms ; t++ )
        {
            l1 = l1terms[t] ;
            l2 = l2terms[t] ;
            m  = mterms [t] ;
            if ( l1 == 0 )
            {
                dij = 0.0e+00 ;
                switch ( i )
                {
                    case 0:
                        pij = idata->po[0] ;
                        if ( c == 1 ) pij = idata->po[8] ;
                        break ;
                    case 1:
                        pij = idata->po[6] ;
                        break ;
                    case 2:
                        pij = idata->po[7] ;
                        break ;
                }
            }
            else
            {
                dij = idata->ddp[lij] ;
                pij = idata->po [lij] ;
            }
            if ( l2 == 0 )
            {
                dkl = 0.0e+00 ;
                switch ( k )
                {
                    case 0:
                        pkl = jdata->po[0] ;
                        if ( c == 2 ) pkl = jdata->po[8] ;
                        break ;
                    case 1:
                        pkl = jdata->po[6] ;
                        break ;
                    case 2:
                        pkl = jdata->po[7] ;
                        break ;
                }
            }
            else
            {
                dkl = jdata->ddp[lkl] ;
                pkl = jdata->po [lkl] ;
            }
            add = pow ( ( pij + pkl ), 2 ) ;
            integral += chijkl[t] * Evaluate ( r, l1, l2, m, dij, dkl, add ) ;
        }
    }
    return integral ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the interaction of two point-charge configurations.
!---------------------------------------------------------------------------------------------------------------------------------*/
/*
! . r      distance.
! . l1,m   quantum numbers for multipole of configuration 1.
! . l2,m   quantum numbers for multipole of configuration 2.
! . da     charge separation of configuration 1.
! . db     charge separation of configuration 2.
! . add    additive term.
*/
static Real MNDOIntegralUtilities_2CChargeInteraction ( const Real r, const Integer l1, const Integer l2, const Integer m, const Real da, const Real db, const Real add )
{
    Real aa, ab, charg, dxdx, dxqxz, dzdz, dzqzz, qqzz, qxzdx, qxzqxz, qzzdz, qzzq, xyxy, zzzz ;
    charg = 0.0e+00 ;
    switch ( l1 )
    {
        case 0:
            switch ( l2 )
            {
                case 0:
                    charg = 1.0e+00 / sqrt( r*r + add ) ;
                break ;
                case 1:
                    charg = 1.0e+00 / sqrt(pow(r + db,2) + add)  - 1.0e+00/sqrt(pow(r - db,2) + add) ; charg /= 2.0e+00 ;
                break ;
                case 2:
                    qqzz  = 1.0e+00/sqrt(pow(r - db,2) + add) - 2.0e+00/sqrt(r*r + db*db + add) + 1.0e+00/sqrt(pow(r + db,2) + add) ;
                    charg = qqzz/4.0e+00 ;
                break ;
            }
            break ;
        case 1:
            switch ( l2 )
            {
                case 0:
                    charg = (-1.0e+00/sqrt(pow(r + da,2) + add)) + 1.0e+00/sqrt(pow(r - da,2) + add) ; charg /= 2.0e+00 ;
                break ;
                case 1:
                if ( m == 0 )
                {
                    dzdz  = 1.0e+00/sqrt(pow(r + da - db,2) + add) + 1.0e+00/sqrt(pow(r - da + db,2) + add) - 1.0e+00/sqrt(pow(r - da - db,2) + add) - 1.0e+00/sqrt(pow(r + da + db,2) + add) ;
                    charg = dzdz/4.0e+00 ;
                }
                else if ( m == 1 )
                {
                    dxdx  = 2.0e+00/sqrt(r*r + pow(da - db,2) + add) - 2.0e+00/sqrt(r*r + pow(da + db,2) + add) ;
                    charg = dxdx/4.0e+00 ;
                }
                break ;
                case 2:
                if ( m == 0 )
                {
                    dzqzz = 1.0e+00/sqrt(pow(r - da - db,2) + add) - 2.0e+00/sqrt(pow(r - da,2) + db*db + add) + 1.0e+00/sqrt(pow(r + db - da,2) + add) -
                            1.0e+00/sqrt(pow(r - db + da,2) + add) + 2.0e+00/sqrt(pow(r + da,2) + db*db + add) - 1.0e+00/sqrt(pow(r + da + db,2) + add) ;
                    charg = dzqzz/8.0e+00 ;
                }
                else if ( m == 1 )
                {
                    ab = db/sqrt(2.0e+00) ;
                    dxqxz = (-2.0e+00/sqrt(pow(r - ab,2) + pow(da - ab,2) + add)) + 2.0e+00/sqrt(pow(r + ab,2) + pow(da - ab,2) + add) + 2.0e+00/sqrt(pow(r - ab,2) + pow(da + ab,2) + add) - 2.0e+00/sqrt(pow(r + ab,2) + pow(da + ab,2) + add) ;
                    charg = dxqxz/8.0e+00 ;
                }
                break ;
            }
            break ;
        case 2:
            switch ( l2 )
            {
                case 0:
                    qzzq  = 1.0e+00/sqrt(pow(r - da,2) + add) - 2.0e+00/sqrt(r*r + da*da + add) + 1.0e+00/sqrt(pow(r + da,2) + add) ;
                    charg = qzzq/4.0e+00 ;
                break ;
                case 1:
                if ( m == 0 )
                {
                    qzzdz = (-1.0e+00/sqrt(pow(r - da - db,2) + add)) + 2.0e+00/sqrt(pow(r - db,2) + da*da + add) - 1.0e+00/sqrt(pow(r + da - db,2) + add) +
                              1.0e+00/sqrt(pow(r - da + db,2) + add) - 2.0e+00/sqrt(pow(r + db,2) + da*da + add) + 1.0e+00/sqrt(pow(r + da + db,2) + add) ;
                    charg = qzzdz/8.0e+00 ;
                }
                else if ( m == 1 )
                {
                    aa = da/sqrt(2.0e+00) ;
                    qxzdx = (-2.0e+00/sqrt(pow(r + aa,2) + pow(aa - db,2) + add)) + 2.0e+00/sqrt(pow(r - aa,2) + pow(aa - db,2) + add) + 2.0e+00/sqrt(pow(r + aa,2) + pow(aa + db,2) + add) - 2.0e+00/sqrt(pow(r - aa,2) + pow(aa + db,2) + add) ;
                    charg = qxzdx/8.0e+00 ;
                }
                break ;
                case 2:
                if ( m == 0 )
                {
                    zzzz  = 1.0e+00/sqrt(pow(r - da - db,2) + add) + 1.0e+00/sqrt(pow(r + da + db,2) + add) + 1.0e+00/sqrt(pow(r - da + db,2) + add) + 1.0e+00/sqrt(pow(r + da - db,2) + add) -
                                                          2.0e+00/sqrt(pow(r - da,2) + db*db + add) - 2.0e+00/sqrt(pow(r - db,2) + da*da + add) - 2.0e+00/sqrt(pow(r + da,2) + db*db + add) -
                                                          2.0e+00/sqrt(pow(r + db,2) + da*da + add) + 2.0e+00/sqrt(r*r + pow(da - db,2) + add) + 2.0e+00/sqrt(r*r + pow(da + db,2) + add) ;
                    xyxy  = 4.0e+00/sqrt(r*r + pow(da - db,2) + add) + 4.0e+00/sqrt(r*r + pow(da + db,2) + add) - 8.0e+00/sqrt(r*r + da*da + db*db + add) ;
                    charg = zzzz/16.0e+00 - xyxy/64.0e+00 ;
                }
                else if ( m == 1 )
                {
                    aa = da/sqrt(2.0e+00) ;
                    ab = db/sqrt(2.0e+00) ;
                    qxzqxz = 2.0e+00/sqrt(pow(r + aa - ab,2) + pow(aa - ab,2) + add) - 2.0e+00/sqrt(pow(r + aa + ab,2) + pow(aa - ab,2) + add) - 2.0e+00/sqrt(pow(r - aa - ab,2) + pow(aa - ab,2) + add) +
                             2.0e+00/sqrt(pow(r - aa + ab,2) + pow(aa - ab,2) + add) - 2.0e+00/sqrt(pow(r + aa - ab,2) + pow(aa + ab,2) + add) + 2.0e+00/sqrt(pow(r + aa + ab,2) + pow(aa + ab,2) + add) +
                             2.0e+00/sqrt(pow(r - aa - ab,2) + pow(aa + ab,2) + add) - 2.0e+00/sqrt(pow(r - aa + ab,2) + pow(aa + ab,2) + add) ;
                    charg  = qxzqxz/16.0e+00 ;
                }
                else if ( m == 2 )
                {
                    xyxy  = 4.0e+00/sqrt(r*r + pow(da - db,2) + add) + 4.0e+00/sqrt(r*r + pow(da + db,2) + add) - 8.0e+00/sqrt(r*r + da*da + db*db + add) ;
                    charg = xyxy/16.0e+00 ;
                }
                break ;
            }
            break ;
    }
    return charg ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the derivative of the interaction of two point-charge configurations.
!---------------------------------------------------------------------------------------------------------------------------------*/
/*
! . r      distance.
! . l1,m   quantum numbers for multipole of configuration 1.
! . l2,m   quantum numbers for multipole of configuration 2.
! . da     charge separation of configuration 1.
! . db     charge separation of configuration 2.
! . add    additive term.
*/
static Real MNDOIntegralUtilities_2CChargeInteractionD ( const Real r, const Integer l1, const Integer l2, const Integer m, const Real da, const Real db, const Real add )
{
    Real aa, ab, dcharg, dxdx, dxqxz, dzdz, dzqzz, fac1, fac2, fac3, fac4, fac5, fac6, fac7, fac8, fac9, fac10, qqzz, qxzdx, qxzqxz, qzzdz, qzzq, xyxy, zzzz ;
    dcharg = 0.0e+00 ;
    switch ( l1 )
    {
        case 0:
            switch ( l2 )
            {
                case 0:
                    fac1   = r*r + add ;
                    dcharg = - r / ( fac1 * sqrt ( fac1 ) ) ;
                break ;
                case 1:
                    fac1    = pow(r + db,2) + add ;
                    fac2    = pow(r - db,2) + add ;
                    dcharg  = ( r + db ) / ( fac1 * sqrt ( fac1 ) ) - ( r - db ) / ( fac2 * sqrt ( fac2 ) ) ;
                    dcharg *= - 0.5e+00 ;
                break ;
                case 2:
                    fac1   = pow(r - db,2) + add ;
                    fac2   = r*r + db*db + add ;
                    fac3   = pow(r + db,2) + add ;
                    qqzz   = ( r - db ) / ( fac1 * sqrt ( fac1 ) ) - 2.0e+00 * r / ( fac2 * sqrt ( fac2 ) ) + ( r + db ) / ( fac3 * sqrt ( fac3 ) ) ;
                    dcharg = - qqzz / 4.0e+00 ;
                break ;
            }
            break ;
        case 1:
            switch ( l2 )
            {
                case 0:
                    fac1    = pow(r + da,2) + add ;
                    fac2    = pow(r - da,2) + add ;
                    dcharg  = - ( r + da ) / ( fac1 * sqrt ( fac1 ) ) + ( r - da ) / ( fac2 * sqrt ( fac2 ) ) ;
                    dcharg *= - 0.5e+00 ;
                break ;
                case 1:
                if ( m == 0 )
                {
                    fac1   = pow(r + da - db,2) + add ;
                    fac2   = pow(r - da + db,2) + add ;
                    fac3   = pow(r - da - db,2) + add ;
                    fac4   = pow(r + da + db,2) + add ;
                    dzdz   = ( r + da - db ) / ( fac1 * sqrt ( fac1 ) ) + ( r - da + db ) / ( fac2 * sqrt ( fac2 ) ) - ( r - da - db ) / ( fac3 * sqrt ( fac3 ) ) - ( r + da + db ) / ( fac4 * sqrt ( fac4 ) ) ;
                    dcharg = - dzdz / 4.0e+00 ;
                }
                else if ( m == 1 )
                {
                    fac1   = r*r + pow(da - db,2) + add ;
                    fac2   = r*r + pow(da + db,2) + add ;
                    dxdx   = 2.0e+00 * r / ( fac1 * sqrt ( fac1 ) ) - 2.0e+00 * r / ( fac2 * sqrt ( fac2 ) ) ;
                    dcharg = - dxdx / 4.0e+00 ;
                }
                break ;
                case 2:
                if ( m == 0 )
                {
                    fac1   = pow(r - da - db,2) + add ;
                    fac2   = pow(r - da,2) + db*db + add ;
                    fac3   = pow(r + db - da,2) + add ;
                    fac4   = pow(r - db + da,2) + add ;
                    fac5   = pow(r + da,2) + db*db + add ;
                    fac6   = pow(r + da + db,2) + add ;
                    dzqzz  = ( r - da - db ) / ( fac1 * sqrt ( fac1 ) ) - 2.0e+00 * ( r - da ) / ( fac2 * sqrt ( fac2 ) ) + ( r + db - da ) / ( fac3 * sqrt ( fac3 ) ) -
                             ( r - db + da ) / ( fac4 * sqrt ( fac4 ) ) + 2.0e+00 * ( r + da ) / ( fac5 * sqrt ( fac5 ) ) - ( r + da + db ) / ( fac6 * sqrt ( fac6 ) ) ;
                    dcharg = - dzqzz / 8.0e+00 ;
                }
                else if ( m == 1 )
                {
                    ab     = db / sqrt ( 2.0e+00 ) ;
                    fac1   = pow(r - ab,2) + pow(da - ab,2) + add ;
                    fac2   = pow(r + ab,2) + pow(da - ab,2) + add ;
                    fac3   = pow(r - ab,2) + pow(da + ab,2) + add ;
                    fac4   = pow(r + ab,2) + pow(da + ab,2) + add ;
                    dxqxz  = - 2.0e+00 * ( r - ab ) / ( fac1 * sqrt ( fac1 ) ) + 2.0e+00 * ( r + ab ) / ( fac2 * sqrt ( fac2 ) ) + 2.0e+00 * ( r - ab ) / ( fac3 * sqrt ( fac3 ) ) - 2.0e+00 * ( r + ab ) / ( fac4 * sqrt ( fac4 ) ) ;
                    dcharg = - dxqxz / 8.0e+00 ;
                }
                break ;
            }
            break ;
        case 2:
            switch ( l2 )
            {
                case 0:
                    fac1   = pow(r - da,2) + add ;
                    fac2   = r*r + da*da   + add ;
                    fac3   = pow(r + da,2) + add ;
                    qzzq   = ( r - da ) / ( fac1 * sqrt ( fac1 ) ) - 2.0e+00 * r / ( fac2 * sqrt ( fac2 ) ) + ( r + da ) / ( fac3 * sqrt ( fac3 ) ) ;
                    dcharg = - qzzq / 4.0e+00 ;
                break ;
                case 1:
                if ( m == 0 )
                {
                    fac1   = pow(r - da - db,2) + add ;
                    fac2   = pow(r - db,2) + da*da + add ;
                    fac3   = pow(r + da - db,2) + add ;
                    fac4   = pow(r - da + db,2) + add ;
                    fac5   = pow(r + db,2) + da*da + add ;
                    fac6   = pow(r + da + db,2) + add ;
                    qzzdz  = - ( r - da - db ) / ( fac1 * sqrt ( fac1 ) ) + 2.0e+00 * ( r - db ) / ( fac2 * sqrt ( fac2 ) ) - ( r + da - db ) / ( fac3 * sqrt ( fac3 ) ) +
                               ( r - da + db ) / ( fac4 * sqrt ( fac4 ) ) - 2.0e+00 * ( r + db ) / ( fac5 * sqrt ( fac5 ) ) + ( r + da + db ) / ( fac6 * sqrt ( fac6 ) ) ;
                    dcharg = - qzzdz / 8.0e+00 ;
                }
                else if ( m == 1 )
                {
                    aa     = da / sqrt ( 2.0e+00 ) ;
                    fac1   = pow(r + aa,2) + pow(aa - db,2) + add ;
                    fac2   = pow(r - aa,2) + pow(aa - db,2) + add ;
                    fac3   = pow(r + aa,2) + pow(aa + db,2) + add ;
                    fac4   = pow(r - aa,2) + pow(aa + db,2) + add ;
                    qxzdx  = -2.0e+00 * ( r + aa ) / ( fac1 * sqrt ( fac1 ) ) + 2.0e+00 * ( r - aa ) / ( fac2 * sqrt ( fac2 ) ) + 2.0e+00 * ( r + aa ) / ( fac3 * sqrt ( fac3 ) ) - 2.0e+00 * ( r - aa ) / ( fac4 * sqrt ( fac4 ) ) ;
                    dcharg = - qxzdx / 8.0e+00 ;
                }
                break ;
                case 2:
                if ( m == 0 )
                {
                    fac1   = pow(r - da - db,2) + add ;
                    fac2   = pow(r + da + db,2) + add ;
                    fac3   = pow(r - da + db,2) + add ;
                    fac4   = pow(r + da - db,2) + add ;
                    fac5   = pow(r - da,2) + db*db + add ;
                    fac6   = pow(r - db,2) + da*da + add ;
                    fac7   = pow(r + da,2) + db*db + add ;
                    fac8   = pow(r + db,2) + da*da + add ;
                    fac9   = r*r + pow(da - db,2) + add ;
                    fac10  = r*r + pow(da + db,2) + add ;
                    zzzz   = ( r - da - db ) / ( fac1 * sqrt ( fac1 ) ) + ( r + da + db ) / ( fac2 * sqrt ( fac2 ) ) +
                             ( r - da + db ) / ( fac3 * sqrt ( fac3 ) ) + ( r + da - db ) / ( fac4 * sqrt ( fac4 ) ) -
                             2.0e+00 * ( r - da ) / ( fac5 * sqrt ( fac5 ) ) - 2.0e+00 * ( r - db ) / ( fac6 * sqrt ( fac6 ) ) -
                             2.0e+00 * ( r + da ) / ( fac7 * sqrt ( fac7 ) ) - 2.0e+00 * ( r + db ) / ( fac8 * sqrt ( fac8 ) ) +
                             2.0e+00 * r / ( fac9 * sqrt ( fac9 ) )          + 2.0e+00 * r / ( fac10 * sqrt ( fac10 ) ) ;
                    fac1   = r*r + pow(da - db,2) + add ;
                    fac2   = r*r + pow(da + db,2) + add ;
                    fac3   = r*r + da*da + db*db  + add ;
                    xyxy   = 4.0e+00 * r / ( fac1 * sqrt ( fac1 ) ) + 4.0e+00 * r / ( fac2 * sqrt ( fac2 ) ) - 8.0e+00 * r / ( fac3 * sqrt ( fac3 ) ) ;
                    dcharg = - ( zzzz / 16.0e+00 - xyxy / 64.0e+00 ) ;
                }
                else if ( m == 1 )
                {
                    aa      = da / sqrt ( 2.0e+00 ) ;
                    ab      = db / sqrt ( 2.0e+00 ) ;
                    fac1    = pow(r + aa - ab,2) + pow(aa - ab,2) + add ;
                    fac2    = pow(r + aa + ab,2) + pow(aa - ab,2) + add ;
                    fac3    = pow(r - aa - ab,2) + pow(aa - ab,2) + add ;
                    fac4    = pow(r - aa + ab,2) + pow(aa - ab,2) + add ;
                    fac5    = pow(r + aa - ab,2) + pow(aa + ab,2) + add ;
                    fac6    = pow(r + aa + ab,2) + pow(aa + ab,2) + add ;
                    fac7    = pow(r - aa - ab,2) + pow(aa + ab,2) + add ;
                    fac8    = pow(r - aa + ab,2) + pow(aa + ab,2) + add ;
                    qxzqxz  = 2.0e+00 * ( r + aa - ab ) / ( fac1 * sqrt ( fac1 ) ) - 2.0e+00 * ( r + aa + ab ) / ( fac2 * sqrt ( fac2 ) ) - 2.0e+00 * ( r - aa - ab ) / ( fac3 * sqrt ( fac3 ) ) +
                              2.0e+00 * ( r - aa + ab ) / ( fac4 * sqrt ( fac4 ) ) - 2.0e+00 * ( r + aa - ab ) / ( fac5 * sqrt ( fac5 ) ) + 2.0e+00 * ( r + aa + ab ) / ( fac6 * sqrt ( fac6 ) ) +
                              2.0e+00 * ( r - aa - ab ) / ( fac7 * sqrt ( fac7 ) ) - 2.0e+00 * ( r - aa + ab ) / ( fac8 * sqrt ( fac8 ) ) ;
                    dcharg  = - qxzqxz / 16.0e+00 ;
                }
                else if ( m == 2 )
                {
                    fac1   = r*r + pow(da - db,2) + add ;
                    fac2   = r*r + pow(da + db,2) + add ;
                    fac3   = r*r + da*da + db*db  + add ;
                    xyxy   = 4.0e+00 * r / ( fac1 * sqrt ( fac1 ) ) + 4.0e+00 * r / ( fac2 * sqrt ( fac2 ) ) - 8.0e+00 * r / ( fac3 * sqrt ( fac3 ) ) ;
                    dcharg = - xyxy / 16.0e+00 ;
                }
                break ;
            }
            break ;
    }
    return dcharg ;
}
# endif
