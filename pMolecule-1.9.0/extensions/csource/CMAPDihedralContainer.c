/*------------------------------------------------------------------------------
! . File      : CMAPDihedralContainer.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
==================================================================================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>

# include "Boolean1DArray.h"
# include "CMAPDihedralContainer.h"
# include "Memory.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Integer CMAPDihedralTerm_Compare ( const void *vterm1, const void *vterm2 ) ;

/*==================================================================================================================================
! . Procedures.
==================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Activate terms.
!---------------------------------------------------------------------------------------------------------------------------------*/
void CMAPDihedralContainer_ActivateTerms ( CMAPDihedralContainer *self )
{
    if ( self != NULL )
    {
        auto Integer i ;
	for ( i = 0 ; i < self->nterms ; i++ ) self->terms[i].isActive = True ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
CMAPDihedralContainer *CMAPDihedralContainer_Allocate ( const Integer nterms, const Integer nparameters )
{
    CMAPDihedralContainer *self = NULL ;
    if ( ( nterms != 0 ) && ( nparameters != 0 ) )
    {
        auto Integer i ;
	self = ( CMAPDihedralContainer * ) Memory_Allocate ( sizeof ( CMAPDihedralContainer ) ) ;
        self->isSorted    = False       ;
	self->nterms      = nterms      ;
	self->nparameters = nparameters ;
	self->terms	  = ( CMAPDihedral   * ) Memory_Allocate_Array ( nterms,      sizeof ( CMAPDihedral    ) ) ;
	self->parameters  = ( BicubicSpline ** ) Memory_Allocate_Array ( nparameters, sizeof ( BicubicSpline * ) ) ;
	/* . Make all terms inactive. */
	for ( i = 0 ; i < nterms ; i++ ) self->terms[i].isActive = False ;
        /* . Nullify parameters. */
        for ( i = 0 ; i < nparameters ; i++ ) self->parameters[i] = NULL ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
CMAPDihedralContainer *CMAPDihedralContainer_Clone ( const CMAPDihedralContainer *self )
{
    CMAPDihedralContainer *new = NULL ;
    if ( self != NULL )
    {
        auto Integer i ;
        new = CMAPDihedralContainer_Allocate ( self->nterms, self->nparameters ) ;
        for ( i = 0 ; i < self->nterms ; i++ )
        {
            new->terms[i].isActive = self->terms[i].isActive ;
            new->terms[i].atom1    = self->terms[i].atom1    ;
            new->terms[i].atom2    = self->terms[i].atom2    ;
            new->terms[i].atom3    = self->terms[i].atom3    ;
            new->terms[i].atom4    = self->terms[i].atom4    ;
            new->terms[i].atom5    = self->terms[i].atom5    ;
            new->terms[i].atom6    = self->terms[i].atom6    ;
            new->terms[i].atom7    = self->terms[i].atom7    ;
            new->terms[i].atom8    = self->terms[i].atom8    ;
            new->terms[i].type     = self->terms[i].type     ;
        }
        for ( i = 0 ; i < self->nparameters ; i++ ) new->parameters[i] = BicubicSpline_Clone ( self->parameters[i], NULL ) ;
        new->isSorted = self->isSorted ;
    }
    return new ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deactivate terms between fixed atoms.
! . Already deactivated terms are not affected.
!---------------------------------------------------------------------------------------------------------------------------------*/
void CMAPDihedralContainer_DeactivateFixedAtomTerms ( CMAPDihedralContainer *self, Selection *fixedatoms )
{
    if ( ( self != NULL ) && ( fixedatoms != NULL ) )
    {
        auto Integer i, n ;
        n = CMAPDihedralContainer_UpperBound ( self ) ;
        Selection_MakeFlags ( fixedatoms, n ) ;
	for ( i = 0 ; i < self->nterms ; i++ )
	{
            if ( self->terms[i].isActive )
            {
                self->terms[i].isActive = ! ( fixedatoms->flags[self->terms[i].atom1] &&
                                              fixedatoms->flags[self->terms[i].atom2] &&
                                              fixedatoms->flags[self->terms[i].atom3] &&
                                              fixedatoms->flags[self->terms[i].atom4] &&
                                              fixedatoms->flags[self->terms[i].atom5] &&
                                              fixedatoms->flags[self->terms[i].atom6] &&
                                              fixedatoms->flags[self->terms[i].atom7] &&
                                              fixedatoms->flags[self->terms[i].atom8] ) ;
            }
	}
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deactivate terms involving QC atoms.
! . qcAtoms is the selection of both pure and boundary QC atoms.
! . Already deactivated terms are not affected.
!---------------------------------------------------------------------------------------------------------------------------------*/
void CMAPDihedralContainer_DeactivateQCAtomTerms ( CMAPDihedralContainer *self, Selection *qcAtoms, Selection *boundaryatoms )
{
    if ( ( self != NULL ) && ( qcAtoms != NULL ) )
    {
        auto Boolean  toExclude ;
        auto Integer  i, n ;
        n = CMAPDihedralContainer_UpperBound ( self ) ;
        Selection_MakeFlags ( qcAtoms, n ) ;
	for ( i = 0 ; i < self->nterms ; i++ )
	{
            if ( self->terms[i].isActive )
            {
                toExclude = ( qcAtoms->flags[self->terms[i].atom1] && qcAtoms->flags[self->terms[i].atom2] &&
                              qcAtoms->flags[self->terms[i].atom3] && qcAtoms->flags[self->terms[i].atom4] &&
                              qcAtoms->flags[self->terms[i].atom5] && qcAtoms->flags[self->terms[i].atom6] &&
                              qcAtoms->flags[self->terms[i].atom7] && qcAtoms->flags[self->terms[i].atom8] ) ;
                self->terms[i].isActive = ! toExclude ;
            }
	}
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void CMAPDihedralContainer_Deallocate ( CMAPDihedralContainer **self )
{
    if ( (*self) != NULL )
    {
        auto Integer i ;
        for ( i = 0 ; i < (*self)->nparameters ; i++ ) BicubicSpline_Deallocate ( &((*self)->parameters[i]) ) ;
        free ( (*self)->terms      ) ;
        free ( (*self)->parameters ) ;
        free ( (*self) ) ;
        (*self) = NULL   ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Energy and gradients.
! . Following Becker, Berendsen and van Gunsteren, JCC 16 p527 (1995).
!---------------------------------------------------------------------------------------------------------------------------------*/
# define LOWCOSPHI 0.5e+00
Real CMAPDihedralContainer_Energy ( const CMAPDihedralContainer *self, const Coordinates3 *coordinates3, Coordinates3 *gradients3 )
{
    Real energy = 0.0e+00 ;
    if ( ( self != NULL ) && ( coordinates3 != NULL ) )
    {
	auto Boolean doGradients ;
        auto Integer i1, i2, j1, j2, k1, k2, l1, l2, n, t ;
        auto Real    cosphi, dotij, dotlk, dtxi, dtyi, dtzi, dtxj, dtyj, dtzj, dtxk, dtyk, dtzk, dtxl, dtyl, dtzl, e, mn, sinphi, sx, sy, sz ;
        auto Real    df1, m21, mx1, my1, mz1, n21, nx1, ny1, nz1, phi1, rkj1, rkj21, xij1, yij1, zij1, xkj1, ykj1, zkj1, xlk1, ylk1, zlk1 ;
        auto Real    df2, m22, mx2, my2, mz2, n22, nx2, ny2, nz2, phi2, rkj2, rkj22, xij2, yij2, zij2, xkj2, ykj2, zkj2, xlk2, ylk2, zlk2 ;
	doGradients = ( gradients3 != NULL ) ;
	for ( n = 0 ; n < self->nterms ; n++ )
	{
	    if ( self->terms[n].isActive )
	    {
                /* . Local data. */
	        i1 = self->terms[n].atom1 ;
	        j1 = self->terms[n].atom2 ;
	        k1 = self->terms[n].atom3 ;
	        l1 = self->terms[n].atom4 ;
	        i2 = self->terms[n].atom5 ;
	        j2 = self->terms[n].atom6 ;
	        k2 = self->terms[n].atom7 ;
	        l2 = self->terms[n].atom8 ;
	        t  = self->terms[n].type  ;
                /* . First dihedral. */
	        /* . Coordinate displacements. */
	        Coordinates3_DifferenceRow ( coordinates3, i1, j1, xij1, yij1, zij1 ) ;
	        Coordinates3_DifferenceRow ( coordinates3, k1, j1, xkj1, ykj1, zkj1 ) ;
	        Coordinates3_DifferenceRow ( coordinates3, l1, k1, xlk1, ylk1, zlk1 ) ;
	        rkj21 = xkj1 * xkj1 + ykj1 * ykj1 + zkj1 * zkj1 ;
	        rkj1  = sqrt ( rkj21 ) ;
                /* . m and n. */
                mx1 = yij1 * zkj1 - zij1 * ykj1 ;
	        my1 = zij1 * xkj1 - xij1 * zkj1 ;
	        mz1 = xij1 * ykj1 - yij1 * xkj1 ;
	        nx1 = ylk1 * zkj1 - zlk1 * ykj1 ;
	        ny1 = zlk1 * xkj1 - xlk1 * zkj1 ;
	        nz1 = xlk1 * ykj1 - ylk1 * xkj1 ;
                m21 = mx1 * mx1 + my1 * my1 + mz1 * mz1 ;
	        n21 = nx1 * nx1 + ny1 * ny1 + nz1 * nz1 ;
	        mn  = sqrt ( m21 * n21 ) ;
                /* . Cosine and sine of the dihedral. */
                cosphi =        (  mx1 * nx1 +  my1 * ny1 +  mz1 * nz1 ) / mn ;
                sinphi = rkj1 * ( xij1 * nx1 + yij1 * ny1 + zij1 * nz1 ) / mn ;
                /* . Dihedral - follow CHARMM. */
                if ( ( cosphi < -LOWCOSPHI ) || ( cosphi > LOWCOSPHI ) )
                {
                    phi1 = asin ( sinphi ) ;
                    if ( cosphi < 0.0e+00 )
                    {
                        if ( phi1 > 0.0e+00 ) phi1 =     M_PI - phi1 ;
                        else                  phi1 = - ( M_PI + phi1 ) ;
                    }
                }
                else
                {
                    phi1 = acos ( cosphi ) ;
                    if ( sinphi < 0.0e+00 ) phi1 *= -1.e+00 ;
                }
                /* . Second dihedral. */
	        /* . Coordinate displacements. */
	        Coordinates3_DifferenceRow ( coordinates3, i2, j2, xij2, yij2, zij2 ) ;
	        Coordinates3_DifferenceRow ( coordinates3, k2, j2, xkj2, ykj2, zkj2 ) ;
	        Coordinates3_DifferenceRow ( coordinates3, l2, k2, xlk2, ylk2, zlk2 ) ;
	        rkj22 = xkj2 * xkj2 + ykj2 * ykj2 + zkj2 * zkj2 ;
	        rkj2  = sqrt ( rkj22 ) ;
                /* . m and n. */
                mx2 = yij2 * zkj2 - zij2 * ykj2 ;
	        my2 = zij2 * xkj2 - xij2 * zkj2 ;
	        mz2 = xij2 * ykj2 - yij2 * xkj2 ;
	        nx2 = ylk2 * zkj2 - zlk2 * ykj2 ;
	        ny2 = zlk2 * xkj2 - xlk2 * zkj2 ;
	        nz2 = xlk2 * ykj2 - ylk2 * xkj2 ;
                m22 = mx2 * mx2 + my2 * my2 + mz2 * mz2 ;
	        n22 = nx2 * nx2 + ny2 * ny2 + nz2 * nz2 ;
	        mn  = sqrt ( m22 * n22 ) ;
                /* . Cosine and sine of the dihedral. */
                cosphi =        (  mx2 * nx2 +  my2 * ny2 +  mz2 * nz2 ) / mn ;
                sinphi = rkj2 * ( xij2 * nx2 + yij2 * ny2 + zij2 * nz2 ) / mn ;
                /* . Dihedral - follow CHARMM. */
                if ( ( cosphi < -LOWCOSPHI ) || ( cosphi > LOWCOSPHI ) )
                {
                    phi2 = asin ( sinphi ) ;
                    if ( cosphi < 0.0e+00 )
                    {
                        if ( phi2 > 0.0e+00 ) phi2 =     M_PI - phi2 ;
                        else                  phi2 = - ( M_PI + phi2 ) ;
                    }
                }
                else
                {
                    phi2 = acos ( cosphi ) ;
                    if ( sinphi < 0.0e+00 ) phi2 *= -1.e+00 ;
                }
                /* . The energy term. */
                BicubicSpline_Evaluate ( self->parameters[t], phi1, phi2, &e, &df1, &df2, NULL ) ;
                energy += e ;
/*
{
auto Real em1, em2, ep1, ep2 ;
# define _STEP 1.0e-4
printf ( "EC> i1, j1, k1, l1, i2, j2, k2, l2, t, phi1, phi2, e, df1, df2 = %d %d %d %d %d %d %d %d %d %20.5f %20.5f  %20.5f %20.5f %20.5f\n", i1, j1, k1, l1, i2, j2, k2, l2, t, phi1, phi2, e, df1, df2 ) ;
BicubicSpline_Evaluate ( self->parameters[t], phi1-_STEP, phi2, &em1, NULL, NULL, NULL ) ;
BicubicSpline_Evaluate ( self->parameters[t], phi1+_STEP, phi2, &ep1, NULL, NULL, NULL ) ;
BicubicSpline_Evaluate ( self->parameters[t], phi1, phi2-_STEP, &em2, NULL, NULL, NULL ) ;
BicubicSpline_Evaluate ( self->parameters[t], phi1, phi2+_STEP, &ep2, NULL, NULL, NULL ) ;
printf ( "\nDerivative 1 = %20.5f %20.5f %20.5f\n", df1, (ep1-em1)/(2.0e+00 * _STEP),fabs ( df1 - (ep1-em1)/(2.0e+00 * _STEP) ) ) ;
printf ( "\nDerivative 2 = %20.5f %20.5f %20.5f\n", df2, (ep2-em2)/(2.0e+00 * _STEP),fabs ( df2 - (ep2-em2)/(2.0e+00 * _STEP) ) ) ;
}
*/
                if ( doGradients )
                {
	            /* . The derivatives - first dihedral. */
                    /* . i and l. */
                    dtxi =   df1 * rkj1 * mx1 / m21 ;
                    dtyi =   df1 * rkj1 * my1 / m21 ;
                    dtzi =   df1 * rkj1 * mz1 / m21 ;
                    dtxl = - df1 * rkj1 * nx1 / n21 ;
                    dtyl = - df1 * rkj1 * ny1 / n21 ;
                    dtzl = - df1 * rkj1 * nz1 / n21 ;
                    /* . j and k. */
                    dotij = xij1 * xkj1 + yij1 * ykj1 + zij1 * zkj1 ;
                    dotlk = xlk1 * xkj1 + ylk1 * ykj1 + zlk1 * zkj1 ;
	            sx    = ( dotij * dtxi + dotlk * dtxl ) / rkj21 ;
	            sy    = ( dotij * dtyi + dotlk * dtyl ) / rkj21 ;
	            sz    = ( dotij * dtzi + dotlk * dtzl ) / rkj21 ;
	            dtxj  =   sx - dtxi ;
	            dtyj  =   sy - dtyi ;
	            dtzj  =   sz - dtzi ;
	            dtxk  = - sx - dtxl ;
	            dtyk  = - sy - dtyl ;
	            dtzk  = - sz - dtzl ;
                    /* . Add in the contributions. */
	            Coordinates3_IncrementRow ( gradients3, i1, dtxi, dtyi, dtzi ) ;
	            Coordinates3_IncrementRow ( gradients3, j1, dtxj, dtyj, dtzj ) ;
	            Coordinates3_IncrementRow ( gradients3, k1, dtxk, dtyk, dtzk ) ;
	            Coordinates3_IncrementRow ( gradients3, l1, dtxl, dtyl, dtzl ) ;
	            /* . The derivatives - second dihedral. */
                    /* . i and l. */
                    dtxi =   df2 * rkj2 * mx2 / m22 ;
                    dtyi =   df2 * rkj2 * my2 / m22 ;
                    dtzi =   df2 * rkj2 * mz2 / m22 ;
                    dtxl = - df2 * rkj2 * nx2 / n22 ;
                    dtyl = - df2 * rkj2 * ny2 / n22 ;
                    dtzl = - df2 * rkj2 * nz2 / n22 ;
                    /* . j and k. */
                    dotij = xij2 * xkj2 + yij2 * ykj2 + zij2 * zkj2 ;
                    dotlk = xlk2 * xkj2 + ylk2 * ykj2 + zlk2 * zkj2 ;
	            sx    = ( dotij * dtxi + dotlk * dtxl ) / rkj22 ;
	            sy    = ( dotij * dtyi + dotlk * dtyl ) / rkj22 ;
	            sz    = ( dotij * dtzi + dotlk * dtzl ) / rkj22 ;
	            dtxj  =   sx - dtxi ;
	            dtyj  =   sy - dtyi ;
	            dtzj  =   sz - dtzi ;
	            dtxk  = - sx - dtxl ;
	            dtyk  = - sy - dtyl ;
	            dtzk  = - sz - dtzl ;
                    /* . Add in the contributions. */
	            Coordinates3_IncrementRow ( gradients3, i2, dtxi, dtyi, dtzi ) ;
	            Coordinates3_IncrementRow ( gradients3, j2, dtxj, dtyj, dtzj ) ;
	            Coordinates3_IncrementRow ( gradients3, k2, dtxk, dtyk, dtzk ) ;
	            Coordinates3_IncrementRow ( gradients3, l2, dtxl, dtyl, dtzl ) ;
        	}
	    }
	}
    }
    return energy ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Merging.
!---------------------------------------------------------------------------------------------------------------------------------*/
CMAPDihedralContainer *CMAPDihedralContainer_Merge ( const CMAPDihedralContainer *self, const CMAPDihedralContainer *other, const Integer atomincrement )
{
    CMAPDihedralContainer *new = NULL ;
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        auto Integer i ;
        new = CMAPDihedralContainer_Allocate ( self->nterms + other->nterms, self->nparameters + other->nparameters ) ;
        for ( i = 0 ; i < self->nterms ; i++ )
        {
            new->terms[i].isActive = self->terms[i].isActive ;
            new->terms[i].atom1    = self->terms[i].atom1    ;
            new->terms[i].atom2    = self->terms[i].atom2    ;
            new->terms[i].atom3    = self->terms[i].atom3    ;
            new->terms[i].atom4    = self->terms[i].atom4    ;
            new->terms[i].atom5    = self->terms[i].atom5    ;
            new->terms[i].atom6    = self->terms[i].atom6    ;
            new->terms[i].atom7    = self->terms[i].atom7    ;
            new->terms[i].atom8    = self->terms[i].atom8    ;
            new->terms[i].type     = self->terms[i].type     ;
        }
        for ( i = 0 ; i < other->nterms ; i++ )
        {
            new->terms[i+self->nterms].isActive = other->terms[i].isActive ;
            new->terms[i+self->nterms].atom1    = other->terms[i].atom1 + atomincrement     ;
            new->terms[i+self->nterms].atom2    = other->terms[i].atom2 + atomincrement     ;
            new->terms[i+self->nterms].atom3    = other->terms[i].atom3 + atomincrement     ;
            new->terms[i+self->nterms].atom4    = other->terms[i].atom4 + atomincrement     ;
            new->terms[i+self->nterms].atom5    = other->terms[i].atom5 + atomincrement     ;
            new->terms[i+self->nterms].atom6    = other->terms[i].atom6 + atomincrement     ;
            new->terms[i+self->nterms].atom7    = other->terms[i].atom7 + atomincrement     ;
            new->terms[i+self->nterms].atom8    = other->terms[i].atom8 + atomincrement     ;
            new->terms[i+self->nterms].type     = other->terms[i].type  + self->nparameters ;
        }
        for ( i = 0 ; i < self->nparameters  ; i++ ) new->parameters[i]                   = BicubicSpline_Clone ( self->parameters [i], NULL ) ;
        for ( i = 0 ; i < other->nparameters ; i++ ) new->parameters[i+self->nparameters] = BicubicSpline_Clone ( other->parameters[i], NULL ) ;
        new->isSorted = ( self->isSorted && other->isSorted ) ;
    }
    return new ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Return the number of inactive terms.
!---------------------------------------------------------------------------------------------------------------------------------*/
int CMAPDihedralContainer_NumberOfInactiveTerms ( const CMAPDihedralContainer *self )
{
    Integer n = 0 ;
    if ( self != NULL )
    {
        auto Integer i ;
	for ( i = 0 ; i < self->nterms ; i++ ) if ( ! self->terms[i].isActive ) n++ ;
    }
    return n ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Pruning.
! . Only terms are pruned.
!---------------------------------------------------------------------------------------------------------------------------------*/
CMAPDihedralContainer *CMAPDihedralContainer_Prune ( CMAPDihedralContainer *self, Selection *selection )
{
    CMAPDihedralContainer *new = NULL ;
    if ( ( self != NULL ) && ( selection != NULL ) )
    {
        auto Boolean1DArray *flags ;
        auto Integer         i, n  ;
        n = CMAPDihedralContainer_UpperBound ( self ) ;
        Selection_MakeFlags     ( selection, n ) ;
        Selection_MakePositions ( selection, n ) ;
        flags = Boolean1DArray_Allocate ( self->nterms, NULL ) ;
	for ( i = 0, n = 0 ; i < self->nterms ; i++ )
	{
            Boolean1DArray_Item ( flags, i ) = ( selection->flags[self->terms[i].atom1] && selection->flags[self->terms[i].atom2] &&
                                                 selection->flags[self->terms[i].atom3] && selection->flags[self->terms[i].atom4] &&
                                                 selection->flags[self->terms[i].atom5] && selection->flags[self->terms[i].atom6] &&
                                                 selection->flags[self->terms[i].atom7] && selection->flags[self->terms[i].atom8] ) ;
            if ( Boolean1DArray_Item ( flags, i ) ) n++ ;
	}
	if ( n > 0 )
	{
            new = CMAPDihedralContainer_Allocate ( n, self->nparameters ) ;
            for ( i = 0 ; i < self->nparameters ; i++ ) new->parameters[i] = BicubicSpline_Clone ( self->parameters[i], NULL ) ;
            for ( i = 0, n = 0 ; i < self->nterms ; i++ )
            {
                if ( Boolean1DArray_Item ( flags, i ) )
                {
        	    new->terms[n].isActive =                      self->terms[i].isActive ;
        	    new->terms[n].atom1    = selection->positions[self->terms[i].atom1]   ;
        	    new->terms[n].atom2    = selection->positions[self->terms[i].atom2]   ;
        	    new->terms[n].atom3    = selection->positions[self->terms[i].atom3]   ;
        	    new->terms[n].atom4    = selection->positions[self->terms[i].atom4]   ;
        	    new->terms[n].atom5    = selection->positions[self->terms[i].atom5]   ;
        	    new->terms[n].atom6    = selection->positions[self->terms[i].atom6]   ;
        	    new->terms[n].atom7    = selection->positions[self->terms[i].atom7]   ;
        	    new->terms[n].atom8    = selection->positions[self->terms[i].atom8]   ;
        	    new->terms[n].type     =                      self->terms[i].type     ;
        	    n++ ;
                }
            }
            new->isSorted = self->isSorted ;
	}
	Boolean1DArray_Deallocate ( &flags ) ;
    }
    return new ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Sorting.
! . Within a harmonic improper, atom2 > atom3.
! . Within the array, ordering is done with increased values of atom2 and then
! . atom3 and then atom1 and then atom4.
! . Duplicates are not removed.
!---------------------------------------------------------------------------------------------------------------------------------*/
void CMAPDihedralContainer_Sort ( CMAPDihedralContainer *self )
{
    if ( ( self != NULL ) && ( ! self->isSorted ) )
    {
        auto Integer atom1, atom2, atom3, atom4, i ;
        /* . Order atom2 and atom3 within each term. */
        for ( i = 0 ; i < self->nterms ; i++ )
        {
            atom1 = self->terms[i].atom1 ;
            atom2 = self->terms[i].atom2 ;
            atom3 = self->terms[i].atom3 ;
            atom4 = self->terms[i].atom4 ;
            if ( atom3 > atom2 )
            {
                self->terms[i].atom1 = atom4 ;
                self->terms[i].atom2 = atom3 ;
                self->terms[i].atom3 = atom2 ;
                self->terms[i].atom4 = atom1 ;
            }
            atom1 = self->terms[i].atom5 ;
            atom2 = self->terms[i].atom6 ;
            atom3 = self->terms[i].atom7 ;
            atom4 = self->terms[i].atom8 ;
            if ( atom3 > atom2 )
            {
                self->terms[i].atom5 = atom4 ;
                self->terms[i].atom6 = atom3 ;
                self->terms[i].atom7 = atom2 ;
                self->terms[i].atom8 = atom1 ;
            }
        }
        /* . Order the terms within the container. */
        qsort ( ( void * ) self->terms, ( size_t ) self->nterms, sizeof ( CMAPDihedral ), ( void * ) CMAPDihedralTerm_Compare ) ;
        self->isSorted = True ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Return the upper bound for the container.
! . This is the value of the largest index plus one.
!---------------------------------------------------------------------------------------------------------------------------------*/
int CMAPDihedralContainer_UpperBound ( CMAPDihedralContainer *self )
{
    Integer upperBound = 0 ;
    if ( ( self != NULL ) && ( self->nterms > 0 ) )
    {
        CMAPDihedralContainer_Sort ( self ) ;
        upperBound  = Maximum ( self->terms[self->nterms-1].atom1, self->terms[self->nterms-1].atom2 ) ;
        upperBound  = Maximum ( upperBound, self->terms[self->nterms-1].atom4 ) ;
        upperBound  = Maximum ( upperBound, self->terms[self->nterms-1].atom5 ) ;
        upperBound  = Maximum ( upperBound, self->terms[self->nterms-1].atom6 ) ;
        upperBound  = Maximum ( upperBound, self->terms[self->nterms-1].atom8 ) ;
        upperBound += 1 ;
    }
    return upperBound ;
}

/*==================================================================================================================================
! . Private procedures.
!=================================================================================================================================*/
static Integer CMAPDihedralTerm_Compare ( const void *vterm1, const void *vterm2 )
{
    CMAPDihedral *term1, *term2 ;
    Integer       i ;
    term1 = ( CMAPDihedral * ) vterm1 ;
    term2 = ( CMAPDihedral * ) vterm2 ;
         if ( term1->atom2 < term2->atom2 ) i = -1 ;
    else if ( term1->atom2 > term2->atom2 ) i =  1 ;
    else if ( term1->atom3 < term2->atom3 ) i = -1 ;
    else if ( term1->atom3 > term2->atom3 ) i =  1 ;
    else if ( term1->atom1 < term2->atom1 ) i = -1 ;
    else if ( term1->atom1 > term2->atom1 ) i =  1 ;
    else if ( term1->atom4 < term2->atom4 ) i = -1 ;
    else if ( term1->atom4 > term2->atom4 ) i =  1 ;
    else if ( term1->atom6 < term2->atom6 ) i = -1 ;
    else if ( term1->atom6 > term2->atom6 ) i =  1 ;
    else if ( term1->atom7 < term2->atom7 ) i = -1 ;
    else if ( term1->atom7 > term2->atom7 ) i =  1 ;
    else if ( term1->atom5 < term2->atom5 ) i = -1 ;
    else if ( term1->atom5 > term2->atom5 ) i =  1 ;
    else if ( term1->atom8 < term2->atom8 ) i = -1 ;
    else if ( term1->atom8 > term2->atom8 ) i =  1 ;
    else if ( term1->type  < term2->type  ) i = -1 ;
    else if ( term1->type  > term2->type  ) i =  1 ;
    else if ( ! term1->isActive && term2->isActive ) i = -1 ;
    else if ( term1->isActive && ! term2->isActive ) i =  1 ;
    else i = 0 ;
    return i ;
}
