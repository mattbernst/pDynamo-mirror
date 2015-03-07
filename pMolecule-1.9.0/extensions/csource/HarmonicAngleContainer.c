/*------------------------------------------------------------------------------
! . File      : HarmonicAngleContainer.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==============================================================================
!=============================================================================*/

# include <math.h>
# include <stdlib.h>

# include "HarmonicAngleContainer.h"
# include "Memory.h"

/* . Use cos(t-t0)/sin(t-t0) in calculation and normalize from c^2+s^2? */

/*------------------------------------------------------------------------------
! . Local procedures.
!-----------------------------------------------------------------------------*/
static Integer HarmonicAngleTerm_Compare ( const void *vterm1, const void *vterm2 ) ;

/*------------------------------------------------------------------------------
! . Definitions.
!-----------------------------------------------------------------------------*/
/* . The tolerance for calculating angle linearity. */
# define DOT_LIMIT 0.999999

/*==============================================================================
! . Procedures.
!=============================================================================*/
/*------------------------------------------------------------------------------
! . Activate terms.
!-----------------------------------------------------------------------------*/
void HarmonicAngleContainer_ActivateTerms ( HarmonicAngleContainer *self )
{
    if ( self != NULL )
    {
        auto Integer i ;
	for ( i = 0 ; i < self->nterms ; i++ ) self->terms[i].QACTIVE = True ;
    }
}

/*------------------------------------------------------------------------------
! . Allocation.
!-----------------------------------------------------------------------------*/
HarmonicAngleContainer *HarmonicAngleContainer_Allocate ( const Integer nterms, const Integer nparameters )
{
    HarmonicAngleContainer *self = NULL ;
    if ( ( nterms != 0 ) && ( nparameters != 0 ) )
    {
        Integer i ;
	self = ( HarmonicAngleContainer * ) Memory_Allocate ( sizeof ( HarmonicAngleContainer ) ) ;
        self->QSORTED     = False       ;
	self->nterms      = nterms      ;
	self->nparameters = nparameters ;
	self->terms	  = ( HarmonicAngle * )          Memory_Allocate_Array ( nterms,      sizeof ( HarmonicAngle )          ) ;
	self->parameters  = ( HarmonicAngleParameter * ) Memory_Allocate_Array ( nparameters, sizeof ( HarmonicAngleParameter ) ) ;
	/* . Make all terms inactive. */
	for ( i = 0 ; i < nterms ; i++ ) self->terms[i].QACTIVE = False ;
    }
    return self ;
}

/*------------------------------------------------------------------------------
! . Cloning.
!-----------------------------------------------------------------------------*/
HarmonicAngleContainer *HarmonicAngleContainer_Clone ( const HarmonicAngleContainer *self )
{
    HarmonicAngleContainer *new = NULL ;
    if ( self != NULL )
    {
        auto Integer i ;
        new = HarmonicAngleContainer_Allocate ( self->nterms, self->nparameters ) ;
        for ( i = 0 ; i < self->nterms ; i++ )
        {
            new->terms[i].QACTIVE = self->terms[i].QACTIVE ;
            new->terms[i].atom1   = self->terms[i].atom1   ;
            new->terms[i].atom2   = self->terms[i].atom2   ;
            new->terms[i].atom3   = self->terms[i].atom3   ;
            new->terms[i].type    = self->terms[i].type    ;
        }
        for ( i = 0 ; i < self->nparameters ; i++ )
        {
            new->parameters[i].eq = self->parameters[i].eq ;
            new->parameters[i].fc = self->parameters[i].fc ;
        }
        new->QSORTED = self->QSORTED ;
    }
    return new ;
}

/*------------------------------------------------------------------------------
! . Deactivate terms between fixed atoms.
! . Already deactivated terms are not affected.
!-----------------------------------------------------------------------------*/
void HarmonicAngleContainer_DeactivateFixedAtomTerms ( HarmonicAngleContainer *self, Selection *fixedatoms )
{
    if ( ( self != NULL ) && ( fixedatoms != NULL ) )
    {
        auto Integer i, n ;
        n = HarmonicAngleContainer_UpperBound ( self ) ;
        Selection_MakeFlags ( fixedatoms, n ) ;
	for ( i = 0 ; i < self->nterms ; i++ )
	{
            if ( self->terms[i].QACTIVE )
            {
                self->terms[i].QACTIVE = ! ( fixedatoms->flags[self->terms[i].atom1] &&
                                             fixedatoms->flags[self->terms[i].atom2] &&
                                             fixedatoms->flags[self->terms[i].atom3] ) ;
            }
	}
    }
}

/*------------------------------------------------------------------------------
! . Deactivate terms involving QC atoms.
! . qcAtoms is the selection of both pure and boundary QC atoms.
! . Already deactivated terms are not affected.
!-----------------------------------------------------------------------------*/
void HarmonicAngleContainer_DeactivateQCAtomTerms ( HarmonicAngleContainer *self, Selection *qcAtoms, Selection *boundaryatoms )
{
    if ( ( self != NULL ) && ( qcAtoms != NULL ) )
    {
        auto Boolean QEXCLUDE ;
        auto Integer  i, n ;
        n = HarmonicAngleContainer_UpperBound ( self ) ;
        Selection_MakeFlags ( qcAtoms, n ) ;
	for ( i = 0 ; i < self->nterms ; i++ )
	{
            if ( self->terms[i].QACTIVE )
            {
                QEXCLUDE = ( qcAtoms->flags[self->terms[i].atom1] && qcAtoms->flags[self->terms[i].atom2] && qcAtoms->flags[self->terms[i].atom3] ) ;
                self->terms[i].QACTIVE = ! QEXCLUDE ;
            }
	}
    }
}

/*------------------------------------------------------------------------------
! . Deallocation.
!-----------------------------------------------------------------------------*/
void HarmonicAngleContainer_Deallocate ( HarmonicAngleContainer **self )
{
    if ( (*self) != NULL )\
    {
        free ( (*self)->terms      ) ;
        free ( (*self)->parameters ) ;
        free ( (*self) ) ;
        (*self) = NULL   ;
    }
}

/*------------------------------------------------------------------------------
! . Energy and gradients.
!-----------------------------------------------------------------------------*/
double HarmonicAngleContainer_Energy ( const HarmonicAngleContainer *self, const Coordinates3 *coordinates3, Coordinates3 *gradients3 )
{
    Real energy = 0.0e+00 ;
    if ( ( self != NULL ) && ( coordinates3 != NULL ) )
    {
	auto Boolean   QGRADIENTS ;
        auto Real df, disp, dot, dtdx, dtxi, dtxk, dtyi, dtyk, dtzi, dtzk, rij, rkj, theta, xij, yij, zij, xkj, ykj, zkj ;
        auto Integer    i, j, k, n, t ;
	QGRADIENTS = ( gradients3 != NULL ) ;
	for ( n = 0 ; n < self->nterms ; n++ )
	{
	    if ( self->terms[n].QACTIVE )
	    {
		i = self->terms[n].atom1 ;
		j = self->terms[n].atom2 ;
	        k = self->terms[n].atom3 ;
		t = self->terms[n].type  ;
	        Coordinates3_DifferenceRow ( coordinates3, i, j, xij, yij, zij ) ;
	        Coordinates3_DifferenceRow ( coordinates3, k, j, xkj, ykj, zkj ) ;
	        rij = sqrt ( xij * xij + yij * yij + zij * zij ) ;
	        rkj = sqrt ( xkj * xkj + ykj * ykj + zkj * zkj ) ;
	        xij /= rij ; yij /= rij ; zij /= rij ;
	        xkj /= rkj ; ykj /= rkj ; zkj /= rkj ;
	        dot   = xij * xkj + yij * ykj + zij * zkj ;
                dot   = Maximum ( - DOT_LIMIT, dot ) ;
	        dot   = Minimum (   DOT_LIMIT, dot ) ;
                theta = acos ( dot ) ;
	        disp  = theta - self->parameters[t].eq ;
	        df    = self->parameters[t].fc * disp  ;
	        energy += ( df * disp ) ;
                if ( QGRADIENTS )
                {
	            dtdx = - 1.0 / sqrt ( 1.0 - dot * dot ) ;
                    df  *= ( 2.0e+00 * dtdx ) ;
                    dtxi = df * ( xkj - dot * xij ) / rij ;
                    dtyi = df * ( ykj - dot * yij ) / rij ;
                    dtzi = df * ( zkj - dot * zij ) / rij ;
                    dtxk = df * ( xij - dot * xkj ) / rkj ;
                    dtyk = df * ( yij - dot * ykj ) / rkj ;
                    dtzk = df * ( zij - dot * zkj ) / rkj ;
	            Coordinates3_IncrementRow ( gradients3, i,   dtxi,            dtyi,            dtzi          ) ;
	            Coordinates3_IncrementRow ( gradients3, k,          dtxk,            dtyk,            dtzk   ) ;
	            Coordinates3_DecrementRow ( gradients3, j, ( dtxi + dtxk ), ( dtyi + dtyk ), ( dtzi + dtzk ) ) ;
        	}
	    }
	}
    }
    return energy ;
}

/*------------------------------------------------------------------------------
! . Merging.
!-----------------------------------------------------------------------------*/
HarmonicAngleContainer *HarmonicAngleContainer_Merge ( const HarmonicAngleContainer *self, const HarmonicAngleContainer *other, const Integer atomincrement )
{
    HarmonicAngleContainer *new = NULL ;
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        auto Integer i ;
        new = HarmonicAngleContainer_Allocate ( self->nterms + other->nterms, self->nparameters + other->nparameters ) ;
        for ( i = 0 ; i < self->nterms ; i++ )
        {
            new->terms[i].QACTIVE = self->terms[i].QACTIVE ;
            new->terms[i].atom1   = self->terms[i].atom1   ;
            new->terms[i].atom2   = self->terms[i].atom2   ;
            new->terms[i].atom3   = self->terms[i].atom3   ;
            new->terms[i].type    = self->terms[i].type    ;
        }
        for ( i = 0 ; i < other->nterms ; i++ )
        {
            new->terms[i+self->nterms].QACTIVE = other->terms[i].QACTIVE ;
            new->terms[i+self->nterms].atom1   = other->terms[i].atom1 + atomincrement     ;
            new->terms[i+self->nterms].atom2   = other->terms[i].atom2 + atomincrement     ;
            new->terms[i+self->nterms].atom3   = other->terms[i].atom3 + atomincrement     ;
            new->terms[i+self->nterms].type    = other->terms[i].type  + self->nparameters ;
        }
        for ( i = 0 ; i < self->nparameters ; i++ )
        {
            new->parameters[i].eq = self->parameters[i].eq ;
            new->parameters[i].fc = self->parameters[i].fc ;
        }
        for ( i = 0 ; i < other->nparameters ; i++ )
        {
            new->parameters[i+self->nparameters].eq = other->parameters[i].eq ;
            new->parameters[i+self->nparameters].fc = other->parameters[i].fc ;
        }
        new->QSORTED = ( self->QSORTED && other->QSORTED ) ;
    }
    return new ;
}

/*------------------------------------------------------------------------------
! . Return the number of inactive terms.
!-----------------------------------------------------------------------------*/
int HarmonicAngleContainer_NumberOfInactiveTerms ( const HarmonicAngleContainer *self )
{
    Integer n = 0 ;
    if ( self != NULL )
    {
        auto Integer i ;
	for ( i = 0 ; i < self->nterms ; i++ ) if ( ! self->terms[i].QACTIVE ) n++ ;
    }
    return n ;
}

/*------------------------------------------------------------------------------
! . Pruning.
! . Only terms are pruned.
!-----------------------------------------------------------------------------*/
HarmonicAngleContainer *HarmonicAngleContainer_Prune ( HarmonicAngleContainer *self, Selection *selection )
{
    HarmonicAngleContainer *new = NULL ;
    if ( ( self != NULL ) && ( selection != NULL ) )
    {
        auto Boolean *flags ;
        auto Integer   i, n  ;
        n = HarmonicAngleContainer_UpperBound ( self ) ;
        Selection_MakeFlags     ( selection, n ) ;
        Selection_MakePositions ( selection, n ) ;
        flags = Memory_Allocate_Array_Boolean ( self->nterms ) ;
	for ( i = 0, n = 0 ; i < self->nterms ; i++ )
	{
            flags[i] = ( selection->flags[self->terms[i].atom1] && selection->flags[self->terms[i].atom2]  && selection->flags[self->terms[i].atom3] ) ;
            if ( flags[i] ) n++ ;
	}
	if ( n > 0 )
	{
            new = HarmonicAngleContainer_Allocate ( n, self->nparameters ) ;
            for ( i = 0 ; i < self->nparameters ; i++ )
            {
                new->parameters[i].eq = self->parameters[i].eq ;
                new->parameters[i].fc = self->parameters[i].fc ;
            }
            for ( i = 0, n = 0 ; i < self->nterms ; i++ )
            {
                if ( flags[i] )
                {
        	    new->terms[n].QACTIVE =                      self->terms[i].QACTIVE ;
        	    new->terms[n].atom1   = selection->positions[self->terms[i].atom1]  ;
        	    new->terms[n].atom2   = selection->positions[self->terms[i].atom2]  ;
        	    new->terms[n].atom3   = selection->positions[self->terms[i].atom3]  ;
        	    new->terms[n].type    =                      self->terms[i].type    ;
        	    n++ ;
                }
            }
            new->QSORTED = self->QSORTED ;
	}
	Memory_Deallocate ( flags ) ;
    }
    return new ;
}

/*------------------------------------------------------------------------------
! . Sorting.
! . Within a harmonicangle, atom1 > atom3.
! . Within the array, ordering is done with increased values of atom2 and then
! . atom1 and then atom3.
! . Duplicate terms are not removed.
!-----------------------------------------------------------------------------*/
void HarmonicAngleContainer_Sort ( HarmonicAngleContainer *self )
{
    if ( ( self != NULL ) && ( ! self->QSORTED ) )
    {
        auto Integer atom1, atom3, i ;
        /* . Order atom1 and atom2 within each term. */
        for ( i = 0 ; i < self->nterms ; i++ )
        {
            atom1 = self->terms[i].atom1 ;
            atom3 = self->terms[i].atom3 ;
            if ( atom3 > atom1 )
            {
                self->terms[i].atom1 = atom3 ;
                self->terms[i].atom3 = atom1 ;
            }
        }
        /* . Order the terms within the container. */
        qsort ( ( void * ) self->terms, ( size_t ) self->nterms, sizeof ( HarmonicAngle ), ( void * ) HarmonicAngleTerm_Compare ) ;
        self->QSORTED = True ;
    }
}

/*------------------------------------------------------------------------------
! . Return the upper bound for the container.
! . This is the value of the largest index plus one.
!-----------------------------------------------------------------------------*/
int HarmonicAngleContainer_UpperBound ( HarmonicAngleContainer *self )
{
    Integer upperBound = 0 ;
    if ( ( self != NULL ) && ( self->nterms > 0 ) )
    {
        HarmonicAngleContainer_Sort ( self ) ;
        upperBound = Maximum ( self->terms[self->nterms-1].atom1, self->terms[self->nterms-1].atom2 ) + 1 ;
    }
    return upperBound ;
}

/*==============================================================================
! . Private procedures.
!============================================================================*/
static Integer HarmonicAngleTerm_Compare ( const void *vterm1, const void *vterm2 )
{
    HarmonicAngle *term1, *term2 ;
    Integer i ;
    term1 = ( HarmonicAngle * ) vterm1 ;
    term2 = ( HarmonicAngle * ) vterm2 ;
         if ( term1->atom2 < term2->atom2 ) i = -1 ;
    else if ( term1->atom2 > term2->atom2 ) i =  1 ;
    else if ( term1->atom1 < term2->atom1 ) i = -1 ;
    else if ( term1->atom1 > term2->atom1 ) i =  1 ;
    else if ( term1->atom3 < term2->atom3 ) i = -1 ;
    else if ( term1->atom3 > term2->atom3 ) i =  1 ;
    else if ( term1->type  < term2->type  ) i = -1 ;
    else if ( term1->type  > term2->type  ) i =  1 ;
    else if ( ! term1->QACTIVE && term2->QACTIVE ) i = -1 ;
    else if ( term1->QACTIVE && ! term2->QACTIVE ) i =  1 ;
    else i = 0 ;
    return i ;
}
