/*------------------------------------------------------------------------------
! . File      : HarmonicBondContainer.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==============================================================================
!=============================================================================*/

# include <math.h>
# include <stdlib.h>

# include "HarmonicBondContainer.h"
# include "Memory.h"

/*------------------------------------------------------------------------------
! . Local procedures.
!-----------------------------------------------------------------------------*/
static Integer HarmonicBondTerm_Compare ( const void *vterm1, const void *vterm2 ) ;

/*==============================================================================
! . Standard procedures.
!=============================================================================*/
/*------------------------------------------------------------------------------
! . Activate terms.
!-----------------------------------------------------------------------------*/
void HarmonicBondContainer_ActivateTerms ( HarmonicBondContainer *self )
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
HarmonicBondContainer *HarmonicBondContainer_Allocate ( const Integer nterms, const Integer nparameters )
{
    HarmonicBondContainer *self = NULL ;
    if ( ( nterms != 0 ) && ( nparameters != 0 ) )
    {
        Integer i ;
	self = ( HarmonicBondContainer * ) Memory_Allocate ( sizeof ( HarmonicBondContainer ) ) ;
        self->QSORTED     = False       ;
	self->nterms      = nterms      ;
	self->nparameters = nparameters ;
	self->terms	  = ( HarmonicBond * )          Memory_Allocate_Array ( nterms,      sizeof ( HarmonicBond )          ) ;
	self->parameters  = ( HarmonicBondParameter * ) Memory_Allocate_Array ( nparameters, sizeof ( HarmonicBondParameter ) ) ;
	/* . Make all terms inactive. */
	for ( i = 0 ; i < nterms ; i++ ) self->terms[i].QACTIVE = False ;
    }
    return self ;
}

/*------------------------------------------------------------------------------
! . Cloning.
!-----------------------------------------------------------------------------*/
HarmonicBondContainer *HarmonicBondContainer_Clone ( const HarmonicBondContainer *self )
{
    HarmonicBondContainer *new = NULL ;
    if ( self != NULL )
    {
        auto Integer i ;
        new = HarmonicBondContainer_Allocate ( self->nterms, self->nparameters ) ;
        for ( i = 0 ; i < self->nterms ; i++ )
        {
            new->terms[i].QACTIVE = self->terms[i].QACTIVE ;
            new->terms[i].atom1   = self->terms[i].atom1   ;
            new->terms[i].atom2   = self->terms[i].atom2   ;
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
void HarmonicBondContainer_DeactivateFixedAtomTerms ( HarmonicBondContainer *self, Selection *fixedatoms )
{
    if ( ( self != NULL ) && ( fixedatoms != NULL ) )
    {
        auto Integer i, n ;
        n = HarmonicBondContainer_UpperBound ( self ) ;
        Selection_MakeFlags ( fixedatoms, n ) ;
	for ( i = 0 ; i < self->nterms ; i++ )
	{
            if ( self->terms[i].QACTIVE )
            {
                self->terms[i].QACTIVE = ! ( fixedatoms->flags[self->terms[i].atom1] && fixedatoms->flags[self->terms[i].atom2] ) ;
            }
	}
    }
}

/*------------------------------------------------------------------------------
! . Deactivate terms involving QC atoms.
! . qcAtoms is the selection of both pure and boundary QC atoms.
! . Already deactivated terms are not affected.
!-----------------------------------------------------------------------------*/
void HarmonicBondContainer_DeactivateQCAtomTerms ( HarmonicBondContainer *self, Selection *qcAtoms, Selection *boundaryatoms )
{
    if ( ( self != NULL ) && ( qcAtoms != NULL ) )
    {
        auto Boolean QBOUNDARY, QEXCLUDE ;
        auto Integer  i, n ;
        n = HarmonicBondContainer_UpperBound ( self ) ;
        Selection_MakeFlags ( qcAtoms,       n ) ;
        Selection_MakeFlags ( boundaryatoms, n ) ;
        QBOUNDARY = ( boundaryatoms != NULL ) ;
	for ( i = 0 ; i < self->nterms ; i++ )
	{
            if ( self->terms[i].QACTIVE )
            {
                QEXCLUDE = ( qcAtoms->flags[self->terms[i].atom1] && qcAtoms->flags[self->terms[i].atom2] ) ;
                if ( QBOUNDARY && QEXCLUDE ) QEXCLUDE = ! ( boundaryatoms->flags[self->terms[i].atom1] || boundaryatoms->flags[self->terms[i].atom2] ) ;
                self->terms[i].QACTIVE = ! QEXCLUDE ;
            }
	}
    }
}

/*------------------------------------------------------------------------------
! . Deallocation.
!-----------------------------------------------------------------------------*/
void HarmonicBondContainer_Deallocate ( HarmonicBondContainer **self )
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
double HarmonicBondContainer_Energy ( const HarmonicBondContainer *self, const Coordinates3 *coordinates3, Coordinates3 *gradients3 )
{
    Real energy = 0.0e+00 ;
    if ( ( self != NULL ) && ( coordinates3 != NULL ) )
    {
	auto Boolean   QGRADIENTS ;
	auto Real df, disp, rij, xij, yij, zij ;
	auto Integer    i, j, n, t ;
	QGRADIENTS = ( gradients3 != NULL ) ;
	for ( n = 0 ; n < self->nterms ; n++ )
	{
	    if ( self->terms[n].QACTIVE )
	    {
		i = self->terms[n].atom1 ;
		j = self->terms[n].atom2 ;
		t = self->terms[n].type  ;
		Coordinates3_DifferenceRow ( coordinates3, i, j, xij, yij, zij ) ;
		rij  = sqrt ( xij * xij + yij * yij + zij * zij ) ;
		disp = rij - self->parameters[t].eq ;
		df   = self->parameters[t].fc * disp ;
		energy += ( df * disp ) ;
        	if ( QGRADIENTS )
        	{
        	    df  *= ( 2.0e+00 / rij ) ;
		    xij *= df ;
		    yij *= df ;
		    zij *= df ;
		    Coordinates3_IncrementRow ( gradients3, i, xij, yij, zij ) ;
		    Coordinates3_DecrementRow ( gradients3, j, xij, yij, zij ) ;
        	}
	    }
	}
    }
    return energy ;
}

/*------------------------------------------------------------------------------
! . Identify boundary atoms.
! . qcAtoms is the selection of pure QC atoms only.
! . No checking/sorting of the returned data need be done.
!-----------------------------------------------------------------------------*/
int HarmonicBondContainer_IdentifyBoundaryAtoms ( HarmonicBondContainer *self, Selection *qcAtoms, Integer **mmboundary, Integer **qcpartners )
{
    Integer n = 0 ;
    if ( ( self != NULL ) && ( qcAtoms != NULL ) )
    {
        auto Integer i ;
        n = HarmonicBondContainer_UpperBound ( self ) ;
        Selection_MakeFlags ( qcAtoms, n ) ;
        /* . Loop to find the number of boundary atoms. */
	for ( i = 0 ; i < self->nterms ; i++ )
	{
            if ( ( qcAtoms->flags[self->terms[i].atom1] && ! qcAtoms->flags[self->terms[i].atom2] ) ||
                 ( ! qcAtoms->flags[self->terms[i].atom1] && qcAtoms->flags[self->terms[i].atom2] ) ) n++ ;
	}
        /* . Fill the boundary atom arrays. */
        if ( ( n > 0 ) && ( mmboundary != NULL ) && ( qcpartners != NULL ) )
        {
            (*mmboundary) = Memory_Allocate_Array_Integer ( n ) ;
            (*qcpartners) = Memory_Allocate_Array_Integer ( n ) ;
            for ( i = 0, n = 0 ; i < self->nterms ; i++ )
	    {
                if ( qcAtoms->flags[self->terms[i].atom1] && ! qcAtoms->flags[self->terms[i].atom2] )
                {
                    (*qcpartners)[n] = self->terms[i].atom1 ;
                    (*mmboundary)[n] = self->terms[i].atom2 ;
                    n++ ;
                }
                else if ( ! qcAtoms->flags[self->terms[i].atom1] && qcAtoms->flags[self->terms[i].atom2] )
                {
                    (*mmboundary)[n] = self->terms[i].atom1 ;
                    (*qcpartners)[n] = self->terms[i].atom2 ;
                    n++ ;
                }
	    }
        }
    }
    return n ;
}

/*------------------------------------------------------------------------------
! . Merging.
!-----------------------------------------------------------------------------*/
HarmonicBondContainer *HarmonicBondContainer_Merge ( const HarmonicBondContainer *self, const HarmonicBondContainer *other, const Integer atomincrement )
{
    HarmonicBondContainer *new = NULL ;
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        auto Integer i ;
        new = HarmonicBondContainer_Allocate ( self->nterms + other->nterms, self->nparameters + other->nparameters ) ;
        for ( i = 0 ; i < self->nterms ; i++ )
        {
            new->terms[i].QACTIVE = self->terms[i].QACTIVE ;
            new->terms[i].atom1   = self->terms[i].atom1   ;
            new->terms[i].atom2   = self->terms[i].atom2   ;
            new->terms[i].type    = self->terms[i].type    ;
        }
        for ( i = 0 ; i < other->nterms ; i++ )
        {
            new->terms[i+self->nterms].QACTIVE = other->terms[i].QACTIVE ;
            new->terms[i+self->nterms].atom1   = other->terms[i].atom1 + atomincrement     ;
            new->terms[i+self->nterms].atom2   = other->terms[i].atom2 + atomincrement     ;
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
int HarmonicBondContainer_NumberOfInactiveTerms ( const HarmonicBondContainer *self )
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
HarmonicBondContainer *HarmonicBondContainer_Prune ( HarmonicBondContainer *self, Selection *selection )
{
    HarmonicBondContainer *new = NULL ;
    if ( ( self != NULL ) && ( selection != NULL ) )
    {
        auto Boolean *flags ;
        auto Integer   i, n  ;
        n = HarmonicBondContainer_UpperBound ( self ) ;
        Selection_MakeFlags     ( selection, n ) ;
        Selection_MakePositions ( selection, n ) ;
        flags = Memory_Allocate_Array_Boolean ( self->nterms ) ;
	for ( i = 0, n = 0 ; i < self->nterms ; i++ )
	{
            flags[i] = ( selection->flags[self->terms[i].atom1] && selection->flags[self->terms[i].atom2] ) ;
            if ( flags[i] ) n++ ;
	}
	if ( n > 0 )
	{
            new = HarmonicBondContainer_Allocate ( n, self->nparameters ) ;
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
! . Within a harmonicbond, atom1 > atom3.
! . Within the array, ordering is done with increased values of atom2 and then
! . atom1 and then atom3.
! . Duplicate terms are not removed.
!-----------------------------------------------------------------------------*/
void HarmonicBondContainer_Sort ( HarmonicBondContainer *self )
{
    if ( ( self != NULL ) && ( ! self->QSORTED ) )
    {
        auto Integer atom1, atom2, i ;
        /* . Order atom1 and atom2 within each term. */
        for ( i = 0 ; i < self->nterms ; i++ )
        {
            atom1 = self->terms[i].atom1 ;
            atom2 = self->terms[i].atom2 ;
            if ( atom2 > atom1 )
            {
                self->terms[i].atom1 = atom2 ;
                self->terms[i].atom2 = atom1 ;
            }
        }
        /* . Order the terms within the container. */
        qsort ( ( void * ) self->terms, ( size_t ) self->nterms, sizeof ( HarmonicBond ), ( void * ) HarmonicBondTerm_Compare ) ;
        self->QSORTED = True ;
    }
}

/*------------------------------------------------------------------------------
! . Return the upper bound for the container.
! . This is the value of the largest index plus one.
!-----------------------------------------------------------------------------*/
int HarmonicBondContainer_UpperBound ( HarmonicBondContainer *self )
{
    Integer upperBound = 0 ;
    if ( ( self != NULL ) && ( self->nterms > 0 ) )
    {
        HarmonicBondContainer_Sort ( self ) ;
        upperBound = self->terms[self->nterms-1].atom1 + 1 ;
    }
    return upperBound ;
}

/*==============================================================================
! . Private procedures.
!============================================================================*/
static Integer HarmonicBondTerm_Compare ( const void *vterm1, const void *vterm2 )
{
    HarmonicBond *term1, *term2 ;
    Integer i ;
    term1 = ( HarmonicBond * ) vterm1 ;
    term2 = ( HarmonicBond * ) vterm2 ;
         if ( term1->atom1 < term2->atom1 ) i = -1 ;
    else if ( term1->atom1 > term2->atom1 ) i =  1 ;
    else if ( term1->atom2 < term2->atom2 ) i = -1 ;
    else if ( term1->atom2 > term2->atom2 ) i =  1 ;
    else if ( term1->type  < term2->type  ) i = -1 ;
    else if ( term1->type  > term2->type  ) i =  1 ;
    else if ( ! term1->QACTIVE && term2->QACTIVE ) i = -1 ;
    else if ( term1->QACTIVE && ! term2->QACTIVE ) i =  1 ;
    else i = 0 ;
    return i ;
}

