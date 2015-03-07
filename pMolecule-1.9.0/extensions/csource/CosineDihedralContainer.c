/*------------------------------------------------------------------------------
! . File      : CosineDihedralContainer.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>

# include "CosineDihedralContainer.h"
# include "Memory.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Real Binomial  ( const Integer n, const Integer k ) ;
static Real Factorial ( const Integer n ) ;

static void CosineDihedralParameter_Clone      ( CosineDihedralParameter *self, const CosineDihedralParameter *other ) ;
static void CosineDihedralParameter_Deallocate ( CosineDihedralParameter *self ) ;
static void CosineDihedralParameter_Initialize ( CosineDihedralParameter *self ) ;
static void CosineDihedralParameter_MakePowers ( CosineDihedralParameter *self ) ;
static Integer  CosineDihedralTerm_Compare         ( const void *vterm1, const void *vterm2 ) ;

/*==================================================================================================================================
! . Procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Activate terms.
!---------------------------------------------------------------------------------------------------------------------------------*/
void CosineDihedralContainer_ActivateTerms ( CosineDihedralContainer *self )
{
    if ( self != NULL )
    {
        auto Integer i ;
	for ( i = 0 ; i < self->nterms ; i++ ) self->terms[i].QACTIVE = True ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
CosineDihedralContainer *CosineDihedralContainer_Allocate ( const Integer nterms, const Integer nparameters )
{
    CosineDihedralContainer *self = NULL ;
    if ( ( nterms != 0 ) && ( nparameters != 0 ) )
    {
        Integer i ;
	self = ( CosineDihedralContainer * ) Memory_Allocate ( sizeof ( CosineDihedralContainer ) ) ;
        self->QSORTED     = False       ;
        self->nperiods    = 0           ;
 	self->nterms      = nterms      ;
	self->nparameters = nparameters ;
	self->terms	  = ( CosineDihedral * )          Memory_Allocate_Array ( nterms,      sizeof ( CosineDihedral )          ) ;
	self->parameters  = ( CosineDihedralParameter * ) Memory_Allocate_Array ( nparameters, sizeof ( CosineDihedralParameter ) ) ;
        /* . Initialize all parameters. */
        for ( i = 0 ; i < nparameters ; i++ ) CosineDihedralParameter_Initialize ( &(self->parameters[i]) ) ;
	/* . Make all terms inactive. */
	for ( i = 0 ; i < nterms ; i++ ) self->terms[i].QACTIVE = False ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
CosineDihedralContainer *CosineDihedralContainer_Clone ( const CosineDihedralContainer *self )
{
    CosineDihedralContainer *new = NULL ;
    if ( self != NULL )
    {
        auto Integer i ;
        new = CosineDihedralContainer_Allocate ( self->nterms, self->nparameters ) ;
        for ( i = 0 ; i < self->nterms ; i++ )
        {
            new->terms[i].QACTIVE = self->terms[i].QACTIVE ;
            new->terms[i].atom1   = self->terms[i].atom1   ;
            new->terms[i].atom2   = self->terms[i].atom2   ;
            new->terms[i].atom3   = self->terms[i].atom3   ;
            new->terms[i].atom4   = self->terms[i].atom4   ;
            new->terms[i].type    = self->terms[i].type    ;
        }
        for ( i = 0 ; i < self->nparameters ; i++ ) CosineDihedralParameter_Clone ( &(new->parameters[i]), &(self->parameters[i]) ) ;
        new->nperiods = self->nperiods ;
        new->QSORTED  = self->QSORTED  ;
    }
    return new ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deactivate terms between fixed atoms.
! . Already deactivated terms are not affected.
!---------------------------------------------------------------------------------------------------------------------------------*/
void CosineDihedralContainer_DeactivateFixedAtomTerms ( CosineDihedralContainer *self, Selection *fixedatoms )
{
    if ( ( self != NULL ) && ( fixedatoms != NULL ) )
    {
        auto Integer i, n ;
        n = CosineDihedralContainer_UpperBound ( self ) ;
        Selection_MakeFlags ( fixedatoms, n ) ;
	for ( i = 0 ; i < self->nterms ; i++ )
	{
            if ( self->terms[i].QACTIVE )
            {
                self->terms[i].QACTIVE = ! ( fixedatoms->flags[self->terms[i].atom1] &&
                                             fixedatoms->flags[self->terms[i].atom2] &&
                                             fixedatoms->flags[self->terms[i].atom3] &&
                                             fixedatoms->flags[self->terms[i].atom4] ) ;
            }
	}
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deactivate terms involving QC atoms.
! . qcAtoms is the selection of both pure and boundary QC atoms.
! . Already deactivated terms are not affected.
!---------------------------------------------------------------------------------------------------------------------------------*/
void CosineDihedralContainer_DeactivateQCAtomTerms ( CosineDihedralContainer *self, Selection *qcAtoms, Selection *boundaryatoms )
{
    if ( ( self != NULL ) && ( qcAtoms != NULL ) )
    {
        auto Boolean QEXCLUDE ;
        auto Integer  i, n ;
        n = CosineDihedralContainer_UpperBound ( self ) ;
        Selection_MakeFlags ( qcAtoms, n ) ;
	for ( i = 0 ; i < self->nterms ; i++ )
	{
            if ( self->terms[i].QACTIVE )
            {
                QEXCLUDE = ( qcAtoms->flags[self->terms[i].atom1] && qcAtoms->flags[self->terms[i].atom2] && qcAtoms->flags[self->terms[i].atom3] && qcAtoms->flags[self->terms[i].atom4] ) ;
                self->terms[i].QACTIVE = ! QEXCLUDE ;
            }
	}
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void CosineDihedralContainer_Deallocate ( CosineDihedralContainer **self )
{
    if ( (*self) != NULL )
    {
        auto Integer i ;
        for ( i = 0 ; i < (*self)->nparameters ; i++ ) CosineDihedralParameter_Deallocate ( &((*self)->parameters[i]) ) ;
        free ( (*self)->terms      ) ;
        free ( (*self)->parameters ) ;
        free ( (*self) ) ;
        (*self) = NULL   ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Energy and gradients.
! . Cosine terms only.
!---------------------------------------------------------------------------------------------------------------------------------*/
double CosineDihedralContainer_Energy ( const CosineDihedralContainer *self, const Coordinates3 *coordinates3, Coordinates3 *gradients3 )
{
    Real energy = 0.0e+00 ;
    if ( ( self != NULL ) && ( coordinates3 != NULL ) )
    {
	auto Boolean   QGRADIENTS ;
        auto Real c, cn, co, cosphi, df, dtxi, dtyi, dtzi, dtxj, dtyj, dtzj, dtxk, dtyk, dtzk, dtxl, dtyl, dtzl, e, m, mx, my, mz, n, nx, ny, nz, sx, sy, sz, xij, yij, zij, xkj, ykj, zkj, xlk, ylk, zlk ;
        auto Integer    i, j, k, l, nt, p, t ;
	QGRADIENTS = ( gradients3 != NULL ) ;
	for ( nt = 0 ; nt < self->nterms ; nt++ )
	{
	    if ( self->terms[nt].QACTIVE )
	    {
                /* . Local data. */
	        i = self->terms[nt].atom1 ;
	        j = self->terms[nt].atom2 ;
	        k = self->terms[nt].atom3 ;
	        l = self->terms[nt].atom4 ;
	        t = self->terms[nt].type  ;
	        /* . Coordinate displacements. */
	        Coordinates3_DifferenceRow ( coordinates3, i, j, xij, yij, zij ) ;
	        Coordinates3_DifferenceRow ( coordinates3, k, j, xkj, ykj, zkj ) ;
	        Coordinates3_DifferenceRow ( coordinates3, l, k, xlk, ylk, zlk ) ;
                /* . m and n. */
                mx = yij * zkj - zij * ykj ;
	        my = zij * xkj - xij * zkj ;
	        mz = xij * ykj - yij * xkj ;
	        nx = ylk * zkj - zlk * ykj ;
	        ny = zlk * xkj - xlk * zkj ;
	        nz = xlk * ykj - ylk * xkj ;
                /* . Normalize. */
                m = sqrt ( mx * mx + my * my + mz * mz ) ;
                mx /= m ; my /= m ; mz /= m ;
	        n = sqrt ( nx * nx + ny * ny + nz * nz ) ;
	        nx /= n ; ny /= n ; nz /= n ;
                /* . Cosine of the dihedral. */
                cosphi = mx * nx +  my * ny +  mz * nz ;
                /* . Loop over powers of the cosine. */
                cn = 1.0e+00 ;
                co = 0.0e+00 ;
                df = 0.0e+00 ;
                e  = 0.0e+00 ;
                for ( p = 0 ; p <= self->parameters[t].npowers ; p++ )
                {
                    c   = self->parameters[t].powercoefficients[p] ;
                    df += c * co * ( Real ) p ;
                    e  += c * cn ;
                    co  = cn ;
                    cn *= cosphi ;
                }
                /* . The energy term. */
                energy += e ;
                if ( QGRADIENTS )
                {
                    /* . rkj ^ m terms. */
                    sx = ykj * mz - zkj * my ;
	            sy = zkj * mx - xkj * mz ;
	            sz = xkj * my - ykj * mx ;
                    dtxi = - cosphi * sx ;
                    dtyi = - cosphi * sy ;
                    dtzi = - cosphi * sz ;
                    dtxl = sx ;
                    dtyl = sy ;
                    dtzl = sz ;
                    /* . rkj ^ n terms. */
                    sx = ykj * nz - zkj * ny ;
	            sy = zkj * nx - xkj * nz ;
	            sz = xkj * ny - ykj * nx ;
                    dtxi += sx ;
                    dtyi += sy ;
                    dtzi += sz ;
                    dtxl -= cosphi * sx ;
                    dtyl -= cosphi * sy ;
                    dtzl -= cosphi * sz ;
                    /* . Finish i and l. */
                    dtxi *= df / m ;
                    dtyi *= df / m ;
                    dtzi *= df / m ;
                    dtxl *= df / n ;
                    dtyl *= df / n ;
                    dtzl *= df / n ;
                    /* . Scale rij. */
                    xij /= m ; yij /= m ; zij /= m ;
                    /* . rij ^ m terms. */
                    sx = yij * mz - zij * my ;
	            sy = zij * mx - xij * mz ;
	            sz = xij * my - yij * mx ;
                    dtxk = cosphi * sx ;
                    dtyk = cosphi * sy ;
                    dtzk = cosphi * sz ;
                    /* . rij ^ n terms. */
                    sx = yij * nz - zij * ny ;
	            sy = zij * nx - xij * nz ;
	            sz = xij * ny - yij * nx ;
                    dtxk -= sx ;
                    dtyk -= sy ;
                    dtzk -= sz ;
                    /* . Scale rlk. */
                    xlk /= n ; ylk /= n ; zlk /= n ;
                    /* . rlk ^ m terms. */
                    sx = ylk * mz - zlk * my ;
	            sy = zlk * mx - xlk * mz ;
	            sz = xlk * my - ylk * mx ;
                    dtxk -= sx ;
                    dtyk -= sy ;
                    dtzk -= sz ;
                    /* . rlk ^ n terms. */
                    sx = ylk * nz - zlk * ny ;
	            sy = zlk * nx - xlk * nz ;
	            sz = xlk * ny - ylk * nx ;
                    dtxk += cosphi * sx ;
                    dtyk += cosphi * sy ;
                    dtzk += cosphi * sz ;
                    /* . Scale k. */
                    dtxk *= df ; dtyk *= df ; dtzk *= df ;
                    /* . Finish j and k. */
	            dtxj  = - dtxk - dtxi ;
	            dtyj  = - dtyk - dtyi ;
	            dtzj  = - dtzk - dtzi ;
	            dtxk -= dtxl ;
	            dtyk -= dtyl ;
	            dtzk -= dtzl ;
                    /* . Add in the contributions. */
	            Coordinates3_IncrementRow ( gradients3, i, dtxi, dtyi, dtzi ) ;
	            Coordinates3_IncrementRow ( gradients3, j, dtxj, dtyj, dtzj ) ;
	            Coordinates3_IncrementRow ( gradients3, k, dtxk, dtyk, dtzk ) ;
	            Coordinates3_IncrementRow ( gradients3, l, dtxl, dtyl, dtzl ) ;
        	}
	    }
	}
    }
    return energy ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Find the maximum period.
!---------------------------------------------------------------------------------------------------------------------------------*/
void CosineDihedralContainer_FindMaximumPeriod ( CosineDihedralContainer *self )
{
    if ( self != NULL )
    {
        auto Integer i, n, p = 0 ;
        for ( i = 0 ; i < self->nparameters ; i++ )
        {
            for ( n = 0 ; n < self->parameters[i].nterms ; n++ )
            {
                p = Maximum ( p, self->parameters[i].periods[n] ) ;
            }
        }
        self->nperiods = p ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make the powers representation of the parameters.
!---------------------------------------------------------------------------------------------------------------------------------*/
void CosineDihedralContainer_MakePowers ( CosineDihedralContainer *self )
{
    if ( self != NULL )
    {
        auto Integer i ;
        for ( i = 0 ; i < self->nparameters ; i++ ) CosineDihedralParameter_MakePowers ( &(self->parameters[i]) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Merging.
!---------------------------------------------------------------------------------------------------------------------------------*/
CosineDihedralContainer *CosineDihedralContainer_Merge ( const CosineDihedralContainer *self, const CosineDihedralContainer *other, const Integer atomincrement )
{
    CosineDihedralContainer *new = NULL ;
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        auto Integer i ;
        new = CosineDihedralContainer_Allocate ( self->nterms + other->nterms, self->nparameters + other->nparameters ) ;
        for ( i = 0 ; i < self->nterms ; i++ )
        {
            new->terms[i].QACTIVE = self->terms[i].QACTIVE ;
            new->terms[i].atom1   = self->terms[i].atom1   ;
            new->terms[i].atom2   = self->terms[i].atom2   ;
            new->terms[i].atom3   = self->terms[i].atom3   ;
            new->terms[i].atom4   = self->terms[i].atom4   ;
            new->terms[i].type    = self->terms[i].type    ;
        }
        for ( i = 0 ; i < other->nterms ; i++ )
        {
            new->terms[i+self->nterms].QACTIVE = other->terms[i].QACTIVE ;
            new->terms[i+self->nterms].atom1   = other->terms[i].atom1 + atomincrement     ;
            new->terms[i+self->nterms].atom2   = other->terms[i].atom2 + atomincrement     ;
            new->terms[i+self->nterms].atom3   = other->terms[i].atom3 + atomincrement     ;
            new->terms[i+self->nterms].atom4   = other->terms[i].atom4 + atomincrement     ;
            new->terms[i+self->nterms].type    = other->terms[i].type  + self->nparameters ;
        }
        for ( i = 0 ; i < self->nparameters  ; i++ ) CosineDihedralParameter_Clone ( &(new->parameters[i]),                   &(self->parameters[i])                   ) ;
        for ( i = 0 ; i < other->nparameters ; i++ ) CosineDihedralParameter_Clone ( &(new->parameters[i+self->nparameters]), &(self->parameters[i+self->nparameters]) ) ;
        new->nperiods = Maximum ( self->nperiods, other->nperiods ) ;
        new->QSORTED  = ( self->QSORTED && other->QSORTED ) ;
    }
    return new ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Return the number of inactive terms.
!---------------------------------------------------------------------------------------------------------------------------------*/
int CosineDihedralContainer_NumberOfInactiveTerms ( const CosineDihedralContainer *self )
{
    Integer n = 0 ;
    if ( self != NULL )
    {
        auto Integer i ;
	for ( i = 0 ; i < self->nterms ; i++ ) if ( ! self->terms[i].QACTIVE ) n++ ;
    }
    return n ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Pruning.
! . Only terms are pruned.
!---------------------------------------------------------------------------------------------------------------------------------*/
CosineDihedralContainer *CosineDihedralContainer_Prune ( CosineDihedralContainer *self, Selection *selection )
{
    CosineDihedralContainer *new = NULL ;
    if ( ( self != NULL ) && ( selection != NULL ) )
    {
        auto Boolean *flags ;
        auto Integer   i, n  ;
        n = CosineDihedralContainer_UpperBound ( self ) ;
        Selection_MakeFlags     ( selection, n ) ;
        Selection_MakePositions ( selection, n ) ;
        flags = Memory_Allocate_Array_Boolean ( self->nterms ) ;
	for ( i = 0, n = 0 ; i < self->nterms ; i++ )
	{
            flags[i] = ( selection->flags[self->terms[i].atom1] && selection->flags[self->terms[i].atom2] && selection->flags[self->terms[i].atom3] && selection->flags[self->terms[i].atom4] ) ;
            if ( flags[i] ) n++ ;
	}
	if ( n > 0 )
	{
            new = CosineDihedralContainer_Allocate ( n, self->nparameters ) ;
            for ( i = 0 ; i < self->nparameters ; i++ ) CosineDihedralParameter_Clone ( &(new->parameters[i]), &(self->parameters[i]) ) ;
            for ( i = 0, n = 0 ; i < self->nterms ; i++ )
            {
                if ( flags[i] )
                {
        	    new->terms[n].QACTIVE =                      self->terms[i].QACTIVE ;
        	    new->terms[n].atom1   = selection->positions[self->terms[i].atom1]  ;
        	    new->terms[n].atom2   = selection->positions[self->terms[i].atom2]  ;
        	    new->terms[n].atom3   = selection->positions[self->terms[i].atom3]  ;
        	    new->terms[n].atom4   = selection->positions[self->terms[i].atom4]  ;
        	    new->terms[n].type    =                      self->terms[i].type    ;
        	    n++ ;
                }
            }
            new->nperiods = self->nperiods ;
            new->QSORTED  = self->QSORTED  ;
	}
	Memory_Deallocate ( flags ) ;
    }
    return new ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Sorting.
! . Within a cosinedihedral, atom2 > atom3.
! . Within the array, ordering is done with increased values of atom2 and then
! . atom3 and then atom1 and then atom4.
! . Duplicates are not removed.
!---------------------------------------------------------------------------------------------------------------------------------*/
void CosineDihedralContainer_Sort ( CosineDihedralContainer *self )
{
    if ( ( self != NULL ) && ( ! self->QSORTED ) )
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
        }
        /* . Order the terms within the container. */
        qsort ( ( void * ) self->terms, ( size_t ) self->nterms, sizeof ( CosineDihedral ), ( void * ) CosineDihedralTerm_Compare ) ;
        self->QSORTED = True ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Return the upper bound for the container.
! . This is the value of the largest index plus one.
!---------------------------------------------------------------------------------------------------------------------------------*/
int CosineDihedralContainer_UpperBound ( CosineDihedralContainer *self )
{
    Integer upperBound = 0 ;
    if ( ( self != NULL ) && ( self->nterms > 0 ) )
    {
        CosineDihedralContainer_Sort ( self ) ;
        /* . atom3 is always less than or equal to atom2. */
        upperBound  = Maximum ( self->terms[self->nterms-1].atom1, self->terms[self->nterms-1].atom2 ) ;
        upperBound  = Maximum ( upperBound, self->terms[self->nterms-1].atom4 ) ;
        upperBound += 1 ;
    }
    return upperBound ;
}

/*==================================================================================================================================
! . Private procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Parameter internal allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void CosineDihedralParameter_Allocate ( CosineDihedralParameter *self, const Integer nterms )
{
    if ( ( self != NULL ) && ( nterms > 0 ) )
    {
        self->nterms            = nterms ;
        self->coefficients      = Memory_Allocate_Array_Real_Initialize    ( nterms, 0.0e+00 ) ;
        self->periods           = Memory_Allocate_Array_Integer_Initialize ( nterms, 0       ) ;
        self->npowers           = -1 ;
        self->powercoefficients = NULL ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Parameter internal clone.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void CosineDihedralParameter_Clone ( CosineDihedralParameter *self, const CosineDihedralParameter *other )
{
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        auto Integer i ;
        CosineDihedralParameter_Allocate ( self, other->nterms ) ;
        for ( i = 0 ; i < other->nterms ; i++ )
        {
            self->coefficients[i] = other->coefficients[i] ;
            self->periods     [i] = other->periods     [i] ;
        }
        CosineDihedralParameter_MakePowers ( self ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Parameter internal deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void CosineDihedralParameter_Deallocate ( CosineDihedralParameter *self )
{
    if ( self != NULL )
    {
        self->npowers = -1 ;
        self->nterms  =  0 ;
        Memory_Deallocate_Real    ( &(self->coefficients     ) ) ;
        Memory_Deallocate_Real    ( &(self->powercoefficients) ) ;
        Memory_Deallocate_Integer ( &(self->periods          ) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Parameter internal initialization.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void CosineDihedralParameter_Initialize ( CosineDihedralParameter *self )
{
    if ( self != NULL )
    {
        self->npowers           = -1   ;
        self->nterms            =  0   ;
        self->coefficients      = NULL ;
        self->periods           = NULL ;
        self->powercoefficients = NULL ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Parameter internal make powers.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Real Binomial ( const Integer n, const Integer k )
{
    return Factorial ( n ) / ( Factorial ( n - k ) * Factorial ( k ) ) ;
}

static Real Factorial ( const Integer n )
{
    auto Real f = 1.0e+00 ;
    auto Integer    i = n ;
    while ( i > 0 )
    {
        f *= ( Real ) i ;
        i -= 1 ;
    }
    return f ;
}

static void CosineDihedralParameter_MakePowers ( CosineDihedralParameter *self )
{
    if ( ( self != NULL ) && ( self->npowers < 0 ) )
    {
        auto Integer i, npowers ;
        /* . Find the largest power. */
        npowers = -1 ;
        for ( i = 0 ; i < self->nterms ; i++ ) npowers = Maximum ( npowers, self->periods[i] ) ;
        /* . Determine the power coefficients. */
        if ( npowers >= 0 )
        {
            auto Real c, f, phase ;
            auto Integer    l, m, n, nby2 ;
            /* . Initialization. */
            self->npowers           = npowers ;
            self->powercoefficients = Memory_Allocate_Array_Real_Initialize  ( npowers + 1, 0.0e+00 ) ;
            /* . Loop over terms. */
            for ( i = 0 ; i < self->nterms ; i++ )
            {
                c     = self->coefficients[i] ;
                n     = self->periods[i]      ;
                nby2  = n / 2 ;
                phase = 1.0e+00 ;
                for ( l = 0 ; l <= nby2 ; l++ )
                {
                    f = 0.0e+00 ;
                    for ( m = l ; m <= nby2 ; m++ ) f += Binomial ( n, 2 * m ) * Binomial ( m, l ) ;
                    self->powercoefficients[n-2*l] += phase * f * c ;
/* printf ( "\nCoefficient for n = %d and power = %d = %10.3f\n", n, n-2*l, phase * f ) ; */
                    phase *= -1.0e+00 ;
                }
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Term comparison function.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Integer CosineDihedralTerm_Compare ( const void *vterm1, const void *vterm2 )
{
    CosineDihedral *term1, *term2 ;
    Integer i ;
    term1 = ( CosineDihedral * ) vterm1 ;
    term2 = ( CosineDihedral * ) vterm2 ;
         if ( term1->atom2 < term2->atom2 ) i = -1 ;
    else if ( term1->atom2 > term2->atom2 ) i =  1 ;
    else if ( term1->atom3 < term2->atom3 ) i = -1 ;
    else if ( term1->atom3 > term2->atom3 ) i =  1 ;
    else if ( term1->atom1 < term2->atom1 ) i = -1 ;
    else if ( term1->atom1 > term2->atom1 ) i =  1 ;
    else if ( term1->atom4 < term2->atom4 ) i = -1 ;
    else if ( term1->atom4 > term2->atom4 ) i =  1 ;
    else if ( term1->type  < term2->type  ) i = -1 ;
    else if ( term1->type  > term2->type  ) i =  1 ;
    else if ( ! term1->QACTIVE && term2->QACTIVE ) i = -1 ;
    else if ( term1->QACTIVE && ! term2->QACTIVE ) i =  1 ;
    else i = 0 ;
    return i ;
}
