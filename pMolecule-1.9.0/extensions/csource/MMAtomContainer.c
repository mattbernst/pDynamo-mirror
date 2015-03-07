/*------------------------------------------------------------------------------
! . File      : MMAtomContainer.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==============================================================================
!=============================================================================*/

# include "MMAtomContainer.h"
# include "Memory.h"
# include "Units.h"

/*==============================================================================
! . Standard procedures.
!=============================================================================*/
/*------------------------------------------------------------------------------
! . Allocation.
!-----------------------------------------------------------------------------*/
MMAtomContainer *MMAtomContainer_Allocate ( const Integer natoms )
{
    MMAtomContainer *self = NULL ;
    if ( natoms >= 0 )
    {
        auto Integer i ;
        self         = ( MMAtomContainer * ) Memory_Allocate ( sizeof ( MMAtomContainer ) ) ;
	self->data   = ( MMAtom * ) Memory_Allocate_Array ( natoms, sizeof ( MMAtom ) ) ;
	self->natoms = natoms ;
	for ( i = 0 ; i < natoms ; i++ )
	{
            self->data[i].QACTIVE  = False   ;
            self->data[i].charge   = 0.0e+00 ;
	    self->data[i].atomtype = -1 ;
	    self->data[i].ljtype   = -1 ;
	}
    }
    return self ;
}

/*------------------------------------------------------------------------------
! . Cloning.
!-----------------------------------------------------------------------------*/
MMAtomContainer *MMAtomContainer_Clone ( const MMAtomContainer *self )
{
    MMAtomContainer *new = NULL ;
    if ( self != NULL )
    {
        auto Integer i ;
        new = MMAtomContainer_Allocate ( self->natoms ) ;
        for ( i = 0 ; i < self->natoms ; i++ )
        {
            new->data[i].QACTIVE  = self->data[i].QACTIVE  ;
            new->data[i].charge   = self->data[i].charge   ;
            new->data[i].atomtype = self->data[i].atomtype ;
            new->data[i].ljtype   = self->data[i].ljtype   ;
        }
    }
    return new ;
}

/*------------------------------------------------------------------------------
! . Deallocation.
!-----------------------------------------------------------------------------*/
void MMAtomContainer_Deallocate ( MMAtomContainer **self )
{
    if ( (*self) != NULL )
    {
        free ( (*self)->data ) ;
        free ( (*self) ) ;
        (*self) = NULL   ;
    }
}

/*------------------------------------------------------------------------------
! . Dipole moment.
!-----------------------------------------------------------------------------*/
# define CONVERSION_FACTOR ( UNITS_LENGTH_ANGSTROMS_TO_BOHRS * UNITS_DIPOLE_ATOMIC_UNITS_TO_DEBYES )

Vector3 *MMAtomContainer_DipoleMoment ( const MMAtomContainer *self, const Coordinates3 *coordinates3, const Vector3 *center )
{
    Vector3 *dipole = NULL;
    if ( ( self != NULL ) && ( coordinates3 != NULL ) )
    {
        auto Real dx, dy, dz, q, x, xc = 0.0e+00, y, yc = 0.0e+00, z, zc = 0.0e+00 ;
        auto Integer    i ;
        /* . Initialization. */
        dipole = Vector3_Allocate ( )   ;
        Vector3_Set ( dipole, 0.0e+00 ) ;
        /* . Get the center. */
        if ( center != NULL )
        {
            xc = Vector3_Item ( center, 0 ) ;
            yc = Vector3_Item ( center, 1 ) ;
            zc = Vector3_Item ( center, 2 ) ;
        }
        /* . Loop over active atoms. */
        for ( i = 0, dx = dy = dz = 0.0e+00 ; i < self->natoms ; i++ )
        {
            if ( self->data[i].QACTIVE )
            {
                /* . Get some information about the atom. */
                q = self->data[i].charge ;
                Coordinates3_GetRow ( coordinates3, i, x, y, z ) ;
                /* . Charge part. */
                dx += q * ( x - xc ) ;
                dy += q * ( y - yc ) ;
                dz += q * ( z - zc ) ;
            }
        }
        /* . Convert to Debyes. */
        dx *= CONVERSION_FACTOR ;
        dy *= CONVERSION_FACTOR ;
        dz *= CONVERSION_FACTOR ;
        /* . Fill the vector. */
        Vector3_Item ( dipole, 0 ) = dx ;
        Vector3_Item ( dipole, 1 ) = dy ;
        Vector3_Item ( dipole, 2 ) = dz ;
    }
    return dipole ;
}

# undef CONVERSION_FACTOR

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get the MM atom charges.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real1DArray *MMAtomContainer_GetCharges ( const MMAtomContainer *self, const Boolean activeOnly, Status *status )
{
    Real1DArray *charges = NULL ;
    if ( self != NULL )
    {
        charges = Real1DArray_Allocate ( self->natoms, status ) ;
        if ( charges != NULL )
        {
            auto Integer i ;
            if ( activeOnly )
            {
                for ( i = 0 ; i < self->natoms ; i++ )
                {
                    if ( self->data[i].QACTIVE ) Real1DArray_Item ( charges, i ) = self->data[i].charge ;
                    else                         Real1DArray_Item ( charges, i ) = 0.0e+00 ;
                }
            }
            else
            {
                for ( i = 0 ; i < self->natoms ; i++ ) Real1DArray_Item ( charges, i ) = self->data[i].charge ;
            }
        }
    }
    return charges ;
}

/*------------------------------------------------------------------------------
! . Merging.
!-----------------------------------------------------------------------------*/
MMAtomContainer *MMAtomContainer_Merge ( const MMAtomContainer *self, const MMAtomContainer *other,
                                                                                      const Integer atomtypeincrement,
									              const Integer   ljtypeincrement )
{
    MMAtomContainer *new = NULL ;
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        auto Integer i ;
        new = MMAtomContainer_Allocate ( self->natoms + other->natoms ) ;
        for ( i = 0 ; i < self->natoms ; i++ )
        {
            new->data[i].QACTIVE  = self->data[i].QACTIVE    ;
            new->data[i].charge   = self->data[i].charge   ;
            new->data[i].atomtype = self->data[i].atomtype ;
            new->data[i].ljtype   = self->data[i].ljtype   ;
        }
        for ( i = 0 ; i < other->natoms ; i++ )
        {
            new->data[i+self->natoms].QACTIVE  = other->data[i].QACTIVE ;
            new->data[i+self->natoms].charge   = other->data[i].charge  ;
            new->data[i+self->natoms].atomtype = other->data[i].atomtype + atomtypeincrement ;
            new->data[i+self->natoms].ljtype   = other->data[i].ljtype   +   ljtypeincrement ;
        }
    }
    return new ;
}

/*------------------------------------------------------------------------------
! . Number of active atoms.
!-----------------------------------------------------------------------------*/
Integer MMAtomContainer_NumberOfActiveAtoms ( const MMAtomContainer *self )
{
    Integer numberActive = 0 ;
    if ( self != NULL )
    {
        auto Integer i ;
	for ( i = 0 ; i < self->natoms ; i++ )
	{
	    if ( self->data[i].QACTIVE ) numberActive += 1 ;
	}
    }
    return numberActive ;
}

/*------------------------------------------------------------------------------
! . Pruning.
!-----------------------------------------------------------------------------*/
MMAtomContainer *MMAtomContainer_Prune ( const MMAtomContainer *self, Selection *selection )
{
    MMAtomContainer *new = NULL ;
    if ( ( self != NULL ) && ( selection != NULL ) )
    {
        auto Integer i, n ;
        Selection_MakeFlags ( selection, self->natoms ) ;
        new = MMAtomContainer_Allocate ( selection->nindices ) ;
        if ( new->natoms > 0 )
        {
            for ( i = 0, n = 0 ; i < self->natoms ; i++ )
            {
                if ( selection->flags[i] )
                {
        	    new->data[n].QACTIVE  = self->data[i].QACTIVE  ;
        	    new->data[n].charge   = self->data[i].charge   ;
        	    new->data[n].atomtype = self->data[i].atomtype ;
        	    new->data[n].ljtype   = self->data[i].ljtype   ;
                    n++ ;
	        }
            }
        }
    }
    return new ;
}

/*------------------------------------------------------------------------------
! . Size.
!-----------------------------------------------------------------------------*/
int MMAtomContainer_Size ( const MMAtomContainer *self )
{
    Integer size = 0 ;
    if ( self != NULL ) size = self->natoms ;
    return size ;
}

/*------------------------------------------------------------------------------
! . Get the total charge of the active MM atoms.
!-----------------------------------------------------------------------------*/
Real MMAtomContainer_TotalCharge ( const MMAtomContainer *self )
{
    Real q = 0.0e+00 ;
    if ( self != NULL )
    {
        auto Integer i ;
        for ( i = 0 ; i < self->natoms ; i++ )
        {
	    if ( self->data[i].QACTIVE ) q += self->data[i].charge ;
	}
    }
    return q ;
}
