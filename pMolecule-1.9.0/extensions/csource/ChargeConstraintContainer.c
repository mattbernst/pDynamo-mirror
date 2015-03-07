/*------------------------------------------------------------------------------
! . File      : ChargeConstraintContainer.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . A container for charge and spin constraints.
!=================================================================================================================================*/

# include "ChargeConstraintContainer.h"
# include "Memory.h"

/*==================================================================================================================================
! . Container procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
ChargeConstraintContainer *ChargeConstraintContainer_Allocate ( const Integer numberOfConstraints, Status *status )
{
    ChargeConstraintContainer *self = NULL ;
    MEMORY_ALLOCATE ( self, ChargeConstraintContainer ) ;
    if ( self != NULL )
    {
        auto Integer n ;
        n = Maximum ( numberOfConstraints, 0 ) ;
        /* . Initialization. */
        self->hasCharges          = False ;
        self->hasSpins            = False ;
        self->highestIndex        = -1    ;
        self->numberOfConstraints = n     ;
        self->constraints         = NULL  ;
        /* . Allocation. */ 
        if ( n > 0 )
        {
            /* . Allocation. */
            MEMORY_ALLOCATEARRAY ( self->constraints, n, ChargeConstraint * ) ;
            if ( self->constraints == NULL ) ChargeConstraintContainer_Deallocate ( &self ) ;
            else
            {
                auto Integer c ;
                for ( c = 0 ; c < n ; c++ ) self->constraints[c] = NULL ;
            }
        }
    }
    if ( self == NULL ) Status_Set ( status, Status_MemoryAllocationFailure ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
ChargeConstraintContainer *ChargeConstraintContainer_Clone ( const ChargeConstraintContainer *self, Status *status )
{
    ChargeConstraintContainer *clone = NULL ;
    if ( self != NULL )
    {
        clone = ChargeConstraintContainer_Allocate ( self->numberOfConstraints, status ) ;
        if ( clone != NULL )
        {
            auto Integer c ;
            /* . Basic data. */
            clone->hasCharges   = self->hasCharges   ;
            clone->hasSpins     = self->hasSpins     ;
            clone->highestIndex = self->highestIndex ;
            /* . Constraints. */
            for ( c = 0 ; c < self->numberOfConstraints ; c++ )
            {
                clone->constraints[c] = ChargeConstraint_Clone ( self->constraints[c], status ) ;
                if ( clone->constraints[c] == NULL ) { ChargeConstraintContainer_Deallocate ( &clone ) ; break ; }
            }
        }
        if ( clone == NULL ) Status_Set ( status, Status_MemoryAllocationFailure ) ;
    }
    return clone ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void ChargeConstraintContainer_Deallocate ( ChargeConstraintContainer **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) )
    {
        auto Integer c ;
        for ( c = 0 ; c < (*self)->numberOfConstraints ; c++ ) ChargeConstraint_Deallocate ( &((*self)->constraints[c]) ) ;
        MEMORY_DEALLOCATE ( (*self)->constraints ) ;
        MEMORY_DEALLOCATE ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Determine the deviations for each charge constraint given arrays of charges and spins.
!---------------------------------------------------------------------------------------------------------------------------------*/
void ChargeConstraintContainer_Deviations ( const ChargeConstraintContainer *self       ,
                                            const Real1DArray               *charges    ,
                                            const Real1DArray               *spins      ,
                                                  Real1DArray               *deviations ,
                                                  Status                    *status     )
{
    if ( ( self != NULL ) && ( deviations != NULL ) && ( deviations->length >= self->numberOfConstraints ) &&
         ( ( ! self->hasCharges ) || ( self->hasCharges && ( charges != NULL ) && ( charges->length > self->highestIndex ) ) ) &&
         ( ( ! self->hasSpins   ) || ( self->hasSpins   && ( spins   != NULL ) && ( spins->length   > self->highestIndex ) ) ) )
    {
        auto Integer           a, c, t    ;
        auto Real              value, w   ;
        auto ChargeConstraint *constraint ;

        /* . Initialization. */
        Real1DArray_Set ( deviations, 0.0e+00 ) ;

        /* . Loop over constraints. */
        for ( c = 0 ; c < self->numberOfConstraints ; c++ )
        {
            constraint =   self->constraints[c] ;
            value      = - constraint->target   ;
            for ( t = 0 ; t < constraint->numberOfCharges ; t++ )
            {
                a = Integer1DArray_Item ( constraint->chargeIndices, t ) ;
                w = Real1DArray_Item    ( constraint->chargeWeights, t ) ;
                value += ( w * Real1DArray_Item ( charges, a ) ) ;
            }
            for ( t = 0 ; t < constraint->numberOfSpins ; t++ )
            {
                a = Integer1DArray_Item ( constraint->spinIndices, t ) ;
                w = Real1DArray_Item    ( constraint->spinWeights, t ) ;
                value += ( w * Real1DArray_Item ( spins, a ) ) ;
            }
            Real1DArray_Item ( deviations, c ) = value ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . The higest index in the container.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer ChargeConstraintContainer_HighestIndex ( const ChargeConstraintContainer *self )
{
    Integer index = 0 ;
    if ( self != NULL ) index = self->highestIndex ;
    return index ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Pruning.
!---------------------------------------------------------------------------------------------------------------------------------*/
ChargeConstraintContainer *ChargeConstraintContainer_Prune ( const ChargeConstraintContainer *self, Selection *selection, Status *status )
{
    ChargeConstraintContainer *new = NULL ;
    if ( ( self != NULL ) && ( self->numberOfConstraints > 0 ) && ( selection != NULL ) )
    {
        auto Status localStatus = Status_Continue ;
        auto ChargeConstraint *constraint, **temporary ;
        MEMORY_ALLOCATEARRAY ( temporary, self->numberOfConstraints, ChargeConstraint * ) ;
        if ( temporary != NULL )
        {
            auto Integer c, n ;
            Selection_MakeFlags     ( selection, self->highestIndex+1 ) ;
            Selection_MakePositions ( selection, self->highestIndex+1 ) ;
            for ( c = n = 0 ; c < self->numberOfConstraints ; c++ )
            {
                constraint = ChargeConstraint_Prune ( self->constraints[c], selection, &localStatus ) ;
                if ( localStatus != Status_Continue ) { n = 0 ; break ; }
                else if ( constraint != NULL ) { temporary[n] = constraint ; n++ ; }
            }
            if ( n > 0 )
            {
                new = ChargeConstraintContainer_Allocate ( n, &localStatus ) ;
                if ( new != NULL ) { for ( c = 0 ; c < n ; c++ ) ChargeConstraintContainer_SetItem ( new, c, &(temporary[c]), NULL ) ; }
            }
            MEMORY_DEALLOCATE ( temporary ) ;
        }
        else localStatus = Status_OutOfMemory ;
        if ( localStatus != Status_Continue )
        {
            ChargeConstraintContainer_Deallocate ( &new ) ;
            Status_Set ( status, localStatus ) ;
        }
    }
    return new ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Add a constraint.
!---------------------------------------------------------------------------------------------------------------------------------*/
void ChargeConstraintContainer_SetItem ( ChargeConstraintContainer *self, const Integer index, ChargeConstraint **constraint, Status *status )
{
    if ( ( self != NULL ) && ( constraint != NULL ) && ( (*constraint) != NULL ) )
    {
        if ( ( index < 0 ) || ( index >= self->numberOfConstraints ) ) Status_Set ( status, Status_IndexOutOfRange ) ;
        else
        {
            /* . The item. */
            self->constraints[index] = (*constraint) ;
            /* . Other data. */
            self->hasCharges = self->hasCharges || ( (*constraint)->numberOfCharges > 0 ) ;
            self->hasSpins   = self->hasSpins   || ( (*constraint)->numberOfSpins   > 0 ) ;
            if ( (*constraint)->numberOfCharges > 0 ) self->highestIndex = Maximum ( self->highestIndex, Integer1DArray_Maximum ( (*constraint)->chargeIndices ) ) ;
            if ( (*constraint)->numberOfSpins   > 0 ) self->highestIndex = Maximum ( self->highestIndex, Integer1DArray_Maximum ( (*constraint)->spinIndices   ) ) ;
            /* . Take ownership. */
            (*constraint) = NULL ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . The size of the container.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer ChargeConstraintContainer_Size ( const ChargeConstraintContainer *self )
{
    Integer size = 0 ;
    if ( self != NULL ) size = self->numberOfConstraints ;
    return size ;
}

/*==================================================================================================================================
! . Constraint procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
ChargeConstraint *ChargeConstraint_Allocate ( const Integer numberOfCharges, const Integer numberOfSpins, Status *status )
{
    ChargeConstraint *self = NULL ;
    MEMORY_ALLOCATE ( self, ChargeConstraint ) ;
    if ( self != NULL )
    {
        auto Boolean isOK = True ;
        auto Integer q, s ;
        q = Maximum ( numberOfCharges, 0 ) ;
        s = Maximum ( numberOfSpins  , 0 ) ;
        self->numberOfCharges = q ;
        self->numberOfSpins   = s ;
        self->target          = 0.0e+00 ;
        self->chargeIndices   = NULL ;
        self->spinIndices     = NULL ;
        self->chargeWeights   = NULL ;
        self->spinWeights     = NULL ;
        if ( q > 0 )
        {
            self->chargeIndices = Integer1DArray_Allocate ( q, status ) ; Integer1DArray_Set ( self->chargeIndices,      -1 ) ;
            self->chargeWeights = Real1DArray_Allocate    ( q, status ) ; Real1DArray_Set    ( self->chargeWeights, 0.0e+00 ) ;
            isOK = isOK && ( self->chargeIndices != NULL ) && ( self->chargeWeights != NULL ) ;
        }
        if ( s > 0 )
        {
            self->spinIndices   = Integer1DArray_Allocate ( q, status ) ; Integer1DArray_Set ( self->spinIndices,      -1 ) ;
            self->spinWeights   = Real1DArray_Allocate    ( q, status ) ; Real1DArray_Set    ( self->spinWeights, 0.0e+00 ) ;
            isOK = isOK && ( self->spinIndices != NULL ) && ( self->spinWeights != NULL ) ;
        }
        if ( ! isOK ) ChargeConstraint_Deallocate ( &self ) ;
    }
    if ( self == NULL ) Status_Set ( status, Status_MemoryAllocationFailure ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
ChargeConstraint *ChargeConstraint_Clone ( const ChargeConstraint *self, Status *status )
{
    ChargeConstraint *clone = NULL ;
    if ( self != NULL )
    {
        clone = ChargeConstraint_Allocate ( self->numberOfCharges, self->numberOfSpins, status ) ;
        if ( clone != NULL )
        {
            clone->target = self->target ;
            Integer1DArray_CopyTo ( self->chargeIndices, clone->chargeIndices, status ) ;
            Integer1DArray_CopyTo ( self->spinIndices  , clone->spinIndices  , status ) ;
            Real1DArray_CopyTo    ( self->chargeWeights, clone->chargeWeights, status ) ;
            Real1DArray_CopyTo    ( self->spinWeights  , clone->spinWeights  , status ) ;
        }
        if ( clone == NULL ) Status_Set ( status, Status_MemoryAllocationFailure ) ;
    }
    return clone ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void ChargeConstraint_Deallocate ( ChargeConstraint **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) )
    {
        Integer1DArray_Deallocate ( &((*self)->chargeIndices) ) ;
        Integer1DArray_Deallocate ( &((*self)->spinIndices  ) ) ;
        Real1DArray_Deallocate    ( &((*self)->chargeWeights) ) ;
        Real1DArray_Deallocate    ( &((*self)->spinWeights  ) ) ;
        MEMORY_DEALLOCATE ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Pruning - all indices must be in the selection.
!---------------------------------------------------------------------------------------------------------------------------------*/
ChargeConstraint *ChargeConstraint_Prune ( const ChargeConstraint *self, const Selection *selection, Status *status )
{
    ChargeConstraint *new = NULL ;
    if ( self != NULL )
    {
        auto Boolean isPrunable = True ;
        auto Integer t ;
        /* . Check to see if all indices are in the selection. */
        for ( t = 0 ; t < self->numberOfCharges ; t++ ) isPrunable = isPrunable && selection->flags[ Integer1DArray_Item ( self->chargeIndices, t ) ] ;
        for ( t = 0 ; t < self->numberOfSpins   ; t++ ) isPrunable = isPrunable && selection->flags[ Integer1DArray_Item ( self->spinIndices  , t ) ] ;
        /* . Clone the constraint and change the indices. */
        if ( isPrunable )
        {
            new = ChargeConstraint_Clone ( self, status ) ;
            if ( new != NULL )
            {
                for ( t = 0 ; t < new->numberOfCharges ; t++ ) Integer1DArray_Item ( new->chargeIndices, t ) = selection->positions[ Integer1DArray_Item ( self->chargeIndices, t ) ] ;
                for ( t = 0 ; t < new->numberOfSpins   ; t++ ) Integer1DArray_Item ( new->spinIndices  , t ) = selection->positions[ Integer1DArray_Item ( self->spinIndices  , t ) ] ;
            }
            else Status_Set ( status, Status_MemoryAllocationFailure ) ;
        }
    }
    return new ;
}
