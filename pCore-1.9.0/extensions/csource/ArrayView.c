/*------------------------------------------------------------------------------
! . File      : ArrayView.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . Array views and associated objects.
!=================================================================================================================================*/

# include "Macros.h"
# include "Memory.h"
# include "ArrayView.h"

/*==================================================================================================================================
! . Array views.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation with rank.
!---------------------------------------------------------------------------------------------------------------------------------*/
ArrayView *ArrayView_Allocate ( const Integer rank, Status *status )
{
    ArrayView *self = NULL ;
    if ( rank >= 0 )
    {
        self = ArrayView_AllocateRaw ( status ) ;
        if ( ( self != NULL ) && ( rank > 0 ) )
        {
            self->rank = rank ;
            MEMORY_ALLOCATEARRAY ( self->extents, rank, Integer ) ;
            MEMORY_ALLOCATEARRAY ( self->strides, rank, Integer ) ;
            if ( ( self->extents == NULL ) || ( self->strides == NULL ) ) ArrayView_Deallocate ( &self ) ;
        }
        if ( self == NULL ) Status_Set ( status, Status_MemoryAllocationFailure ) ;
    }
    else Status_Set ( status, Status_InvalidRank ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Raw allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
ArrayView *ArrayView_AllocateRaw ( Status *status )
{
    ArrayView *self = NULL ;
    MEMORY_ALLOCATE ( self, ArrayView ) ;
    if ( self != NULL )
    {
        self->extents = NULL  ;
        self->strides = NULL  ;
        self->rank    = 0 ;
        self->offset  = 0 ;
        self->size    = 0 ;
    }
    else Status_Set ( status, Status_MemoryAllocationFailure ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation with extents.
!---------------------------------------------------------------------------------------------------------------------------------*/
ArrayView *ArrayView_AllocateWithExtents ( const Integer rank, const Integer *extents, Status *status )
{
    ArrayView *self = NULL ;
    self = ArrayView_Allocate ( rank, status ) ;
    /* . Determine extents and strides. */
    if ( self != NULL )
    {
        auto Integer d, extent, n ;
        for ( d = rank - 1, n = 1 ; d >= 0 ; d-- )
        {
            extent = extents[d] ;
            if ( extent < 0 )
            {
                extent = 0 ;
                Status_Set ( status, Status_InvalidArrayExtent ) ;
            }
            self->extents[d] = extent ;
            self->strides[d] = n ;
            n *= extent ;
        }
        self->size = n ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check the view versus a capacity.
! . Both maximum and minimum indices are calculated and checked.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean ArrayView_CheckCapacity ( const ArrayView *self, const Integer capacity )
{
    Boolean isOK = False ;
    if ( self != NULL )
    {
        auto Integer d, m = 0, n = 0, p ;
        if ( self->size > 0 )
        {
            m = self->offset ;
            n = self->offset ;
            for ( d = 0 ; d < self->rank ; d++ )
            {
                p = ( self->extents[d] - 1 ) * self->strides[d] ;
                if ( p < 0 ) m += p ;
                else         n += p ;
            }
        }
        isOK = ( ( m >= 0 ) && ( n < capacity ) ) ;
/*printf ( "Check Capacity> %d %d %d", capacity, m, n ) ;*/
    }
    return isOK ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
ArrayView *ArrayView_Clone ( const ArrayView *self, Status *status )
{
    ArrayView *clone = NULL ;
    if ( self != NULL )
    {
        clone = ArrayView_AllocateWithExtents ( self->rank, self->extents, status ) ;
        if ( clone != NULL )
        {
            auto Integer d ;
            clone->offset = self->offset ;
            for ( d = 0 ; d < self->rank ; d++ ) clone->strides[d] = self->strides[d] ;
        }
    }
    return clone ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void ArrayView_Deallocate ( ArrayView **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) )
    {
        MEMORY_DEALLOCATE ( (*self)->extents ) ;
        MEMORY_DEALLOCATE ( (*self)->strides ) ;
        MEMORY_DEALLOCATE ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Extent along a dimension.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer ArrayView_Extent ( const ArrayView *self, const Integer dimension, Status *status )
{
    Integer extent = 0 ;
    if ( self != NULL )
    {
        if ( ( dimension >= 0 ) && ( dimension < self->rank ) ) extent = self->extents[dimension] ; 
        else Status_Set ( status, Status_InvalidDimension ) ;
    }
    return extent ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check for compactness (minimal spacing of items).
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean ArrayView_IsCompact ( const ArrayView *self )
{
    Boolean isCompact = False ;
    if ( ( self != NULL ) && ( self->rank > 0 ) ) isCompact = ArrayView_IsUniform ( self ) && ( self->strides[self->rank-1] == 1 ) ;
    return isCompact ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check two views for conformability.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean ArrayView_IsConformable ( const ArrayView *self, const ArrayView *other )
{
    Boolean isConformable = False ;
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        if ( ( self->rank == other->rank ) && ( self->size == other->size ) )
        {
            auto Integer d ;
            isConformable = True ;
            for ( d = 0 ; d < self->rank ; d++ )
            {
                if ( self->extents[d] != other->extents[d] ) { isConformable = False ; break ; }
            }
        }
    }
    return isConformable ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check for uniformity (equal spacing of items).
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean ArrayView_IsUniform ( const ArrayView *self )
{
    Boolean isUniform = False ;
    if ( self != NULL )
    {
        auto Integer d ;
        isUniform = True ;
        for ( d = self->rank - 1 ; d >= 1 ; d-- )
        {
            if ( self->strides[d-1] != self->extents[d] * self->strides[d] ) { isUniform = False ; break ; }
        }
    }
    return isUniform ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get the index of an item.
! . The input rank can be less than the view rank in which case zero is assumed for the remaining indices.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer ArrayView_ItemIndex ( const ArrayView *self, const Integer rank, const Integer *indices )
{
    Integer index = -1 ;
    if ( self != NULL )
    {
        auto Integer d ;
        index = self->offset ;
        for ( d = 0 ; d < Minimum ( rank, self->rank ) ; d++ )
        {
            if ( ( indices[d] < 0 ) || ( indices[d] >= self->extents[d] ) ) { index = -1 ; break ; }
            index += indices[d] * self->strides[d] ;
        }
    }
    return index ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Rank.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer ArrayView_Rank ( const ArrayView *self )
{
    Integer rank = 0 ;
    if ( self != NULL ) rank = self->rank ;
    return rank ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Reshape to 1-D view.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean ArrayView_Reshape1D ( const ArrayView *self, Integer *offset, Integer *extent, Integer *stride, Status *status )
{
    Boolean isOK = False ;
    if ( ( self != NULL ) && ArrayView_IsUniform ( self ) )
    {
        (*extent) = self->size   ;
        (*offset) = self->offset ;
        (*stride) = self->strides[self->rank-1] ; /* . OK if row major. */
        isOK = True ;
    }
    else Status_Set ( status, Status_NonReshapeableArray ) ;
    return isOK ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Size.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer ArrayView_Size ( const ArrayView *self )
{
    Integer size = 0 ;
    if ( self != NULL ) size = self->size ;
    return size ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Slice a view.
!---------------------------------------------------------------------------------------------------------------------------------*/
Status ArrayView_Slice ( const ArrayView *self, const MultiSliceX *multiSlice, ArrayView **view )
{
    Status status = Status_Null ;
    if ( ( self != NULL ) && ( multiSlice != NULL ) && ( view != NULL ) && ( (*view) == NULL ) )
    {
        if ( self->rank != multiSlice->size ) status = Status_InvalidRank ;
        else
        {
            ArrayView *new = NULL ;
            new = ArrayView_Allocate ( self->rank, &status ) ;
            if ( new != NULL )
            {
                auto Integer d ;
                new->offset = self->offset ;
                new->rank   = 0 ;
                new->size   = 1 ;
                for ( d = 0 ; d < self->rank ; d++ )
                {
                    new->offset += ( multiSlice->items[d].start * self->strides[d] ) ;
                    if ( ! multiSlice->items[d].isScalar )
                    {
                        new->extents[new->rank] = multiSlice->items[d].extent ;
                        new->strides[new->rank] = multiSlice->items[d].stride * self->strides[d] ;
                        new->rank              += 1 ;
                        new->size              *= multiSlice->items[d].extent ;
                    }
                }
                status  = Status_Success ;
                (*view) = new ;
            }
        }
    }
    return status ;
}

/*==================================================================================================================================
! . Array view item iterators.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
ArrayViewItemIterator *ArrayViewItemIterator_Allocate ( ArrayView *target, Status *status )
{
    ArrayViewItemIterator *self = NULL ;
    if ( target != NULL )
    {
        MEMORY_ALLOCATE ( self, ArrayViewItemIterator ) ;
        if ( self != NULL )
        {
            self->current = -1     ;
            self->indices = NULL   ;
            self->target  = target ;
            if ( target->rank > 0 )
            {
                MEMORY_ALLOCATEARRAY ( self->indices, target->rank, Integer ) ;
                if ( self->indices == NULL ) ArrayViewItemIterator_Deallocate ( &self ) ;
            }
        }
        if ( self == NULL ) Status_Set ( status, Status_MemoryAllocationFailure ) ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void ArrayViewItemIterator_Deallocate ( ArrayViewItemIterator **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) )
    {
        MEMORY_DEALLOCATE ( (*self)->indices ) ;
        MEMORY_DEALLOCATE ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialization.
!---------------------------------------------------------------------------------------------------------------------------------*/
void ArrayViewItemIterator_Initialize ( ArrayViewItemIterator *self )
{
    if ( self != NULL )
    {
        auto Integer d ;
        self->current = self->target->offset ;
        for ( d = 0 ; d < self->target->rank ; d++ ) self->indices[d] = 0 ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Iteration.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer ArrayViewItemIterator_Next ( ArrayViewItemIterator *self )
{
    Integer current = -1 ;
    if ( ( self != NULL ) && ( self->current >= 0 ) )
    {
        Integer d, index ;
        current = self->current ;
        for ( d = self->target->rank-1 ; d >= 0 ; d-- )
        {
            index = self->indices[d] + 1 ;
            if ( index >= self->target->extents[d] )
            {
                if ( d == 0 ) self->current = -1 ;
                else
                {
                    self->indices[d] = 0 ;
                    self->current -= ( self->target->extents[d] - 1 ) * self->target->strides[d] ;
                }
            }
            else
            {
                self->indices[d] = index ;
                self->current += self->target->strides[d] ;
                break ;
            }
        }
    }
    return current ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Iteration - returning indices too.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer ArrayViewItemIterator_NextWithIndices ( ArrayViewItemIterator *self, Integer *indices )
{
    Integer current = -1 ;
    if ( ( self != NULL ) && ( self->current >= 0 ) )
    {
        Integer d, index ;
        current = self->current ;
        if ( indices != NULL )
        {
            for ( d = 0 ; d < self->target->rank ; d++ ) indices[d] = self->indices[d] ;
        }
        for ( d = self->target->rank-1 ; d >= 0 ; d-- )
        {
            index = self->indices[d] + 1 ;
            if ( index >= self->target->extents[d] )
            {
                if ( d == 0 ) self->current = -1 ;
                else
                {
                    self->indices[d] = 0 ;
                    self->current -= ( self->target->extents[d] - 1 ) * self->target->strides[d] ;
                }
            }
            else
            {
                self->indices[d] = index ;
                self->current += self->target->strides[d] ;
                break ;
            }
        }
    }
    return current ;
}
