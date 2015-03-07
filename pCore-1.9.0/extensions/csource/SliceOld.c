/*------------------------------------------------------------------------------
! . File      : SliceOld.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . Slice handling.
!=================================================================================================================================*/

# include "Memory.h"
# include "SliceOld.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Slice indices:
!
! . start  - must always be a valid index [0,extent-1].
! . stop   - within the range [-1,extent].
! . stride - any non-zero integer.
!
! . It is tried to mimic Python as much as possible.
! . Unfortunately this has a problem when stop is -1
! . and stride is also negative.
!
!---------------------------------------------------------------------------------------------------------------------------------*/

/*==================================================================================================================================
! . Individual slices.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
Status SliceX_Allocate ( SliceX **self )
{
    SliceX *new = NULL ;
    if ( ( self != NULL ) && ( (*self) == NULL ) )
    {
        MEMORY_ALLOCATE  ( new, SliceX ) ;
        SliceX_Initialize ( new ) ;
        (*self) = new ;
    }
    return ( new == NULL ) ? Status_MemoryAllocationFailure : Status_Success ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SliceX_Deallocate ( SliceX **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) )
    {
        MEMORY_DEALLOCATE ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialization.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SliceX_Initialize ( SliceX *self )
{
    if ( self != NULL )
    {
        self->isScalar = False ;
        self->extent   = 0 ;
        self->start    = 0 ;
        self->stop     = 0 ;
        self->stride   = 1 ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set a slice from a scalar.
!---------------------------------------------------------------------------------------------------------------------------------*/
Status SliceX_SetFromScalar ( SliceX *self, Integer index, const Integer rawExtent )
{
    Status status = Status_Null ;
    if ( self != NULL )
    {
        if ( index < 0 ) index += rawExtent ;
        if ( ( index >= 0 ) && ( index < rawExtent ) )
        {
            self->isScalar = True      ;
            self->extent   = 1         ;
            self->start    = index     ;
            self->stop     = index + 1 ;
            self->stride   = 1         ;
            status = Status_Success ;
        }
        else status = Status_IndexOutOfRange ;
    }
    return status ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set a slice from a slice.
!---------------------------------------------------------------------------------------------------------------------------------*/
Status SliceX_SetFromSlice ( SliceX *self, Integer start, Integer stop, const Integer stride, const Integer rawExtent )
{
    Status status = Status_Null ;
    if ( self != NULL )
    {
        auto Integer extent ;
        /* . Adjust start and stop. */
        if ( start < 0 ) start += rawExtent ;
        if ( stop  < 0 )
        {
            if ( ! ( ( stop == -1 ) && ( stride < 0 ) ) ) stop += rawExtent ;
        }
        /* . Find extent. */
        extent = SliceIndices_GetExtent ( start, stop, stride, rawExtent, &status ) ;
        if ( extent >= 0 )
        {
            self->isScalar = False  ;
            self->extent   = extent ;
            self->start    = start  ;
            self->stop     = stop   ;
            self->stride   = stride ;
            status = Status_Success ;
        }
    }
    return status ;
}

/*==================================================================================================================================
! . Multislices.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
Status MultiSliceX_Allocate ( MultiSliceX **self, const Integer size )
{
    Status status ;
    if ( size >= 0 )
    {
        MultiSliceX *new = NULL ;
        if ( ( self != NULL ) && ( (*self) == NULL ) )
        {
            if ( Status_IsOK ( MultiSliceX_AllocateRaw ( &new ) ) )
            {
                if ( size > 0 )
                {
                    MEMORY_ALLOCATEARRAY ( new->items, size, SliceX ) ;
                    if ( new->items == NULL ) MultiSliceX_Deallocate ( &new ) ;
                    else
                    {
                        new->size = size ;
                        auto Integer i ;
                        for ( i = 0 ; i < size ; i++ ) SliceX_Initialize ( &(new->items[i]) ) ;
                    }
                }
                (*self) = new ;
            }
        }
        status = ( new == NULL ) ? Status_MemoryAllocationFailure : Status_Success ;
    }
    else status = Status_InvalidArraySize ;
    return status ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Raw allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
Status MultiSliceX_AllocateRaw ( MultiSliceX **self )
{
    MultiSliceX *new = NULL ;
    if ( ( self != NULL ) && ( (*self) == NULL ) )
    {
        MEMORY_ALLOCATE ( new, MultiSliceX ) ;
        if ( new != NULL )
        {
            new->items = NULL ;
            new->rank  = 0    ;
            new->size  = 0    ;
            (*self) = new ;
        }
    }
    return ( new == NULL ) ? Status_MemoryAllocationFailure : Status_Success ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MultiSliceX_Deallocate ( MultiSliceX **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) )
    {
        MEMORY_DEALLOCATE ( (*self)->items ) ;
        MEMORY_DEALLOCATE ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set the rank from the constituent slices.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MultiSliceX_SetRank ( MultiSliceX *self )
{
    if ( self != NULL )
    {
        auto Integer i, rank = 0 ;
        for ( i = 0 ; i < self->size ; i++ )
        {
            if ( ! self->items[i].isScalar ) rank += 1 ;
        }
        self->rank = rank ;
    }
}

/*==================================================================================================================================
! . Utility functions.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Determine the extent of a slice given the relevant indices.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer SliceIndices_GetExtent ( const Integer start, const Integer stop, const Integer stride, const Integer extent, Status *status )
{
    Integer n = -1 ;
    if ( ( start >= 0 ) && ( start < extent ) && ( stop >= -1 ) && ( stop <= extent ) && ( stride != 0 ) )
    {
        if      ( ( stride > 0 ) && ( stop > start ) ) n = ( stop - start - 1 ) / stride + 1 ;
        else if ( ( stride < 0 ) && ( stop < start ) ) n = ( stop - start + 1 ) / stride + 1 ;
        else n = 0 ;
        Status_Set ( status, Status_Continue ) ;
    }
    else Status_Set ( status, Status_InvalidSlice ) ;
    return n ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Determine the length of a slice given the relevant indices.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer SliceIndices_GetLength ( const Integer start, const Integer stop, const Integer stride, const Integer length, Status *status )
{
    Integer n = 0 ;
    if ( ( start >= 0 ) && ( start < length ) && ( stop >= -1 ) && ( stop <= length ) && ( stride != 0 ) )
    {
        if      ( ( stride > 0 ) && ( stop > start ) ) n = ( stop - start - 1 ) / stride + 1 ;
        else if ( ( stride < 0 ) && ( stop < start ) ) n = ( stop - start + 1 ) / stride + 1 ;
        Status_Set ( status, Status_Continue ) ;
    }
    else Status_Set ( status, Status_InvalidSliceIndices ) ;
    return n ;
}
