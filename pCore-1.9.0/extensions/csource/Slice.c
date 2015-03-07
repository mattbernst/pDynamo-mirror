/*------------------------------------------------------------------------------
! . File      : Slice.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . Slice handling.
!=================================================================================================================================*/

# include "MemoryBlock.h"
# include "Slice.h"

/*==================================================================================================================================
! . Individual slices.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
Slice *Slice_Allocate ( Status *status )
{
    Slice *self = NULL ;
    Memory_AllocateObject ( self, Slice ) ;
    if ( self == NULL ) Status_Set ( status, Status_MemoryAllocationFailure ) ;
    else Slice_Initialize ( self ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Apply an extent to the slice.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Start and stop indices are adjusted by + extent if they are negative, except in the case of stop when the stride is negative. */
Boolean Slice_ApplyExtent ( Slice *self, const Integer extent, Status *status )
{
    Boolean isOK = False ;
    if ( self != NULL )
    {
        /* . Start position. */
        if ( self->start < 0 ) self->start += extent ;
        if ( ( self->start >= 0 ) && ( self->start < extent ) )
        {
            /* . Scalar. */
            if ( self->isScalar )
            {
                isOK         = True            ;
                self->extent = 1               ;
                self->stop   = self->start + 1 ;
                self->stride = 1               ;
            }
            /* . Slice. */
            else if ( self->stride != 0 )
            {
                /* . Stop position. */
                if ( ( self->stop  < 0 ) && ( self->stride > 0 ) ) self->stop += extent ;
                /* . Force start and stop to their appropriate limits. */
                                          if ( self->start <  0 ) self->start =  0 ; else if ( self->start > ( extent - 1 ) ) self->start = ( extent - 1 ) ;
                if ( self->stride > 0 ) { if ( self->stop  <  0 ) self->stop  =  0 ; else if ( self->stop  >   extent       ) self->stop  =   extent       ; }
                else                    { if ( self->stop  < -1 ) self->stop  = -1 ; else if ( self->stop  > ( extent - 1 ) ) self->stop  = ( extent - 1 ) ; }
                /* . Find the extent of the slice. */
                     if ( ( self->stride > 0 ) && ( self->stop > self->start ) ) self->extent = ( self->stop - self->start - 1 ) / self->stride + 1 ;
                else if ( ( self->stride < 0 ) && ( self->stop < self->start ) ) self->extent = ( self->stop - self->start + 1 ) / self->stride + 1 ;
                else self->extent = 0 ;
                isOK = True ;
            }
        }
        if ( ! isOK ) Status_Set ( status, Status_ArrayInvalidSlice ) ;
    }
    return isOK ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Slice_Deallocate ( Slice **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) ) Memory_DeallocateObject ( (*self) ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialization.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Slice_Initialize ( Slice *self )
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
void Slice_SetFromScalar ( Slice *self, Integer index )
{
    if ( self != NULL ) { self->isScalar = True ; self->start = index ; }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set a slice from a set of slice indices (with defaults).
!---------------------------------------------------------------------------------------------------------------------------------*/
void Slice_SetFromSliceIndices ( Slice *self, Integer *start, Integer *stop, Integer *stride, Status *status )
{
    if ( self != NULL )
    {
        self->isScalar = False ;
        if ( stride == NULL ) self->stride = 1 ;
        else                  self->stride = (*stride) ;
        if ( self->stride == 0 ) Status_Set ( status, Status_ArrayInvalidSlice ) ;
        else
        {
            auto Integer d = 0 ;
            if ( self->stride < 0 ) d = -1 ;
            if ( start == NULL ) self->start =  d ; else self->start = (*start) ;
            if ( stop  == NULL ) self->stop  = -1 ; else self->stop  = (*stop ) ;
        }
    }
}

/*==================================================================================================================================
! . Multislices.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
MultiSlice *MultiSlice_Allocate ( const Integer size, Status *status )
{
    MultiSlice *self = NULL ;
    if ( size >= 0 )
    {
        self = MultiSlice_AllocateRaw ( status ) ;
        if ( ( self != NULL ) && ( size > 0 ) )
        {
            Memory_AllocateArray ( self->items, size, Slice ) ;
            if ( self->items == NULL )
            {
                MultiSlice_Deallocate ( &self ) ;
                Status_Set ( status, Status_MemoryAllocationFailure ) ;
            }
            else
            {
                self->size = size ;
                auto Integer i ;
                for ( i = 0 ; i < size ; i++ ) Slice_Initialize ( &(self->items[i]) ) ;
            }
        }
    }
    else Status_Set ( status, Status_ArrayInvalidExtent ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Raw allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
MultiSlice *MultiSlice_AllocateRaw ( Status *status )
{
    MultiSlice *self = NULL ;
    Memory_AllocateObject ( self, MultiSlice ) ;
    if ( self == NULL ) Status_Set ( status, Status_MemoryAllocationFailure ) ;
    else { self->items = NULL ; self->rank = 0 ; self->size = 0 ; }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MultiSlice_Deallocate ( MultiSlice **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) )
    {
        Memory_DeallocateArray  ( (*self)->items ) ;
        Memory_DeallocateObject ( (*self)        ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check for conformability with a shape.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean MultiSlice_IsConformable ( MultiSlice *self, const Integer rank, const Integer *extents, Status *status )
{
    Boolean isOK = False ;
    if ( self != NULL )
    {
        if ( rank == self->size )
        {
            auto Integer d ;
            isOK = True ;
            for ( d = 0 ; d < rank ; d++ ) { isOK = isOK && Slice_ApplyExtent ( &(self->items[d]), extents[d], NULL ) ; }
        }
        if ( ! isOK ) Status_Set ( status, Status_ArrayInvalidSlice ) ;
    }
    return isOK ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set the rank from the constituent slices.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MultiSlice_SetRank ( MultiSlice *self )
{
    if ( self != NULL )
    {
        auto Integer i, rank = 0 ;
        for ( i = 0 ; i < self->size ; i++ ) { if ( ! self->items[i].isScalar ) rank += 1 ; }
        self->rank = rank ;
    }
}
