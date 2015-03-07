/*------------------------------------------------------------------------------
! . File      : MemoryBlock.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . Memory allocation, deallocation and reallocation.
!=================================================================================================================================*/

# include <stdio.h>
# include <stdlib.h>

# include "MemoryBlock.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
MemoryBlock *MemoryBlock_Allocate ( const CSize capacity, const CSize itemSize, Status *status )
{
    MemoryBlock *self = NULL ;
    if ( ( capacity >= 0 ) && ( itemSize > 0 ) )
    {
        Memory_AllocateObject ( self, MemoryBlock ) ;
        if ( self != NULL )
        {
            self->capacity = capacity ;
            self->itemSize = itemSize ;
            Memory_AllocateRawArray ( self->data, capacity, itemSize ) ;
            if ( self->data == NULL ) MemoryBlock_Deallocate ( &self ) ;
        }
    }
    if ( self == NULL ) Status_Set ( status, Status_MemoryAllocationFailure ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Byte size.
!---------------------------------------------------------------------------------------------------------------------------------*/
CSize MemoryBlock_ByteSize ( const MemoryBlock *self )
{
    CSize size = 0 ;
    if ( self != NULL ) size = self->capacity * self->itemSize ;
    return size ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
MemoryBlock *MemoryBlock_Clone ( const MemoryBlock *self, Status *status )
{
    MemoryBlock *clone = NULL ;
    if ( self != NULL )
    {
        clone = MemoryBlock_Allocate ( self->capacity, self->itemSize, status ) ;
        if ( clone != NULL ) Memory_CopyTo ( self->data, clone->data, MemoryBlock_ByteSize ( self ) ) ;
    }
    return clone ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MemoryBlock_Deallocate ( MemoryBlock **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) )
    {
        Memory_DeallocateRawArray ( (*self)->data ) ;
        Memory_DeallocateObject   ( (*self) ) ;
    }
}
