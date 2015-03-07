/*------------------------------------------------------------------------------
! . File      : BlockStorage.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==============================================================================
! . The block storage module handles indexed real data stored in blocks.
!=============================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "BlockStorage.h"
# include "Memory.h"

/*------------------------------------------------------------------------------
! . Definitions.
!-----------------------------------------------------------------------------*/
# define BLOCKSTORAGE_DEFAULTSIZE 1024

/*------------------------------------------------------------------------------
! . Local procedures.
!-----------------------------------------------------------------------------*/
static Block *Block_Allocate   ( const Integer blocksize, const Integer nindices16, const Integer nindices32 ) ;
static void   Block_Deallocate ( void *vblock ) ;

/*==============================================================================
! . Public procedures.
!=============================================================================*/
/*------------------------------------------------------------------------------
! . Allocation.
!-----------------------------------------------------------------------------*/
BlockStorage *BlockStorage_Allocate ( void )
{
   BlockStorage *blockstorage ;
   blockstorage = ( BlockStorage * ) Memory_Allocate ( sizeof ( BlockStorage ) ) ;
   blockstorage->blocksize     = BLOCKSTORAGE_DEFAULTSIZE ;
   blockstorage->ndata         = 0 ;
   blockstorage->nindices16    = 0 ;
   blockstorage->nindices32    = 0 ;
   blockstorage->QUNDERFLOW    = False ;
   blockstorage->underflow     = 0.0e+00 ;
   blockstorage->blocks        = List_Allocate ( ) ;
   blockstorage->blocks->Element_Deallocate = Block_Deallocate ;
   return blockstorage ;
}

/*------------------------------------------------------------------------------
! . Add data to the block storage.
! . The block storage must already be allocated and indices must be present
! . if the data is indexed.
!-----------------------------------------------------------------------------*/
void BlockStorage_Data_Add ( BlockStorage *blockstorage, const Integer ndata, const Real *data, const Integer16 *indices16, const Integer32 *indices32 )
{
   if ( ( blockstorage != NULL ) && ( ndata > 0 ) && ( data != NULL ) )
   {
      auto Block *block ;
      auto Integer    i, j ;
      /* . Get the current block. */
      if ( blockstorage->blocks->last == NULL )
      {
         block = Block_Allocate ( blockstorage->blocksize, blockstorage->nindices16, blockstorage->nindices32 ) ;
         List_Element_Append ( blockstorage->blocks, ( void * ) block ) ;
      }
      else block = ( Block * ) blockstorage->blocks->last->node ;
      /* . Copy all data. */
      if ( ! ( blockstorage->QUNDERFLOW ) )
      {
         for ( i = 0 ; i < ndata ; i++ )
         {
            if ( block->ndata >= blockstorage->blocksize )
            {
               block = Block_Allocate ( blockstorage->blocksize, blockstorage->nindices16, blockstorage->nindices32 ) ;
               List_Element_Append ( blockstorage->blocks, ( void * ) block ) ;
            }
            block->data[block->ndata] = data[i] ;
            for ( j = 0 ; j < blockstorage->nindices16 ; j++ ) block->indices16[blockstorage->nindices16*block->ndata+j] = indices16[blockstorage->nindices16*i+j] ;
            for ( j = 0 ; j < blockstorage->nindices32 ; j++ ) block->indices32[blockstorage->nindices32*block->ndata+j] = indices32[blockstorage->nindices32*i+j] ;
            block->ndata++ ;
         }
         blockstorage->ndata += ndata ;
      }
      /* . Copy data of a certain size only. */
      else
      {
         for ( i = 0 ; i < ndata ; i++ )
         {
            if ( block->ndata >= blockstorage->blocksize )
            {
               block = Block_Allocate ( blockstorage->blocksize, blockstorage->nindices16, blockstorage->nindices32 ) ;
               List_Element_Append ( blockstorage->blocks, ( void * ) block ) ;
            }
            if ( fabs ( data[i] ) > blockstorage->underflow )
            {
               block->data[block->ndata] = data[i] ;
               for ( j = 0 ; j < blockstorage->nindices16 ; j++ ) block->indices16[blockstorage->nindices16*block->ndata+j] = indices16[blockstorage->nindices16*i+j] ;
               for ( j = 0 ; j < blockstorage->nindices32 ; j++ ) block->indices32[blockstorage->nindices32*block->ndata+j] = indices32[blockstorage->nindices32*i+j] ;
               block->ndata++ ;
               blockstorage->ndata++ ;
            }
         }
      }
   }
}

/*------------------------------------------------------------------------------
! . Deallocation.
!-----------------------------------------------------------------------------*/
void BlockStorage_Deallocate ( BlockStorage **blockstorage )
{
   if ( (*blockstorage) != NULL )
   {
      List_Deallocate ( &((*blockstorage)->blocks) ) ;
      Memory_Deallocate ( (*blockstorage) ) ;
   }
}

/*------------------------------------------------------------------------------
! . Iterate over the interactions in a list.
!-----------------------------------------------------------------------------*/
Block *BlockStorage_Iterate ( BlockStorage *blockstorage )
{
   if ( blockstorage == NULL ) return NULL ;
   else                        return ( Block * ) List_Iterate ( blockstorage->blocks ) ;
}

/*------------------------------------------------------------------------------
! . Size of the block storage in (decimal) GB.
!-----------------------------------------------------------------------------*/
# define TOGIGABYTES 1.0e+09 /* . Assumes sizeof returns bytes. */
Real BlockStorage_Size ( BlockStorage *self )
{
    Real size = 0.0e+00 ;
    if  ( self != NULL )
    {
        auto Block *block = NULL ;

        /* . Initialization. */
        size = sizeof ( BlockStorage ) ;

        /* . Blocks. */
        if ( self->blocks != NULL )
        {
            size += sizeof ( List ) ;
            List_Iterate_Initialize ( self->blocks ) ;
            while ( ( block = BlockStorage_Iterate ( self ) ) != NULL ) size += block->ndata * ( self->nindices16 * sizeof ( Integer16 ) + self->nindices32 * sizeof ( Integer32 ) + sizeof ( Real ) ) ;
        }
    }
    return ( size / TOGIGABYTES ) ;
}
# undef TOGIGABYTES

/*==============================================================================
! . Private procedures.
!=============================================================================*/
/*------------------------------------------------------------------------------
! . Allocate a block.
!-----------------------------------------------------------------------------*/
static Block *Block_Allocate ( const Integer blocksize, const Integer nindices16, const Integer nindices32 )
{
   Block *block ;
   block = ( Block * ) Memory_Allocate ( sizeof ( Block ) ) ;
   block->ndata = 0 ;
   block->data  = Memory_Allocate_Array_Real ( blocksize ) ;
   if ( nindices16 > 0 ) block->indices16 = Memory_Allocate_Array_Integer16 ( blocksize * nindices16 ) ;
   else                  block->indices16 = NULL ;
   if ( nindices32 > 0 ) block->indices32 = Memory_Allocate_Array_Integer32 ( blocksize * nindices32 ) ;
   else                  block->indices32 = NULL ;
   return block ;
}

/*------------------------------------------------------------------------------
! . Deallocate a block.
!-----------------------------------------------------------------------------*/
static void Block_Deallocate ( void *vblock )
{
   Block *block ;
   block = ( Block * ) vblock ;
   free ( block->data      ) ;
   free ( block->indices16 ) ;
   free ( block->indices32 ) ;
   free ( block ) ;
}

