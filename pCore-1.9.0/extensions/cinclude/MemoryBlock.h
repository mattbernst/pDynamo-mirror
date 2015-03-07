/*------------------------------------------------------------------------------
! . File      : MemoryBlock.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _MEMORYBLOCK
# define _MEMORYBLOCK

# include <stdlib.h>
# include <string.h>

# include "Definitions.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The memory block type. */
typedef struct {
    CSize  capacity ;
    CSize  itemSize ;
    void  *data     ;
} MemoryBlock ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Objects. */
# define Memory_AllocateObject(   object, type )               { object = ( type * ) malloc ( sizeof ( type ) ) ; }
# define Memory_DeallocateObject( object       )               { free ( ( void * ) object ) ; object = NULL ; }

/* . Object arrays. */
# define Memory_AllocateArray(    object, extent, type )       { object = ( type * ) calloc ( ( CSize ) extent, sizeof ( type ) ) ; }
# define Memory_DeallocateArray                                Memory_DeallocateObject
# define Memory_ReallocateArray(  new, object, extent, type )  { new = ( type * ) realloc ( ( void * ) object, ( ( CSize ) extent ) * sizeof ( type ) ) ; }

/* . Raw arrays. */
# define Memory_AllocateRawArray(   object, extent, itemSize ) { object = calloc ( ( CSize ) extent, ( CSize ) itemSize ) ; }
# define Memory_DeallocateRawArray                             Memory_DeallocateObject

/* . Copying. */
# define Memory_CopyTo( self, other, size ) memcpy ( other, self, size )

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern MemoryBlock *MemoryBlock_Allocate   ( const CSize capacity, const CSize itemSize, Status *status ) ;
extern CSize        MemoryBlock_ByteSize   ( const MemoryBlock  *self ) ;
extern MemoryBlock *MemoryBlock_Clone      ( const MemoryBlock  *self, Status *status ) ;
extern void         MemoryBlock_Deallocate (       MemoryBlock **self ) ;

# endif
