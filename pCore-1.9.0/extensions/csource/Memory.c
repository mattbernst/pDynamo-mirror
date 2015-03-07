/*------------------------------------------------------------------------------
! . File      : Memory.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==============================================================================
!=============================================================================*/

# include <stdio.h>
# include <stdlib.h>

# include "Memory.h"

/*------------------------------------------------------------------------------
! . Allocation.
!-----------------------------------------------------------------------------*/
void *Memory_Allocate ( const CSize nsize )
{
   return malloc ( nsize ) ;
}

/*------------------------------------------------------------------------------
! . Array allocation.
!-----------------------------------------------------------------------------*/
/* . General array. */
void *Memory_Allocate_Array ( const CSize nelements, const CSize nsize )
{
   return calloc ( nelements, nsize ) ;
}

/* . Booleans. */
Boolean *Memory_Allocate_Array_Boolean ( const CSize nelements )
{
   Boolean *array ;
   array = ( Boolean * ) Memory_Allocate_Array ( nelements, sizeof ( Boolean ) ) ;
   return array ;
}

Boolean *Memory_Allocate_Array_Boolean_Initialize ( const CSize nelements, const Boolean value )
{
   Boolean   *array ;
   Integer i ;
   array = Memory_Allocate_Array_Boolean ( nelements ) ;
   for ( i = 0 ; i < nelements ; i++ ) array[i] = value ;
   return array ;
}

/* . Cardinals. */
Cardinal *Memory_Allocate_Array_Cardinal ( const CSize nelements )
{
   Cardinal *array ;
   array = ( Cardinal * ) Memory_Allocate_Array ( nelements, sizeof ( Cardinal ) ) ;
   return array ;
}

Cardinal *Memory_Allocate_Array_Cardinal_Initialize ( const CSize nelements, const Cardinal value )
{
   Cardinal *array ;
   Cardinal i ;
   array = Memory_Allocate_Array_Cardinal ( nelements ) ;
   for ( i = 0 ; i < nelements ; i++ ) array[i] = value ;
   return array ;
}

/* . Integer16 and Integer32. */
Integer16 *Memory_Allocate_Array_Integer16 ( const CSize nelements )
{
   Integer16 *array ;
   array = ( Integer16 * ) Memory_Allocate_Array ( nelements, sizeof ( Integer16 ) ) ;
   return array ;
}

Integer16 *Memory_Allocate_Array_Integer16_Initialize ( const CSize nelements, const Integer16 value )
{
   Integer i ;
   Integer16 *array ;
   array = Memory_Allocate_Array_Integer16 ( nelements ) ;
   for ( i = 0 ; i < nelements ; i++ ) array[i] = value ;
   return array ;
}

Integer32 *Memory_Allocate_Array_Integer32 ( const CSize nelements )
{
   Integer32 *array ;
   array = ( Integer32 * ) Memory_Allocate_Array ( nelements, sizeof ( Integer32 ) ) ;
   return array ;
}

Integer32 *Memory_Allocate_Array_Integer32_Initialize ( const CSize nelements, const Integer32 value )
{
   Integer i ;
   Integer32 *array ;
   array = Memory_Allocate_Array_Integer32 ( nelements ) ;
   for ( i = 0 ; i < nelements ; i++ ) array[i] = value ;
   return array ;
}

/* . Integers. */
int *Memory_Allocate_Array_Integer ( const CSize nelements )
{
   Integer *array ;
   array = ( Integer * ) Memory_Allocate_Array ( nelements, sizeof ( Integer ) ) ;
   return array ;
}

int *Memory_Allocate_Array_Integer_Initialize ( const CSize nelements, const Integer value )
{
   Integer *array ;
   Integer i ;
   array = Memory_Allocate_Array_Integer ( nelements ) ;
   for ( i = 0 ; i < nelements ; i++ ) array[i] = value ;
   return array ;
}

/* . Real. */
Real *Memory_Allocate_Array_Real ( const CSize nelements )
{
   Real *array ;
   array = ( Real * ) Memory_Allocate_Array ( nelements, sizeof ( Real ) ) ;
   return array ;
}

Real *Memory_Allocate_Array_Real_Initialize ( const CSize nelements, const Real value )
{
   Integer i ;
   Real *array ;
   array = Memory_Allocate_Array_Real ( nelements ) ;
   for ( i = 0 ; i < nelements ; i++ ) array[i] = value ;
   return array ;
}

/*------------------------------------------------------------------------------
! . Deallocation.
!-----------------------------------------------------------------------------*/
void Memory_Deallocate_Boolean ( Boolean **address )
{
   if ( (*address) != NULL )
   {
      free ( (*address) ) ;
      (*address) = NULL ;
   }
}

void Memory_Deallocate_Real ( Real **address )
{
   if ( (*address) != NULL )
   {
      free ( (*address) ) ;
      (*address) = NULL ;
   }
}

void Memory_Deallocate_Integer ( Integer **address )
{
   if ( (*address) != NULL )
   {
      free ( (*address) ) ;
      (*address) = NULL ;
   }
}
