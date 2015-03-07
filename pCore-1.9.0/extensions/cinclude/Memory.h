/*------------------------------------------------------------------------------
! . File      : Memory.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . Memory allocation, deallocation and reallocation.
!=================================================================================================================================*/
# ifndef _MEMORY
# define _MEMORY

# include <stdlib.h>
# include "Definitions.h"

/* . Should put Integer, Real, etc. procedures into relevant modules. E.g. Real_Allocate and Real_Deallocate. Also Real_Initialize. */

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Allocate an object. */
# define MEMORY_ALLOCATE( object, type ) { object = ( type * ) malloc ( sizeof ( type ) ) ; }

/* . Allocate an array. */
# define MEMORY_ALLOCATEARRAY( object, length, type ) { object = ( type * ) calloc ( ( CSize ) length, sizeof ( type ) ) ; }

/* . Deallocate an array or object. */
# define MEMORY_DEALLOCATE( object ) { free ( ( void * ) object ) ; object = NULL ; }

/* . Reallocate an array. */
# define MEMORY_REALLOCATEARRAY( new, object, length, type ) { new = ( type * ) realloc ( ( void * ) object, ( ( CSize ) length ) * sizeof ( type ) ) ; }

/* . Obsolete. */

/* . Reallocate an array. */
# define MEMORY_REALLOCATEARRAYUNSAFE( object, length, type ) { object = ( type * ) realloc ( ( void * ) object, ( ( CSize ) length ) * sizeof ( type ) ) ; }

/* . Memory deallocation. */
# define Memory_Deallocate( argument ) { free ( argument ) ; argument = NULL ; }

# define Memory_AllocateArray Memory_Allocate_Array

/* . Real procedure aliases. */
# define Real_Allocate   Memory_Allocate_Array_Real
# define Real_Deallocate Memory_Deallocate_Real

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void      *Memory_Allocate                            ( const CSize nsize ) ;
extern void      *Memory_Allocate_Array                      ( const CSize nelements, const CSize     nsize ) ;
extern Boolean   *Memory_Allocate_Array_Boolean              ( const CSize nelements ) ;
extern Boolean   *Memory_Allocate_Array_Boolean_Initialize   ( const CSize nelements, const Boolean   value ) ;
extern Cardinal  *Memory_Allocate_Array_Cardinal             ( const CSize nelements ) ;
extern Cardinal  *Memory_Allocate_Array_Cardinal_Initialize  ( const CSize nelements, const Cardinal  value ) ;
extern Integer16 *Memory_Allocate_Array_Integer16            ( const CSize nelements ) ;
extern Integer16 *Memory_Allocate_Array_Integer16_Initialize ( const CSize nelements, const Integer16 value ) ;
extern Integer32 *Memory_Allocate_Array_Integer32            ( const CSize nelements ) ;
extern Integer32 *Memory_Allocate_Array_Integer32_Initialize ( const CSize nelements, const Integer32 value ) ;
extern Integer   *Memory_Allocate_Array_Integer              ( const CSize nelements ) ;
extern Integer   *Memory_Allocate_Array_Integer_Initialize   ( const CSize nelements, const Integer   value ) ;
extern Real      *Memory_Allocate_Array_Real                 ( const CSize nelements ) ;
extern Real      *Memory_Allocate_Array_Real_Initialize      ( const CSize nelements, const Real      value ) ;
extern void       Memory_Deallocate_Boolean                  ( Boolean **address ) ;
extern void       Memory_Deallocate_Integer                  ( Integer **address ) ;
extern void       Memory_Deallocate_Real                     ( Real    **address ) ;

# endif
