/*------------------------------------------------------------------------------
! . File      : Boolean1DArray.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _BOOLEAN1DARRAY
# define _BOOLEAN1DARRAY

# include "Boolean.h"
# include "Integer.h"
# include "Integer1DArray.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The array type. */
typedef struct {
    Boolean  isOwner ;
    Boolean  isView  ;
    Integer  length  ;
    Integer  offset  ;
    Integer  size    ;
    Integer  stride  ;
    Boolean *data    ;
} Boolean1DArray ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . A pointer to the start of the data. */
# define Boolean1DArray_Data( self ) ( &((self)->data[(self)->offset]) )

/* . An item. */
# define Boolean1DArray_Item( self, i ) ( (self)->data[(self)->offset+(i)*(self)->stride] )

/* . A pointer to an item. */
# define Boolean1DArray_ItemPointer( self, i ) ( &((self)->data[(self)->offset+(i)*(self)->stride]) )

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern Boolean1DArray *Boolean1DArray_Allocate   ( const Integer length, Status *status ) ;
extern Boolean1DArray *Boolean1DArray_Clone      ( const Boolean1DArray  *self, Status *status ) ;
extern void            Boolean1DArray_CopyTo     ( const Boolean1DArray  *self, Boolean1DArray *other, Status *status ) ;
extern void            Boolean1DArray_Deallocate (       Boolean1DArray **self ) ;
extern Integer         Boolean1DArray_Length     ( const Boolean1DArray  *self ) ;
extern void            Boolean1DArray_Resize     (       Boolean1DArray  *self, const Integer length, const Boolean *initializer, Status *status ) ;
extern void            Boolean1DArray_Reverse    (       Boolean1DArray  *self ) ;
extern void            Boolean1DArray_Set        (       Boolean1DArray  *self, Boolean value ) ;
extern void            Boolean1DArray_Sort       (       Boolean1DArray  *self ) ;
extern void            Boolean1DArray_SortIndex  ( const Boolean1DArray  *self, Integer1DArray *indices, Status *status ) ;
extern void            Boolean1DArray_SortUnique (       Boolean1DArray  *self ) ;

# endif
